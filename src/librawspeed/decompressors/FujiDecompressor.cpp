/*
    RawSpeed - RAW file decoder.

    Copyright (C) 2016 Alexey Danilchenko
    Copyright (C) 2016 Alex Tutubalin
    Copyright (C) 2017 Uwe Müssel
    Copyright (C) 2017 Roman Lebedev

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/

#include "rawspeedconfig.h"
#include "decompressors/FujiDecompressor.h"
#include "MemorySanitizer.h"
#include "adt/Array1DRef.h"
#include "adt/Array2DRef.h"
#include "adt/Casts.h"
#include "adt/CroppedArray1DRef.h"
#include "adt/CroppedArray2DRef.h"
#include "adt/Invariant.h"
#include "adt/Optional.h"
#include "adt/Point.h"
#include "bitstreams/BitStreamerMSB.h"
#include "common/BayerPhase.h"
#include "common/Common.h"
#include "common/RawImage.h"
#include "common/XTransPhase.h"
#include "decoders/RawDecoderException.h"
#include "io/Buffer.h"
#include "io/ByteStream.h"
#include "io/Endianness.h"
#include "metadata/ColorFilterArray.h"
#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace rawspeed {

namespace {

struct BayerTag;
struct XTransTag;

template <typename T> constexpr iPoint2D MCU;

template <> constexpr iPoint2D MCU<BayerTag> = {2, 2};

template <> constexpr iPoint2D MCU<XTransTag> = {6, 6};

struct int_pair final {
  int value1;
  int value2;
};

struct FujiQTable final {
  std::vector<int8_t> q_table;
  std::array<int, 5> q_point = {};

  int q_base = 0;
  int max_bits = 0;
  int min_value = 0x40;
  int raw_bits = 0;
  int total_values = 0;
  int maxDiff = 0;
  int max_grad = 0;
  int q_grad_mult = 9;

  [[nodiscard]] int8_t lookup(int delta) const {
    return q_table[maxValue() + delta];
  }

  [[nodiscard]] int maxValue() const { return q_point[4]; }

  [[nodiscard]] int initialMaxDiff() const { return maxDiff; }
};

enum xt_lines : uint8_t {
  R0 = 0,
  R1,
  R2,
  R3,
  R4,
  G0,
  G1,
  G2,
  G3,
  G4,
  G5,
  G6,
  G7,
  B0,
  B1,
  B2,
  B3,
  B4,
  ltotal
};

struct FujiGradients final {
  std::array<int_pair, 41> main;
  std::array<std::array<int_pair, 5>, 3> lossyStatic;
};

struct FujiPrediction final {
  int v1 = 0;
  int v2 = 0;
  int interp_val = 0;
  int edge_metric = 0;
};

int RAWSPEED_READNONE fujiLog2Ceil(int value) {
  invariant(value > 0);
  return implicit_cast<int>(std::bit_width(static_cast<unsigned>(value - 1)));
}

int RAWSPEED_READNONE maxValueForRawBits(int raw_bits) {
  invariant(raw_bits > 0);
  return (1 << raw_bits) - 1;
}

void populateQTable(FujiQTable& qt) {
  const int maxValue = qt.maxValue();
  qt.q_table.resize(2 * (maxValue + 1));

  for (int delta = -maxValue; delta <= maxValue; ++delta) {
    int8_t grad;
    if (delta <= -qt.q_point[3])
      grad = -4;
    else if (delta <= -qt.q_point[2])
      grad = -3;
    else if (delta <= -qt.q_point[1])
      grad = -2;
    else if (delta < -qt.q_point[0])
      grad = -1;
    else if (delta <= qt.q_point[0])
      grad = 0;
    else if (delta < qt.q_point[1])
      grad = 1;
    else if (delta < qt.q_point[2])
      grad = 2;
    else if (delta < qt.q_point[3])
      grad = 3;
    else
      grad = 4;

    qt.q_table[maxValue + delta] = grad;
  }
}

struct fuji_compressed_params final {
  explicit fuji_compressed_params(
      const FujiDecompressor::FujiHeader& h,
      Optional<Array1DRef<const Array1DRef<const uint8_t>>> qbaseStrips);

  FujiQTable qtable;
  struct FujiLossyQTableCache final {
    std::bitset<256> observed_qbases;
    std::array<int, 256> lossy_main_qtable_indices = {};
    std::vector<FujiQTable> lossy_main_qtables;
    std::array<FujiQTable, 3> lossy_static_qtables;

    FujiLossyQTableCache() { lossy_main_qtable_indices.fill(-1); }

    void init(const FujiDecompressor::FujiHeader& h,
              Array1DRef<const Array1DRef<const uint8_t>> qbaseStrips);

    [[nodiscard]] const FujiQTable& mainLossyQTable(uint8_t qbase) const {
      const int idx = lossy_main_qtable_indices[qbase];
      invariant(idx >= 0);
      invariant(idx < implicit_cast<int>(lossy_main_qtables.size()));
      return lossy_main_qtables[idx];
    }

    [[nodiscard]] const FujiQTable& staticLossyQTable(int idx) const {
      invariant(idx >= 0 &&
                idx < implicit_cast<int>(lossy_static_qtables.size()));
      return lossy_static_qtables[idx];
    }
  };

  FujiLossyQTableCache lossy_qtables;
  uint16_t line_width = 0;
};

struct FujiStrip final {
  // part of which 'image' this block is
  const FujiDecompressor::FujiHeader& h;

  // which strip is this, 0 .. h.blocks_in_row-1
  const int n;

  // the compressed data of this strip
  const Array1DRef<const uint8_t> input;

  // the active q_base entries for this strip, one per encoded line
  const Optional<Array1DRef<const uint8_t>> qbase;

  FujiStrip() = delete;
  FujiStrip(const FujiStrip&) = delete;
  FujiStrip(FujiStrip&&) noexcept = delete;
  FujiStrip& operator=(const FujiStrip&) noexcept = delete;
  FujiStrip& operator=(FujiStrip&&) noexcept = delete;

  FujiStrip(const FujiDecompressor::FujiHeader& h_, int block,
            Array1DRef<const uint8_t> input_,
            Optional<Array1DRef<const uint8_t>> qbase_)
      : h(h_), n(block), input(input_), qbase(qbase_) {
    invariant(n >= 0 && n < h.blocks_in_row);
    invariant(!h.isLossy() || qbase);
  }

  // each strip's line corresponds to 6 output lines.
  static int RAWSPEED_READONLY lineHeight() { return 6; }

  // how many vertical lines does this block encode?
  [[nodiscard]] int RAWSPEED_READONLY height() const { return h.total_lines; }

  // how many horizontal pixels does this block encode?
  [[nodiscard]] int RAWSPEED_READONLY width() const {
    // if this is not the last block, we are good.
    if ((n + 1) != h.blocks_in_row)
      return h.block_size;

    // ok, this is the last block...

    invariant(h.block_size * h.blocks_in_row >= h.raw_width);
    return h.raw_width - offsetX();
  }

  // how many horizontal pixels does this block encode?
  [[nodiscard]] iPoint2D numMCUs(iPoint2D MCU) const {
    invariant(width() % MCU.x == 0);
    invariant(lineHeight() % MCU.y == 0);
    return {width() / MCU.x, lineHeight() / MCU.y};
  }

  // where vertically does this block start?
  [[nodiscard]] int offsetY(int line = 0) const {
    (void)height(); // A note for NDEBUG builds that *this is used.
    invariant(line >= 0 && line < height());
    return lineHeight() * line;
  }

  // where horizontally does this block start?
  [[nodiscard]] int offsetX() const { return h.block_size * n; }

  [[nodiscard]] uint8_t qbaseForLine(int line) const {
    invariant(qbase);
    invariant(line >= 0 && line < height());
    return (*qbase)(line);
  }
};

FujiQTable makeLosslessQTable(int raw_bits) {
  FujiQTable qt;

  qt.q_point[0] = 0;
  qt.q_point[1] = 0x12;
  qt.q_point[2] = 0x43;
  qt.q_point[3] = 0x114;
  qt.q_point[4] = maxValueForRawBits(raw_bits);

  qt.raw_bits = raw_bits;

  populateQTable(qt);

  if (qt.maxValue() == 0xFFFF) { // (1 << raw_bits) - 1
    qt.total_values = 0x10000;   // 1 << raw_bits
    qt.raw_bits = 16;            // raw_bits
    qt.max_bits = 64;            // raw_bits * (64 / raw_bits)
    qt.maxDiff = 1024;           // 1 << (raw_bits - 6)
  } else if (qt.maxValue() == 0x3FFF) {
    qt.total_values = 0x4000;
    qt.raw_bits = 14;
    qt.max_bits = 56;
    qt.maxDiff = 256;
  } else if (qt.maxValue() == 0xFFF) {
    qt.total_values = 4096;
    qt.raw_bits = 12;
    qt.max_bits = 48; // out-of-family, there's greater pattern at play.
    qt.maxDiff = 64;

    ThrowRDE("Aha, finally, a 12-bit compressed RAF! Please consider providing "
             "samples on <https://raw.pixls.us/>, thanks!");
  } else {
    ThrowRDE("FUJI q_point");
  }

  return qt;
}

FujiQTable makeLossyMainQTable(int raw_bits, uint8_t qbase) {
  FujiQTable qt;

  const int qbaseInt = qbase;
  const int maxValue = maxValueForRawBits(raw_bits);
  const int maxValuePlusOne = maxValue + 1;

  qt.q_point[0] = qbaseInt;
  qt.q_point[1] = (3 * qbaseInt) + 0x12;
  qt.q_point[2] = (5 * qbaseInt) + 0x43;
  qt.q_point[3] = (7 * qbaseInt) + 0x114;
  qt.q_point[4] = maxValue;

  if (qt.q_point[1] >= maxValuePlusOne || qt.q_point[1] < qbaseInt + 1)
    qt.q_point[1] = qbaseInt + 1;
  if (qt.q_point[2] < qt.q_point[1] || qt.q_point[2] >= maxValuePlusOne)
    qt.q_point[2] = qt.q_point[1];
  if (qt.q_point[3] < qt.q_point[2] || qt.q_point[3] >= maxValuePlusOne)
    qt.q_point[3] = qt.q_point[2];

  qt.q_base = qbaseInt;
  qt.total_values = (maxValue + (2 * qbaseInt)) / ((2 * qbaseInt) + 1) + 1;
  qt.raw_bits = fujiLog2Ceil(qt.total_values);
  qt.max_bits = 4 * fujiLog2Ceil(maxValuePlusOne);
  qt.maxDiff = std::max(2, (qt.total_values + 0x20) >> 6);
  qt.q_grad_mult = 9;

  populateQTable(qt);
  return qt;
}

FujiQTable makeLossyStaticQTable(int raw_bits, int idx) {
  FujiQTable qt;

  const int maxValue = maxValueForRawBits(raw_bits);

  qt.q_point[4] = maxValue;
  qt.q_base = idx;
  qt.max_grad = 5 + idx;
  qt.q_grad_mult = 3;
  qt.max_bits = 4 * fujiLog2Ceil(maxValue + 1);

  switch (idx) {
  case 0:
    qt.q_point[0] = 0;
    qt.q_point[1] = maxValue >= 0x12 ? 0x12 : qt.q_point[0] + 1;
    qt.q_point[2] = maxValue >= 0x43 ? 0x43 : qt.q_point[1];
    qt.q_point[3] = maxValue >= 0x114 ? 0x114 : qt.q_point[2];
    qt.total_values = maxValue + 1;
    break;
  case 1:
    qt.q_point[0] = 1;
    qt.q_point[1] = maxValue >= 0x15 ? 0x15 : qt.q_point[0] + 1;
    qt.q_point[2] = maxValue >= 0x48 ? 0x48 : qt.q_point[1];
    qt.q_point[3] = maxValue >= 0x11B ? 0x11B : qt.q_point[2];
    qt.total_values = (maxValue + 2) / 3 + 1;
    break;
  case 2:
    qt.q_point[0] = 2;
    qt.q_point[1] = maxValue >= 0x18 ? 0x18 : qt.q_point[0] + 1;
    qt.q_point[2] = maxValue >= 0x4D ? 0x4D : qt.q_point[1];
    qt.q_point[3] = maxValue >= 0x122 ? 0x122 : qt.q_point[2];
    qt.total_values = (maxValue + 4) / 5 + 1;
    break;
  default:
    __builtin_unreachable();
  }

  qt.raw_bits = fujiLog2Ceil(qt.total_values);
  qt.maxDiff = std::max(2, (qt.total_values + 0x20) >> 6);

  populateQTable(qt);
  return qt;
}

std::array<FujiQTable, 3> makeLossyStaticQTables(int raw_bits) {
  return {makeLossyStaticQTable(raw_bits, 0),
          makeLossyStaticQTable(raw_bits, 1),
          makeLossyStaticQTable(raw_bits, 2)};
}

void fuji_compressed_params::FujiLossyQTableCache::init(
    const FujiDecompressor::FujiHeader& h,
    Array1DRef<const Array1DRef<const uint8_t>> qbaseStrips) {
  invariant(h.isLossy());
  invariant(qbaseStrips.size() == h.blocks_in_row);

  for (int strip = 0; strip != qbaseStrips.size(); ++strip) {
    const auto qbases = qbaseStrips(strip);
    invariant(qbases.size() == h.total_lines);
    for (int line = 0; line != qbases.size(); ++line)
      observed_qbases.set(qbases(line));
  }

  lossy_static_qtables = makeLossyStaticQTables(h.raw_bits);
  lossy_main_qtables.reserve(observed_qbases.count());

  for (int qbase = 0; qbase != 256; ++qbase) {
    if (!observed_qbases[qbase])
      continue;

    lossy_main_qtable_indices[qbase] =
        implicit_cast<int>(lossy_main_qtables.size());
    lossy_main_qtables.emplace_back(
        makeLossyMainQTable(h.raw_bits, implicit_cast<uint8_t>(qbase)));
  }
}

fuji_compressed_params::fuji_compressed_params(
    const FujiDecompressor::FujiHeader& h,
    Optional<Array1DRef<const Array1DRef<const uint8_t>>> qbaseStrips) {
  if ((h.block_size % 3 && h.raw_type == 16) ||
      (h.block_size & 1 && h.raw_type == 0)) {
    ThrowRDE("fuji_block_checks");
  }

  if (h.raw_type == 16) {
    line_width = (h.block_size * 2) / 3;
  } else {
    line_width = h.block_size >> 1;
  }

  if (h.isLossless()) {
    qtable = makeLosslessQTable(h.raw_bits);
  } else {
    if (!qbaseStrips)
      ThrowRDE("Fujifilm lossy RAF q_base value missing");
    lossy_qtables.init(h, *qbaseStrips);
  }
}

struct fuji_compressed_block final {
  const Array2DRef<uint16_t> img;
  const FujiDecompressor::FujiHeader& header;
  const fuji_compressed_params& common_info;

  fuji_compressed_block(Array2DRef<uint16_t> img,
                        const FujiDecompressor::FujiHeader& header,
                        const fuji_compressed_params& common_info);

  void reset();

  Optional<BitStreamerMSB> pump;

  // tables of gradients
  std::array<FujiGradients, 3> grad_even;
  std::array<FujiGradients, 3> grad_odd;

  const FujiQTable* current_lossy_main_qtable = nullptr;
  int current_lossy_qbase = -1;

  std::vector<uint16_t> linealloc;
  Array2DRef<uint16_t> lines;

  void fuji_decode_strip(const FujiStrip& strip);

  template <typename Tag, typename T>
  void copy_line(const FujiStrip& strip, int cur_line, T idx) const;

  void copy_line_to_xtrans(const FujiStrip& strip, int cur_line) const;
  void copy_line_to_bayer(const FujiStrip& strip, int cur_line) const;

  static int fuji_zerobits(BitStreamerMSB& pump);
  static int bitDiff(int value1, int value2);

  void resetMainGradients(const FujiQTable& qtable);
  void setLossyQBase(uint8_t qbase);

  template <std::size_t NumGradients>
  [[nodiscard]] int
  fuji_decode_sample(const FujiQTable& qtable, int grad, int interp_val,
                     std::array<int_pair, NumGradients>& grads, bool lossy);
  [[nodiscard]] int fuji_decode_sample_even(xt_lines c, int col,
                                            std::array<int_pair, 41>& grads);
  [[nodiscard]] int fuji_decode_sample_odd(xt_lines c, int col,
                                           std::array<int_pair, 41>& grads);
  [[nodiscard]] int fuji_decode_lossy_sample(FujiPrediction pred,
                                             FujiGradients& gradients);
  [[nodiscard]] int fuji_decode_lossy_sample_even(xt_lines c, int col,
                                                  FujiGradients& gradients);
  [[nodiscard]] int fuji_decode_lossy_sample_odd(xt_lines c, int col,
                                                 FujiGradients& gradients);

  [[nodiscard]] int fuji_quant_gradient(const FujiQTable& qtable, int v1,
                                        int v2) const;

  [[nodiscard]] FujiPrediction
  fuji_decode_interpolation_even_inner(xt_lines c, int col) const;
  [[nodiscard]] FujiPrediction
  fuji_decode_interpolation_odd_inner(xt_lines c, int col) const;
  [[nodiscard]] int fuji_decode_interpolation_even(xt_lines c, int col) const;

  void fuji_extend_generic(int start, int end) const;
  void fuji_extend_red() const;
  void fuji_extend_green() const;
  void fuji_extend_blue() const;

  template <typename TEven, typename TOdd>
  void fuji_decode_block(TEven func_even, TOdd func_odd, int cur_line);
  void xtrans_decode_block(int cur_line);
  void fuji_bayer_decode_block(int cur_line);
};

fuji_compressed_block::fuji_compressed_block(
    Array2DRef<uint16_t> img_, const FujiDecompressor::FujiHeader& header_,
    const fuji_compressed_params& common_info_)
    : img(img_), header(header_), common_info(common_info_),
      linealloc(ltotal * (common_info.line_width + 2), 0),
      lines(&linealloc[0], common_info.line_width + 2, ltotal) {}

void fuji_compressed_block::reset() {
  MSan::Allocated(CroppedArray2DRef(lines));
  current_lossy_main_qtable = nullptr;
  current_lossy_qbase = -1;

  // Zero-initialize first two (read-only, carry-in) lines of each color,
  // including first and last helper columns of the second row.
  // This is needed for correctness.
  for (xt_lines color : {R0, G0, B0}) {
    memset(&lines(color, 0), 0, 2 * sizeof(uint16_t) * lines.width());

    // On the first row, we don't need to zero-init helper columns.
    MSan::Allocated(lines(color, 0));
    MSan::Allocated(lines(color, lines.width() - 1));
  }

  // And the first (real, uninitialized) line of each color gets the content
  // of the last helper column from the last decoded sample of previous
  // line of that color.
  // Again, this is needed for correctness.
  for (xt_lines color : {R2, G2, B2})
    lines(color, lines.width() - 1) = lines(color - 1, lines.width() - 2);

  if (header.isLossless()) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 41; i++) {
        grad_even[j].main[i].value1 = common_info.qtable.initialMaxDiff();
        grad_even[j].main[i].value2 = 1;
        grad_odd[j].main[i].value1 = common_info.qtable.initialMaxDiff();
        grad_odd[j].main[i].value2 = 1;
      }
    }
  }

  if (header.isLossy()) {
    for (int table = 0; table != 3; ++table) {
      const FujiQTable& qtable =
          common_info.lossy_qtables.staticLossyQTable(table);
      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 5; ++i) {
          grad_even[j].lossyStatic[table][i].value1 = qtable.initialMaxDiff();
          grad_even[j].lossyStatic[table][i].value2 = 1;
          grad_odd[j].lossyStatic[table][i].value1 = qtable.initialMaxDiff();
          grad_odd[j].lossyStatic[table][i].value2 = 1;
        }
      }
    }
  }
}

void fuji_compressed_block::resetMainGradients(const FujiQTable& qtable) {
  for (int j = 0; j != 3; ++j) {
    for (int i = 0; i != 41; ++i) {
      grad_even[j].main[i].value1 = qtable.initialMaxDiff();
      grad_even[j].main[i].value2 = 1;
      grad_odd[j].main[i].value1 = qtable.initialMaxDiff();
      grad_odd[j].main[i].value2 = 1;
    }
  }
}

void fuji_compressed_block::setLossyQBase(uint8_t qbase) {
  if (current_lossy_qbase == implicit_cast<int>(qbase)) {
    invariant(current_lossy_main_qtable != nullptr);
    return;
  }

  current_lossy_qbase = qbase;
  current_lossy_main_qtable = &common_info.lossy_qtables.mainLossyQTable(qbase);
  resetMainGradients(*current_lossy_main_qtable);
}

template <typename Tag, typename T>
void fuji_compressed_block::copy_line(const FujiStrip& strip, int cur_line,
                                      T idx) const {
  std::array<CFAColor, MCU<Tag>.x * MCU<Tag>.y> CFAData;
  if constexpr (std::is_same_v<XTransTag, Tag>)
    CFAData = getAsCFAColors(XTransPhase(0, 0));
  else if constexpr (std::is_same_v<BayerTag, Tag>)
    CFAData = getAsCFAColors(BayerPhase::RGGB);
  else
    __builtin_unreachable();
  const Array2DRef<const CFAColor> CFA(CFAData.data(), MCU<Tag>.x, MCU<Tag>.y);

  iPoint2D MCUIdx;
  assert(MCU<Tag> == strip.h.MCU);
  const iPoint2D NumMCUs = strip.numMCUs(MCU<Tag>);
  for (MCUIdx.x = 0; MCUIdx.x != NumMCUs.x; ++MCUIdx.x) {
    for (MCUIdx.y = 0; MCUIdx.y != NumMCUs.y; ++MCUIdx.y) {
      const auto out =
          CroppedArray2DRef(img, strip.offsetX() + (MCU<Tag>.x * MCUIdx.x),
                            strip.offsetY(cur_line) + (MCU<Tag>.y * MCUIdx.y),
                            MCU<Tag>.x, MCU<Tag>.y);
      for (int MCURow = 0; MCURow != MCU<Tag>.y; ++MCURow) {
        for (int MCUCol = 0; MCUCol != MCU<Tag>.x; ++MCUCol) {
          int imgRow = (MCU<Tag>.y * MCUIdx.y) + MCURow;
          int imgCol = (MCU<Tag>.x * MCUIdx.x) + MCUCol;

          int row;

          switch (CFA(MCURow, MCUCol)) {
            using enum CFAColor;
          case RED: // red
            row = R2 + (imgRow >> 1);
            break;

          case GREEN: // green
            row = G2 + imgRow;
            break;

          case BLUE: // blue
            row = B2 + (imgRow >> 1);
            break;

          default:
            __builtin_unreachable();
          }

          out(MCURow, MCUCol) = lines(row, 1 + idx(imgCol));
        }
      }
    }
  }
}

void fuji_compressed_block::copy_line_to_xtrans(const FujiStrip& strip,
                                                int cur_line) const {
  auto index = [](int imgCol) {
    return (((imgCol * 2 / 3) & 0x7FFFFFFE) | ((imgCol % 3) & 1)) +
           ((imgCol % 3) >> 1);
  };

  copy_line<XTransTag>(strip, cur_line, index);
}

void fuji_compressed_block::copy_line_to_bayer(const FujiStrip& strip,
                                               int cur_line) const {
  auto index = [](int imgCol) { return imgCol >> 1; };

  copy_line<BayerTag>(strip, cur_line, index);
}

inline int fuji_compressed_block::fuji_zerobits(BitStreamerMSB& pump) {
  int count = 0;

  // Count-and-skip all the leading `0`s.
  while (true) {
    constexpr int batchSize = 32;
    pump.fill(batchSize);
    uint32_t batch = pump.peekBitsNoFill(batchSize);
    int numZerosInThisBatch = std::countl_zero(batch);
    count += numZerosInThisBatch;
    bool allZeroes = numZerosInThisBatch == batchSize;
    int numBitsToSkip = numZerosInThisBatch;
    if (!allZeroes)
      numBitsToSkip += 1; // Also skip the first `1`.
    pump.skipBitsNoFill(numBitsToSkip);
    if (!allZeroes)
      break; // We're done!
  }

  return count;
}

// Given two non-negative numbers, how many times must the second number
// be multiplied by 2, for it to become not smaller than the first number?
// We are operating on arithmetical numbers here, without overflows.
int RAWSPEED_READNONE fuji_compressed_block::bitDiff(int value1, int value2) {
  invariant(value1 >= 0);
  invariant(value2 > 0);

  int lz1 = std::countl_zero(static_cast<unsigned>(value1));
  int lz2 = std::countl_zero(static_cast<unsigned>(value2));
  int decBits = std::max(lz2 - lz1, 0);
  if ((value2 << decBits) < value1)
    ++decBits;
  return std::min(decBits, 15);
}

template <std::size_t NumGradients>
__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_sample(
    const FujiQTable& qtable, int grad, int interp_val,
    std::array<int_pair, NumGradients>& grads, bool lossy) {
  const int gradient = std::abs(grad);
  invariant(gradient < implicit_cast<int>(grads.size()));

  int sampleBits = fuji_zerobits(*pump);

  int codeBits;
  int codeDelta;
  if (sampleBits < qtable.max_bits - qtable.raw_bits - 1) {
    codeBits = bitDiff(grads[gradient].value1, grads[gradient].value2);
    codeDelta = sampleBits << codeBits;
  } else {
    codeBits = qtable.raw_bits;
    codeDelta = 1;
  }

  int code = 0;
  pump->fill(32);
  if (codeBits)
    code = pump->getBitsNoFill(codeBits);
  code += codeDelta;

  if (code < 0 || code >= qtable.total_values) {
    ThrowRDE("fuji_decode_sample");
  }

  if (code & 1) {
    code = -1 - code / 2;
  } else {
    code /= 2;
  }

  grads[gradient].value1 += std::abs(code);

  if (grads[gradient].value2 == qtable.min_value) {
    grads[gradient].value1 >>= 1;
    grads[gradient].value2 >>= 1;
  }

  grads[gradient].value2++;

  const int scale = lossy ? (2 * qtable.q_base) + 1 : 1;
  if (grad < 0) {
    interp_val -= code * scale;
  } else {
    interp_val += code * scale;
  }

  if (lossy) {
    if (interp_val < -qtable.q_base) {
      interp_val += qtable.total_values * scale;
    } else if (interp_val > qtable.q_base + qtable.maxValue()) {
      interp_val -= qtable.total_values * scale;
    }
  } else {
    if (interp_val < 0) {
      interp_val += qtable.total_values;
    } else if (interp_val > qtable.maxValue()) {
      interp_val -= qtable.total_values;
    }
  }

  if (interp_val < 0)
    return 0;

  return std::min(interp_val, qtable.maxValue());
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_quant_gradient(const FujiQTable& qtable, int v1,
                                           int v2) const {
  return (qtable.q_grad_mult * qtable.lookup(v1)) + qtable.lookup(v2);
}

__attribute__((always_inline)) inline FujiPrediction
fuji_compressed_block::fuji_decode_interpolation_even_inner(xt_lines c,
                                                            int col) const {
  int Rb = lines(c - 1, 1 + (2 * (col + 0)) + 0);
  int Rc = lines(c - 1, 1 + (2 * (col - 1)) + 1);
  int Rd = lines(c - 1, 1 + (2 * (col + 0)) + 1);
  int Rf = lines(c - 2, 1 + (2 * (col + 0)) + 0);

  int diffRcRb = std::abs(Rc - Rb);
  int diffRfRb = std::abs(Rf - Rb);
  int diffRdRb = std::abs(Rd - Rb);

  int Term0 = 2 * Rb;
  int Term1;
  int Term2;
  if (diffRcRb > std::max(diffRfRb, diffRdRb)) {
    Term1 = Rf;
    Term2 = Rd;
  } else {
    if (diffRdRb > std::max(diffRcRb, diffRfRb)) {
      Term1 = Rf;
    } else {
      Term1 = Rd;
    }
    Term2 = Rc;
  }

  int interp_val = Term0 + Term1 + Term2;
  interp_val >>= 2;

  return {.v1 = Rb - Rf,
          .v2 = Rc - Rb,
          .interp_val = interp_val,
          .edge_metric = diffRfRb + diffRcRb};
}

__attribute__((always_inline)) inline FujiPrediction
fuji_compressed_block::fuji_decode_interpolation_odd_inner(xt_lines c,
                                                           int col) const {
  int Ra = lines(c + 0, 1 + (2 * (col + 0)) + 0);
  int Rb = lines(c - 1, 1 + (2 * (col + 0)) + 1);
  int Rc = lines(c - 1, 1 + (2 * (col + 0)) + 0);
  int Rd = lines(c - 1, 1 + (2 * (col + 1)) + 0);
  int Rg = lines(c + 0, 1 + (2 * (col + 1)) + 0);

  int interp_val = (Ra + Rg);
  if (auto [min, max] = std::minmax(Rc, Rd); Rb < min || Rb > max) {
    interp_val += 2 * Rb;
    interp_val >>= 1;
  }
  interp_val >>= 1;

  return {.v1 = Rb - Rc,
          .v2 = Rc - Ra,
          .interp_val = interp_val,
          .edge_metric = std::abs(Rb - Rc) + std::abs(Rc - Ra)};
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_sample_even(
    xt_lines c, int col, std::array<int_pair, 41>& grads) {
  const FujiPrediction prediction =
      fuji_decode_interpolation_even_inner(c, col);
  const int grad =
      fuji_quant_gradient(common_info.qtable, prediction.v1, prediction.v2);
  return fuji_decode_sample(common_info.qtable, grad, prediction.interp_val,
                            grads, /*lossy=*/false);
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_sample_odd(xt_lines c, int col,
                                              std::array<int_pair, 41>& grads) {
  const FujiPrediction prediction = fuji_decode_interpolation_odd_inner(c, col);
  const int grad =
      fuji_quant_gradient(common_info.qtable, prediction.v1, prediction.v2);
  return fuji_decode_sample(common_info.qtable, grad, prediction.interp_val,
                            grads, /*lossy=*/false);
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_lossy_sample(FujiPrediction prediction,
                                                FujiGradients& gradients) {
  invariant(current_lossy_main_qtable != nullptr);

  const FujiQTable* qtable = current_lossy_main_qtable;
  const int qbase = qtable->q_base;

  for (int i = 1; qbase >= i && i < 4; ++i) {
    const FujiQTable& candidate =
        common_info.lossy_qtables.staticLossyQTable(i - 1);
    if (prediction.edge_metric <= candidate.max_grad) {
      const int grad =
          fuji_quant_gradient(candidate, prediction.v1, prediction.v2);
      return fuji_decode_sample(candidate, grad, prediction.interp_val,
                                gradients.lossyStatic[i - 1],
                                /*lossy=*/true);
    }
  }

  const int grad = fuji_quant_gradient(*qtable, prediction.v1, prediction.v2);
  return fuji_decode_sample(*qtable, grad, prediction.interp_val,
                            gradients.main, /*lossy=*/true);
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_lossy_sample_even(xt_lines c, int col,
                                                     FujiGradients& gradients) {
  return fuji_decode_lossy_sample(fuji_decode_interpolation_even_inner(c, col),
                                  gradients);
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_lossy_sample_odd(xt_lines c, int col,
                                                    FujiGradients& gradients) {
  return fuji_decode_lossy_sample(fuji_decode_interpolation_odd_inner(c, col),
                                  gradients);
}

__attribute__((always_inline)) inline int
fuji_compressed_block::fuji_decode_interpolation_even(xt_lines c,
                                                      int col) const {
  return fuji_decode_interpolation_even_inner(c, col).interp_val;
}

void fuji_compressed_block::fuji_extend_generic(int start, int end) const {
  for (int i = start; i <= end; i++) {
    lines(i, 0) = lines(i - 1, 1);
    lines(i, lines.width() - 1) = lines(i - 1, lines.width() - 2);
  }
}

void fuji_compressed_block::fuji_extend_red() const {
  fuji_extend_generic(R2, R4);
}

void fuji_compressed_block::fuji_extend_green() const {
  fuji_extend_generic(G2, G7);
}

void fuji_compressed_block::fuji_extend_blue() const {
  fuji_extend_generic(B2, B4);
}

template <typename TEven, typename TOdd>
__attribute__((always_inline)) inline void
fuji_compressed_block::fuji_decode_block(TEven func_even, TOdd func_odd,
                                         [[maybe_unused]] int cur_line) {
  invariant(common_info.line_width % 2 == 0);
  const int line_width = common_info.line_width / 2;

  auto pass = [this, &line_width, func_even,
               func_odd](std::array<xt_lines, 2> c, int row) {
    int grad = row % 3;

    struct ColorPos final {
      int even = 0;
      int odd = 0;
    };

    std::array<ColorPos, 2> pos;
    for (int i = 0; i != line_width + 4; ++i) {
      if (i < line_width) {
        for (int comp = 0; comp != 2; comp++) {
          int& col = pos[comp].even;
          int sample = func_even(c[comp], col, grad_even[grad], row, i, comp);
          lines(c[comp], 1 + (2 * col) + 0) = implicit_cast<uint16_t>(sample);
          ++col;
        }
      }

      if (i >= 4) {
        for (int comp = 0; comp != 2; comp++) {
          int& col = pos[comp].odd;
          int sample = func_odd(c[comp], col, grad_odd[grad]);
          lines(c[comp], 1 + (2 * col) + 1) = implicit_cast<uint16_t>(sample);
          ++col;
        }
      }
    }
  };

  using Tag = BayerTag;
  const std::array<CFAColor, MCU<Tag>.x * MCU<Tag>.y> CFAData =
      getAsCFAColors(BayerPhase::RGGB);
  const Array2DRef<const CFAColor> CFA(CFAData.data(), MCU<Tag>.x, MCU<Tag>.y);

  std::array<int, 3> PerColorCounter;
  std::fill(PerColorCounter.begin(), PerColorCounter.end(), 0);
  auto ColorCounter = [&PerColorCounter](CFAColor c) -> int& {
    switch (c) {
      using enum CFAColor;
    case RED:
    case GREEN:
    case BLUE:
      return PerColorCounter[static_cast<uint8_t>(c)];
    default:
      __builtin_unreachable();
    }
  };

  auto CurLineForColor = [&ColorCounter](CFAColor c) {
    xt_lines res;
    switch (c) {
      using enum CFAColor;
    case RED:
      res = R2;
      break;
    case GREEN:
      res = G2;
      break;
    case BLUE:
      res = B2;
      break;
    default:
      __builtin_unreachable();
    }
    int& off = ColorCounter(c);
    res = static_cast<xt_lines>(res + off);
    ++off;
    return res;
  };

  for (int row = 0; row != 6; ++row) {
    CFAColor c0 = CFA(row % CFA.height(), /*col=*/0);
    CFAColor c1 = CFA(row % CFA.height(), /*col=*/1);
    pass({CurLineForColor(c0), CurLineForColor(c1)}, row);
    for (CFAColor c : {c0, c1}) {
      switch (c) {
      case CFAColor::RED:
        fuji_extend_red();
        break;
      case CFAColor::GREEN:
        fuji_extend_green();
        break;
      case CFAColor::BLUE:
        fuji_extend_blue();
        break;
      default:
        __builtin_unreachable();
      }
    }
  }
}

void fuji_compressed_block::xtrans_decode_block(int cur_line) {
  fuji_decode_block(
      [this](xt_lines c, int col, FujiGradients& grads, int row, int i,
             int comp) {
        if ((comp == 0 && (row == 0 || (row == 2 && i % 2 == 0) ||
                           (row == 4 && i % 2 != 0) || row == 5)) ||
            (comp == 1 && (row == 1 || row == 2 || (row == 3 && i % 2 != 0) ||
                           (row == 5 && i % 2 == 0))))
          return fuji_decode_interpolation_even(c, col);
        invariant((comp == 0 && (row == 1 || (row == 2 && i % 2 != 0) ||
                                 row == 3 || (row == 4 && i % 2 == 0))) ||
                  (comp == 1 && (row == 0 || (row == 3 && i % 2 == 0) ||
                                 row == 4 || (row == 5 && i % 2 != 0))));
        if (header.isLossy())
          return fuji_decode_lossy_sample_even(c, col, grads);
        return fuji_decode_sample_even(c, col, grads.main);
      },
      [this](xt_lines c, int col, FujiGradients& grads) {
        if (header.isLossy())
          return fuji_decode_lossy_sample_odd(c, col, grads);
        return fuji_decode_sample_odd(c, col, grads.main);
      },
      cur_line);
}

void fuji_compressed_block::fuji_bayer_decode_block(int cur_line) {
  fuji_decode_block(
      [this](xt_lines c, int col, FujiGradients& grads,
             [[maybe_unused]] int row, [[maybe_unused]] int i,
             [[maybe_unused]] int comp) {
        if (header.isLossy())
          return fuji_decode_lossy_sample_even(c, col, grads);
        return fuji_decode_sample_even(c, col, grads.main);
      },
      [this](xt_lines c, int col, FujiGradients& grads) {
        if (header.isLossy())
          return fuji_decode_lossy_sample_odd(c, col, grads);
        return fuji_decode_sample_odd(c, col, grads.main);
      },
      cur_line);
}

void fuji_compressed_block::fuji_decode_strip(const FujiStrip& strip) {
  const unsigned line_size = sizeof(uint16_t) * (common_info.line_width + 2);

  struct i_pair final {
    int a;
    int b;
  };

  const std::array<i_pair, 3> colors = {{{R0, 5}, {G0, 8}, {B0, 5}}};

  for (int cur_line = 0; cur_line < strip.height(); cur_line++) {
    if (header.isLossy())
      setLossyQBase(strip.qbaseForLine(cur_line));

    if (header.raw_type == 16) {
      xtrans_decode_block(cur_line);
    } else {
      fuji_bayer_decode_block(cur_line);
    }

    if (header.raw_type == 16) {
      copy_line_to_xtrans(strip, cur_line);
    } else {
      copy_line_to_bayer(strip, cur_line);
    }

    if (cur_line + 1 == strip.height())
      break;

    // Last two lines of each color become the first two lines.
    for (auto i : colors) {
      memcpy(&lines(i.a, 0), &lines(i.a + i.b - 2, 0), 2 * line_size);
    }

    for (auto i : colors) {
      const auto out = CroppedArray2DRef(
          lines, /*offsetCols=*/0, /*offsetRows=*/i.a + 2,
          /*croppedWidth=*/lines.width(), /*croppedHeight=*/i.b - 2);

      // All other lines of each color become uninitialized.
      MSan::Allocated(out);

      // And the first (real, uninitialized) line of each color gets the content
      // of the last helper column from the last decoded sample of previous
      // line of that color.
      lines(i.a + 2, lines.width() - 1) = lines(i.a + 2 - 1, lines.width() - 2);
    }
  }
}

class FujiDecompressorImpl final {
  RawImage mRaw;
  const Array1DRef<const Array1DRef<const uint8_t>> strips;
  const Optional<Array1DRef<const Array1DRef<const uint8_t>>> qbaseStrips;

  const FujiDecompressor::FujiHeader& header;

  const fuji_compressed_params common_info;

  void decompressThread() const noexcept;

public:
  FujiDecompressorImpl(
      RawImage mRaw, Array1DRef<const Array1DRef<const uint8_t>> strips,
      Optional<Array1DRef<const Array1DRef<const uint8_t>>> qbaseStrips_,
      const FujiDecompressor::FujiHeader& h);

  void decompress();
};

FujiDecompressorImpl::FujiDecompressorImpl(
    RawImage mRaw_, Array1DRef<const Array1DRef<const uint8_t>> strips_,
    Optional<Array1DRef<const Array1DRef<const uint8_t>>> qbaseStrips_,
    const FujiDecompressor::FujiHeader& h_)
    : mRaw(std::move(mRaw_)), strips(strips_), qbaseStrips(qbaseStrips_),
      header(h_), common_info(header, qbaseStrips_) {}

void FujiDecompressorImpl::decompressThread() const noexcept {
  fuji_compressed_block block_info(mRaw->getU16DataAsUncroppedArray2DRef(),
                                   header, common_info);

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
  for (int block = 0; block < header.blocks_in_row; ++block) {
    try {
      Optional<Array1DRef<const uint8_t>> qbase;
      if (qbaseStrips)
        qbase = (*qbaseStrips)(block);
      FujiStrip strip(header, block, strips(block), qbase);
      block_info.reset();
      block_info.pump = BitStreamerMSB(strip.input);
      block_info.fuji_decode_strip(strip);
    } catch (const RawspeedException& err) {
      // Propagate the exception out of OpenMP magic.
      mRaw->setError(err.what());
    } catch (...) {
      // We should not get any other exception type here.
      __builtin_unreachable();
    }
  }
}

void FujiDecompressorImpl::decompress() {
#ifdef HAVE_OPENMP
#pragma omp parallel default(none)                                             \
    num_threads(rawspeed_get_number_of_processor_cores())
#endif
  decompressThread();

  std::string firstErr;
  if (mRaw->isTooManyErrors(1, &firstErr)) {
    ThrowRDE("Too many errors encountered. Giving up. First Error:\n%s",
             firstErr.c_str());
  }
}

} // namespace

FujiDecompressor::FujiDecompressor(RawImage img, ByteStream input_)
    : mRaw(std::move(img)), input(input_) {
  if (mRaw->getCpp() != 1 || mRaw->getDataType() != RawImageType::UINT16 ||
      mRaw->getBpp() != sizeof(uint16_t))
    ThrowRDE("Unexpected component count / data type");

  input.setByteOrder(Endianness::big);

  header = FujiHeader(input);
  if (!header)
    ThrowRDE("compressed RAF header check");

  if (mRaw->dim != iPoint2D(header.raw_width, header.raw_height))
    ThrowRDE("RAF header specifies different dimensions!");

  if (header.isLossless() && 12 == header.raw_bits) {
    ThrowRDE("Aha, finally, a 12-bit compressed RAF! Please consider providing "
             "samples on <https://raw.pixls.us/>, thanks!");
  }

  if (mRaw->cfa.getSize() == iPoint2D(6, 6)) {
    Optional<XTransPhase> p = getAsXTransPhase(mRaw->cfa);
    if (!p)
      ThrowRDE("Invalid X-Trans CFA");
    if (p != iPoint2D(0, 0))
      ThrowRDE("Unexpected X-Trans phase: {%i,%i}. Please file a bug!", p->x,
               p->y);
  } else if (mRaw->cfa.getSize() == iPoint2D(2, 2)) {
    Optional<BayerPhase> p = getAsBayerPhase(mRaw->cfa);
    if (!p)
      ThrowRDE("Invalid Bayer CFA");
    if (p != BayerPhase::RGGB)
      ThrowRDE("Unexpected Bayer phase: %i. Please file a bug!",
               static_cast<int>(*p));
  } else {
    ThrowRDE("Unexpected CFA size");
  }

  // read block sizes
  std::vector<uint32_t> block_sizes;
  block_sizes.resize(header.blocks_in_row);
  for (auto& block_size : block_sizes)
    block_size = input.getU32();

  // some padding?
  if (const uint64_t raw_offset = sizeof(uint32_t) * header.blocks_in_row;
      raw_offset & 0xC) {
    const int padding = 0x10 - (raw_offset & 0xC);
    input.skipBytes(padding);
  }

  if (header.isLossy()) {
    const int qbaseStride = implicit_cast<int>(roundUp(header.total_lines, 16));
    const Buffer::size_type qbaseBytes =
        implicit_cast<Buffer::size_type>(header.blocks_in_row) *
        implicit_cast<Buffer::size_type>(qbaseStride);
    const auto qbaseAll = input.getStream(qbaseBytes).getAsArray1DRef();

    qbaseStrips.reserve(header.blocks_in_row);
    for (int strip = 0; strip != header.blocks_in_row; ++strip) {
      qbaseStrips.emplace_back(
          qbaseAll.getCrop(strip * qbaseStride, header.total_lines)
              .getAsArray1DRef());
    }
  }

  // calculating raw block offsets
  strips.reserve(header.blocks_in_row);

  for (const auto& block_size : block_sizes)
    strips.emplace_back(input.getStream(block_size).getAsArray1DRef());
}

void FujiDecompressor::decompress() const {
  Optional<Array1DRef<const Array1DRef<const uint8_t>>> qbaseStripRefs;
  if (header.isLossy()) {
    invariant(qbaseStrips.size() == strips.size());
    qbaseStripRefs = Array1DRef<const Array1DRef<const uint8_t>>(
        qbaseStrips.data(), implicit_cast<int>(qbaseStrips.size()));
  }

  FujiDecompressorImpl impl(
      mRaw,
      Array1DRef<const Array1DRef<const uint8_t>>(
          strips.data(), implicit_cast<Buffer::size_type>(strips.size())),
      qbaseStripRefs, header);

  impl.decompress();
}

FujiDecompressor::FujiHeader::FujiHeader(ByteStream& bs)
    : signature(bs.getU16()), version(bs.getByte()), raw_type(bs.getByte()),
      raw_bits(bs.getByte()), raw_height(bs.getU16()),
      raw_rounded_width(bs.getU16()), raw_width(bs.getU16()),
      block_size(bs.getU16()), blocks_in_row(bs.getByte()),
      total_lines(bs.getU16()),
      MCU(raw_type == 16 ? ::rawspeed::MCU<XTransTag>
                         : ::rawspeed::MCU<BayerTag>) {}

FujiDecompressor::FujiHeader::operator bool() const {
  // general validation
  const bool invalid =
      (signature != 0x4953 || (!isLossless() && !isLossy()) ||
       raw_height > 0x3000 || raw_height < FujiStrip::lineHeight() ||
       raw_height % FujiStrip::lineHeight() || raw_width > 0x3000 ||
       raw_width < 0x300 || raw_width % 24 || raw_rounded_width > 0x3000 ||
       block_size != 0x300 || raw_rounded_width < block_size ||
       raw_rounded_width % block_size ||
       raw_rounded_width - raw_width >= block_size || blocks_in_row > 0x10 ||
       blocks_in_row == 0 || blocks_in_row != raw_rounded_width / block_size ||
       blocks_in_row != roundUpDivisionSafe(raw_width, block_size) ||
       total_lines > 0x800 || total_lines == 0 ||
       total_lines != raw_height / FujiStrip::lineHeight() ||
       (raw_bits != 12 && raw_bits != 14 && raw_bits != 16) ||
       (raw_type != 16 && raw_type != 0));

  return !invalid;
}

} // namespace rawspeed
