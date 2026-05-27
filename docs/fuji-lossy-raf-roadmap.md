# Fujifilm lossy compressed RAF roadmap

This is a fork-local working roadmap for implementing Fujifilm lossy / non-lossless compressed RAF support in RawSpeed.

It is intended only for `EffortlessSteven/rawspeed` fork-local implementation work. It does not direct upstream work. Any future upstreaming decision is outside this roadmap and will be handled separately by Steven.

## Goal

Implement RawSpeed support for Fujifilm compressed RAF header `version == 0` files in `FujiDecompressor`, while preserving existing behavior for:

- Fujifilm lossless compressed RAF files, header `version == 1`;
- uncompressed RAF routing;
- existing supported X-Trans and Bayer files.

The durable fix belongs in RawSpeed’s Fuji decoder. It should not be implemented as a darktable-side LibRaw fallback.

## Non-goals

Do not include the following in this work:

- darktable fallback routing;
- LibRaw-as-product-path behavior;
- camera XML changes unless separately required;
- raw sample files committed to git;
- broad decoder rewrites;
- formatting-only churn mixed with decoder changes;
- opening upstream PRs from this roadmap;
- upstreaming this roadmap document;
- asking Codex to perform upstream submission work.

## Format summary

Current RawSpeed supports Fuji compressed RAF header `version == 1`, which corresponds to Fuji lossless compressed RAF files.

The unsupported files use Fuji compressed RAF header `version == 0`, which corresponds to Fuji lossy / non-lossless compressed RAF files.

The container layout differs by version.

Version 1 layout:

```text
Fuji compressed header
block-size table
alignment padding
compressed block streams
```

Version 0 layout:

```text
Fuji compressed header
block-size table
alignment padding
q_base table
compressed block streams
```

The version-0 q-base table layout is:

```cpp
qbase_stride = roundUp(total_lines, 16);
qbase_bytes  = blocks_in_row * qbase_stride;
```

The q-base value must be treated as `uint8_t`. Observed samples include q-base values above signed-byte range.

Known decoder paths:

```text
raw_type == 16  # X-Trans
raw_type == 0   # Bayer
```

Both paths must be validated separately.

## Sample coverage

Do not commit raw samples to git. Track samples by SHA256 and parsed metadata, not filename.

The local fixture set should cover these categories:

### X-Trans lossy compressed, version 0

- X-H2S 14-bit
- X-T4 14-bit
- X-T5 14-bit / 40MP
- X-E5 14-bit / 40MP

### X-Trans lossless compressed, version 1 regression

- X-H2S
- X-T4
- X-T5
- X-E5
- X-T20

### Bayer lossy compressed, version 0

- GFX50S II 14-bit, including narrow final strip
- GFX100S 14-bit
- GFX100S 16-bit, including high q-base values
- GFX100S II 16-bit
- GFX100 II 14-bit issue sample
- GFX100 II 16-bit issue sample

### Bayer lossless compressed, version 1 regression

- GFX50S
- GFX50S II
- GFX100
- GFX100S 14-bit
- GFX100S 16-bit
- GFX100S II 14-bit
- GFX100S II 16-bit

### Routing baselines

- X-series uncompressed RAF files
- GFX uncompressed RAF files
- older packed 12-bit / 14-bit RAF files
- legacy RAF paths as available

## Validation principles

A file opening is not sufficient proof.

Decoder work should be validated with:

- parsed header and layout checks;
- raw-plane comparison against a known-good reference decoder where practical;
- RawSpeed sample hashes once output is trusted;
- lossless regression checks;
- OpenMP determinism checks.

For threaded decode paths, compare output under multiple thread counts:

```bash
OMP_NUM_THREADS=1
OMP_NUM_THREADS=2
OMP_NUM_THREADS=8
```

Expected result:

```text
same decoded output hash
```

## Local branch strategy

Use a fork-local integration branch:

```text
fuji-lossy-raf
```

All implementation PRs should target:

```text
EffortlessSteven/rawspeed:fuji-lossy-raf
```

not `develop`.

This roadmap governs fork-local work only. Codex should not open, prepare, or modify upstream PRs from this roadmap.

After the fork-local integration branch is complete and proven, Steven will decide separately whether and how to upstream the work. Any future upstreaming will be handled manually by Steven, one carefully scoped PR at a time.

## PR stack

### PR 1 — distinguish Fuji compressed RAF versions

Title:

```text
FujiDecompressor: distinguish compressed RAF versions
```

Goal:

Recognize Fuji compressed RAF header versions explicitly.

Expected behavior:

```text
version == 1  # lossless compressed, existing path
version == 0  # lossy / non-lossless compressed, known but unsupported initially
```

Allowed changes:

- add `FujiHeader::isLossless()`;
- add `FujiHeader::isLossy()`;
- adjust header validation to recognize `version == 0` as a known Fuji compressed RAF version;
- reject version-0 decode with a precise unsupported error before existing version-1 logic consumes the payload.

Do not include:

- q-base parsing;
- q-table refactors;
- lossy decode;
- LibRaw code;
- sample archive changes;
- darktable changes.

Definition of done:

- version-1 behavior remains unchanged;
- version-0 files fail with a precise unsupported message;
- unknown versions still fail;
- diff is limited to `FujiDecompressor` files unless a small test/diagnostic change is justified.

### PR 2 — parse version-0 q_base table

Title:

```text
FujiDecompressor: parse lossy RAF q_base table
```

Goal:

Parse the version-0 q-base side table after the block-size table and padding, before compressed block streams.

Expected behavior:

- compute `qbase_stride = roundUp(total_lines, 16)`;
- compute `qbase_bytes = blocks_in_row * qbase_stride`;
- store q-base views or equivalent data per strip;
- record observed q-base values as unsigned bytes;
- still reject version-0 decode until the decoder is implemented.

Definition of done:

- stream offsets are computed, not hard-coded;
- version-1 stream layout remains unchanged;
- version-0 q-base table is consumed safely;
- malformed or truncated q-base tables fail cleanly.

### PR 3 — refactor q-table setup without behavior change

Title:

```text
FujiDecompressor: factor Fuji q-table setup
```

Goal:

Introduce explicit q-table objects while preserving existing version-1 output.

Allowed changes:

- factor q-table setup into a small internal helper;
- make q-table lookup table-specific;
- preserve current lossless behavior exactly.

Definition of done:

- version-1 decoded output remains hash-identical;
- no version-0 decode is enabled;
- no shared mutable decode state is introduced.

### PR 4 — add lossy q-table and gradient state

Title:

```text
FujiDecompressor: add lossy q-table state
```

Goal:

Add the state needed for lossy decode without enabling full pixel output yet.

Required pieces:

- dynamic main q-tables keyed by `uint8_t q_base`;
- three static lossy q-tables;
- explicit main vs static lossy gradient state;
- immutable shared q-table cache for OpenMP-safe decode.

Definition of done:

- q-base values above signed-byte range are handled correctly;
- common decode state remains immutable during threaded decode;
- version-1 output remains unchanged.

### PR 5 — decode X-Trans lossy RAF

Title:

```text
FujiDecompressor: decode X-Trans lossy RAF
```

Goal:

Enable version-0 lossy compressed RAF decoding for X-Trans files.

Required pieces:

- q-table selection per sample;
- q-base-scaled reconstruction;
- q-base-aware wrapping;
- even-prediction handling without double-shifting;
- OpenMP deterministic output.

Definition of done:

- X-H2S, X-T4, X-T5, and X-E5 lossy version-0 files decode;
- output matches raw-plane reference where practical;
- X-Trans version-1 lossless files remain unchanged.

### PR 6 — decode Bayer lossy RAF

Title:

```text
FujiDecompressor: decode Bayer lossy RAF
```

Goal:

Enable version-0 lossy compressed RAF decoding for Bayer / GFX files.

Required pieces:

- Bayer lossy decode path;
- narrow final-strip handling;
- 14-bit and 16-bit validation;
- high-q-base validation.

Definition of done:

- GFX50S II narrow final-strip file decodes;
- GFX100S 14-bit and 16-bit files decode;
- GFX100S II 16-bit file decodes;
- GFX100 II issue samples decode;
- Bayer version-1 lossless files remain unchanged;
- output is deterministic under OpenMP.

### PR 7 — validation cleanup

Title:

```text
FujiDecompressor: validate lossy RAF support
```

Goal:

Clean up the local integration branch after both decoder paths work.

Allowed changes:

- remove temporary diagnostics;
- document validation results in PR body or local notes;
- ensure sample hashes / raw-plane comparisons are recorded outside git as appropriate;
- confirm no raw sample files are committed.

Definition of done:

- normal build passes;
- normal tests pass;
- lossless regressions remain unchanged;
- lossy output matches reference expectations;
- GFX100 II issue samples open through RawSpeed;
- branch is ready for Steven’s final review; no upstream action is implied.

## Upstream boundary

This roadmap is fork-local only.

Codex must not open upstream PRs, prepare upstream branches, or modify `darktable-org/rawspeed` from this roadmap. All implementation PRs in this roadmap target only:

```text
EffortlessSteven/rawspeed:fuji-lossy-raf
```

When the fork-local branch is complete and proven, Steven will decide separately whether to upstream the work. If upstreaming happens, Steven will do it manually, one carefully scoped PR at a time.

No upstream action is part of this roadmap.
