Test build of [PR #1](https://github.com/M0LTE/wsjtx/pull/1) (WSPR Max decode depth), built from commit `@SHA@`. **The PR is not merged.** This release exists so the branch can be tried on air before that decision is made.

WSJT-X 3.0.2 plus the PR. Everything except WSPR decoding is stock, and Fast, Normal and Deep are unchanged. The new **Max** depth appears on the Decode menu in WSPR mode only.

## Downloads

| Platform | File |
|---|---|
| Windows x64 | `wsjtx-@VERSION@-win64.exe` |
| Linux x64 | `wsjtx-@VERSION@-linux-x86_64.AppImage`, plus .deb and .rpm |
| Linux ARM64 | `wsjtx-@VERSION@-linux-aarch64.AppImage`, plus .deb and .rpm |

The AppImages are the portable option: `chmod +x` the file and run it. Nothing is installed and an existing WSJT-X is left alone. The .deb and .rpm install into /usr and will replace a packaged WSJT-X.

Help > About reports `@VERSION@`, so a report against this build can be told apart from one against stock 3.0.2 or against a later test build.

## Not identical to the PR

The branch as submitted does not compile on Windows. It makes the OSD pattern boxes per-thread with `!$omp threadprivate` on a COMMON block, and gfortran emits that as a common TLS symbol, which the assembler used on Windows has no directive for. This build holds those arrays in a module instead: same three arrays, same per-thread storage, and it assembles on all three platforms. Linux binaries are otherwise the PR's code.

## Known limitations

- The Windows installer is signed with a throwaway certificate that chains to nothing, so SmartScreen will warn. That is expected for a test build.
- Max is CPU hungry by design: roughly 10 to 30 s of a 2-minute cycle on a busy band across 12 cores.

Checksums are in `SHA256SUMS.txt`.
