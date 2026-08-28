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

## Relationship to the PR

The branch as originally submitted did not compile on Windows: it made the OSD pattern boxes per-thread with `!$omp threadprivate` on a COMMON block, and gfortran emits that as a common TLS symbol, which the assembler used on Windows has no directive for. Those arrays now live in a module instead, which assembles everywhere and was separately confirmed to produce byte-identical decodes (stdout and ALL_WSPR.TXT, deep search included) against the previously validated binary over off-air corpus slots at full Max depth.

That fix is now part of PR #1 itself, so the decoder in these installers is the PR's code. The only thing here that is not in the PR is the build plumbing: the workflow that produces these files, and the version suffix that puts `@VERSION@` in the title bar.

One thing to keep in mind when reading the PR evidence: every off-air number quoted there was produced by binaries built before that fix.

## Known limitations

- The Windows installer is signed with a throwaway certificate that chains to nothing, so SmartScreen will warn. That is expected for a test build.
- Max is CPU hungry by design: roughly 10 to 30 s of a 2-minute cycle on a busy band across 12 cores.

Checksums are in `SHA256SUMS.txt`.
