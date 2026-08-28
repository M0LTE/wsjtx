Test build of [PR #1](https://github.com/M0LTE/wsjtx/pull/1) (WSPR Max decode depth), built from commit `@SHA@`. **The PR is not merged.** This release exists so the branch can be tried on air before that decision is made.

Built from WSJT-X 3.0.2 plus the PR. Everything except WSPR decoding is stock; Fast, Normal and Deep are unchanged. The new **Max** depth appears on the Decode menu in WSPR mode only.

## Downloads

| Platform | File |
|---|---|
| Windows x64 | `wsjtx-@VERSION@-win64.exe` |
| Linux x64 | `wsjtx-@VERSION@-linux-x86_64.AppImage`, plus .deb and .rpm |
| Linux ARM64 | `wsjtx-@VERSION@-linux-aarch64.AppImage`, plus .deb and .rpm |

The AppImages are the portable option: `chmod +x` the file and run it. Nothing is installed and an existing WSJT-X is left alone. The .deb and .rpm install into /usr and will replace a packaged WSJT-X.

## Known limitations

- The Windows installer is signed with a throwaway certificate that chains to nothing, so SmartScreen will warn. That is expected for a test build.
- Max is CPU hungry by design: roughly 10 to 30 s of a 2-minute cycle on a busy band across 12 cores.
- The version reads as 3.0.2-devel in Help > About. The git hash beside it identifies this build.

Checksums are in `SHA256SUMS.txt`.
