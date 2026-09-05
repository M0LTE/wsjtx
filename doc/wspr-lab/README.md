# The WSPR receive-sensitivity lab record

This directory is the written record of the measurement campaign behind the WSPR changes on this branch: what was built, every number that led to a decision, and, deliberately, everything that was tried and did not survive. The raw per-run JSON, the channel simulator, the harnesses and the corpus tooling live in the external `wspr-lab` working directory they were run from; these files are the findings.

Reading order:

- `SUMMARY.txt` is the running scoreboard of the whole campaign, oldest first.
- `NOTES-rounds-two-three.md` is the detailed record of the second and third rounds: the sync gates, the fade-weighted retries, the extra pass, the blanker-sweep fabrications and their fix, the wrong-codeword guard, return-visitor recall, and the off-air validation over six corpuses ending with the four unseen night bands.
- `NOTES-fano-budget.md` closes the Fano cycle-budget question: no acceptance-side gate rescues a larger `-C`, and `-C 500` recorded zero false accepts in 9,600 paired trials.

The remaining files are single-topic studies, most of them negative results kept so nobody spends the time twice:

- `NOTES-coherent-sync-metric.md`: a block-coherent sync statistic (0.10 dB, AWGN only; not shipped, but its control arm exposed the refinement-gate finding that later shipped as `-S`).
- `NOTES-adaptive-coherent.md`: measuring the channel's coherence time and using the matched smoothing window (works, still loses to trying eight windows).
- `NOTES-interference-excision.md`: the narrowband excision study behind `-X`.
- `NOTES-fork-and-literature-review.md`: every WSPR fork and downstream stack examined; no one else had changed the decoder since 2020.
- `NOTES-fst4w-comparison.md`: FST4W-120 measured on the same harness, the yardstick for what a clean-sheet redesign buys.
- `NOTES-jtdx-candidate-window.md`: JTDX's candidate-detection change, measured to do nothing here.

Off-air method, used throughout: recordings from a GPSDO-locked receiver are cut into two-minute slots and decoded in time order, and every printed spot is checked against the worldwide WSPR logs; a spot some other receiver also logged in the same two minutes is corroborated, and a second tier separates calls active on the band that day from the truly unknown. Fabrication counts in these notes are spots the whole world failed to corroborate, minus the ones shown to be real but locally-unique stations.
