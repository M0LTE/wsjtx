# The Fano budget, revisited: can acceptance-side compensation buy back -C?

Date: 2026-08-29.  Question from Tom: -C is pinned at 500 because raising it
was measured to trade tenths of a dB for fabricated callsigns.  Can the
sensitivity of a larger budget be captured with an acceptance-side gate that
keeps the false rate at the branch standard (zero observed at scale)?

**Verdict: no, twice over, and both negatives are well-measured.  A wrong
codeword accepted after a big search is statistically the same object as a
real marginal decode gained by that search, so no threshold separates them;
and behind the hash gate, where acceptance IS safe, the deep-search role is
already saturated by OSD and recall, so a 100x Fano budget adds one decode in
twelve hundred trials.  -C 500 stays.  Do not ship anything.**

## What was built

src/wsprd-cf: the branch tip plus one default-off flag, verified byte-identical
to the tip build on a corpus slot with full Max arguments before any run.

    -E n   after the normal Fano budget fails, retry with budget n and accept
           the result only if the hash table vouches for the station it claims
           to be (call in hashtab, and for type 1 the stored grid must match),
           marked decodetype 2.  The same standard OSD and -r are held to.

tools/cfrun.py: paired characterization.  One AWGN recording per trial,
decoded by -C 500 / 2000 / 10000 / 50000 in isolated scratch dirs, every
ACCEPTED decode logged with its ALL_WSPR diagnostics (cycles consumed, metric,
nhardmin, decode type) and an exact real/false verdict.  Args throughout were
the shipping Max set minus the blanker sweep and recall:
-o 4 -d -N 20 -X 1 -A -S 0,-1 -Y 4 -G, single signal, fresh random callsign
per trial (hash cold).  9,600 trials, SNR -31.5 to -34.0 in 0.5 dB steps,
split across this box and hyperv-gha; 18,196 accepted decodes.

## Characterization

Real decodes and false accepts, totals over the grid:

    arm        real   false      note
    -C   500   4183       0      the shipping setting: clean at scale
    -C  2000   4462       1
    -C 10000   4723      11      +12.9% real at these SNRs (~+0.3 dB), as history said
    -C 50000   4733      83      no further net gain, see below

Every false accept was a plain Fano result (decodetype 0); -G was on, so these
95 are the ones that pass the type-2 guard.

-C 50000 gains almost nothing over 10000 overall and LOSES real decodes at the
easier SNRs (-32.0: 1229 to 1193; -31.5: 1467 to 1421).  New observation, and
the mechanism follows from the pipeline: with a huge budget an early attempt at
slightly-wrong parameters can accept a wrong codeword; if it survives the
guard it is subtracted, and subtracting a signal that is not there (or is
there at other parameters) damages the band for every later attempt.  A big
Fano budget is not merely unsafe at the output; it is self-harming inside the
search.

## The separation hypothesis: refuted

Marginal real decodes (gained at the big budget, not decoded at 500) against
false accepts, pooled over C10000 and C50000:

    population              n     cycles/bit (med)   metric (med)   nhardmin (med)
    marginal real        1296     4.7k / 18.5k        -191 / -233     31 / 32
    false accepts          94     5.9k / 30.1k        -264 / -289     36 / 36

The medians differ; the distributions lie on top of each other.  ROC over
every simple cut and combination:

    cut                     keeps marginal real    admits false
    metric >= 0                     3 / 1296            0 / 94
    cycles <= 100/bit               0 / 1296            0 / 94
    nhardmin <= 24                 30 / 1296            0 / 94
    nhardmin<=24 and metric>=0      0 / 1296            0 / 94

The best zero-false cut keeps 2.3% of the gained decodes.  Both populations
are what operating at the code's limit looks like: large negative metrics,
about a third of the symbols in disagreement, most of the budget consumed.
This is the GLRT study's conclusion arrived at from a third direction: near
the limit the true codeword is not distinguishable from a well-fitted wrong
one by any per-decode statistic, and the only working discriminators are
structural (the message class, in -G) or external (the hash table).

## Hash-gated escalation: saturated

-E 50000 on top of the full Max argument set, known-station protocol (station
decoded at -25 dB two cycles earlier, so hash table and priors warm; 300
trials per point):

    SNR      full Max     full Max -E 50000
    -33       288/300         288/300
    -34       241/300         242/300
    -35       162/300         162/300
    -36        81/300          81/300

One extra decode in twelve hundred trials.  OSD at -o 4 already performs a
far deeper search than Fano at any budget (ordered statistics over the soft
symbols rather than sequential tree probing), behind the same gate, and
recall aims it at the right frequencies; there is nothing left for escalated
Fano to find.  Zero false accepts in either arm.

## CPU, for the record

Median seconds per DECODED trial at one thread: 0.47 (C500), 0.59 (C2000),
1.21 (C10000), 4.29 (C50000); failing trials, which dominate a real band,
cost several times more at the big budgets.  Moot given the verdict, but it
underlines the design point: -C 500 with many cheap attempts beats few
expensive ones on every axis measured here.

## Conclusion

-C 500 is not a compromise; it is the right operating point.  The 0.3 dB a
20x budget buys arrives glued to a false-decode population that no
acceptance-side statistic can strip off (best case keeps 1/43 of the gain),
and where safe acceptance exists the sensitivity is already taken by OSD plus
recall.  The escalation flag stays in src/wsprd-cf as a research knob;
nothing lands in the wsjtx tree.
