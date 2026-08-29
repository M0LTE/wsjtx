# A block-coherent sync metric for wsprd

Date: 2026-08-26.  Lead from the last section of `results/jtw/NOTES.md`.

**Verdict: it works, it is small, and most of what it appears to buy is
something else.  0.29 dB on a steady path, nothing at all on any of three
fading channels, no new false decodes in 1200 busy sub-band files and 4000
pure-noise files, +10% CPU on a busy band.  Of that 0.29 dB, 0.17 dB is NOT the
coherent metric: it is the side effect of leaving `minsync1` open so that every
candidate gets the fine lag and frequency search.  The coherent metric itself is
worth 0.10 dB, on AWGN only.  I would not ship the metric.  The `minsync1`
by-product is worth a separate look, because it is a one-line change, it is
worth more than the thing it was a control for, and it produced no false
decodes here.**

Everything is generated from `src/wsprd-stage` by script.  Nothing in
`/home/tf/wsjtx` was touched, no existing `src/wsprd-*` directory was modified,
and `tools/harness.py`, `tools/wsprchan.c`, `tools/coherent.c.inc` and
`results/SUMMARY.txt` were read only.

## 0. The problem, restated

`sync_and_demodulate()` scores a candidate noncoherently.  Per symbol it takes
the magnitudes p0..p3 of the four tone correlators and accumulates

    cmet  = (p1+p3) - (p0+p2)
    ss   += +-cmet, the sign taken from the sync vector pr3[i]
    totp += p0+p1+p2+p3
    ss    = ss/totp

compared against `minsync1` = 0.10 and `minsync2` = 0.12, the latter dropping to
0.10 in pass 2.

Write A for the sum over symbols of the two tone magnitudes the sync vector
ALLOWS at that symbol and B for the two it RULES OUT.  The four tones partition
into exactly those two pairs, so the stock metric is

    ss = (A - B) / (A + B)

Taking the magnitude symbol by symbol throws away the phase relationship between
symbols, and at these SNRs that is expensive: the useful part of `ss` is
quadratic in the signal amplitude rather than linear.  A one-line model gives
`ss ~ gamma/8` with gamma the per-symbol SNR, which puts `minsync2` = 0.12 at
gamma = 0.96.  At WSPR's symbol rate that is -32.3 dB re 2500 Hz, i.e. exactly
the decode threshold, which is why the gate bites where it does.

## 1. The second idea in the brief is a no-op, and that is provable

"Normalise by the two ruled-out tones instead of by all four" is not a different
statistic.  With A and B as above,

    stock     u = (A-B)/(A+B)
    off-tone  v = (A-B)/(2B) = u/(1-u)

which is strictly increasing.  The two rank every candidate identically and gate
identically once 0.12 is mapped to 0.12/0.88 = 0.1364.  `tools/csynctheory.py`
section A checks it numerically: the orderings are identical and
`max |v - u/(1-u)|` over 2000 trials is 2.8e-17.

The same argument disposes of every variant: anything that is a monotone
function of A/B is the stock metric with a relabelled threshold.  Only a
per-symbol or per-block REWEIGHTING can differ, which is what the block-coherent
metric below is.

It also disposes of the interference argument from `results/fex/NOTES.md`.  A
carrier that lifts all four tones by c takes u to (A-B)/(A+B+4c) and v to
(A-B)/(2B+2c), and since A and B are close near the threshold those collapse by
the same factor.  Nothing was measured on this because nothing needed to be.

## 2. What was built

`tools/mkcsync.py` generates `src/wsprd-csync` from `src/wsprd-stage`.

Within a block of L symbols, combine the correlator outputs coherently and take
the maximum over the 2^L data hypotheses, exactly as
`noncoherent_sequence_detection()` does for demodulation.  Do the same over the
two tones the sync vector rules out, which are a free noise reference with the
same distribution, and form

    ss_coh = (A' - B') / (A' + B')
    A' = sum over blocks of  max over hypotheses |coherent sum, allowed tones|
    B' = sum over blocks of  max over hypotheses |coherent sum, ruled-out tones|

Same shape as the stock metric, so it is normalised the same way, has mean zero
on noise, and the max-over-hypotheses bias cancels between the two arms.

The derotation is what makes the coherent sum legal.  WSPR's tone spacing is
exactly 1/T, so over one symbol every tone advances the phase by the same amount
modulo 2*pi; multiplying by the running product of the per-symbol phase advances
lines the symbols of a block up in phase whatever data was sent.  That is the
argument `coherent_sequence_detection()` already uses in the demodulator, where
it is worth 1.4 dB, applied to the sync vector instead of to the data.

### Switches, all off by default

    -K <L>          block length, 0 = off (the default, and byte-identical)
    -M <mode>       0 = replace: the coherent metric drives the drift
                        refinement, the fine lag and frequency search and both
                        gates
                    1 = additive: the stock chain runs untouched and the
                        coherent metric is evaluated once more at the refined
                        parameters; a candidate passes if EITHER gate passes
                    2 = the coherent metric steers the search but the STOCK
                        metric still gates.  Not a shipping mode; it splits
                        "the estimate got better" from "the gate got better"
    -T <t1,t2[,t3]> thresholds replacing minsync1 and minsync2, t3 being the
                    pass-2 value because stock minsync2 drops to 0.10 there
    -Y <hz>         read-only diagnostics: the true signal offset from 1500 Hz

The first two calls to `sync_and_demodulate` are always left on the stock
metric.  They only have to get into the right basin, and at that point the
frequency is known only to the 0.73 Hz candidate grid, far too coarse for a
block of symbols to add coherently.  After the coarse frequency search the
estimate is good to +-0.125 Hz and after the fine search to +-0.025 Hz.

Cost: one extra correlator pass per candidate, against the twenty-five the
refinement chain already spends, plus 2 * 2^L * L multiply-adds per block, which
is nothing.

### Arms measured

    base    the branch decoder, -K absent
    add3    -K 3 -M 1 -T -1,0.1865,0.1435     the main arm
    add2    -K 2 -M 1 -T -1,0.1702,0.1326
    add3c   -K 3 -M 1 -T 99,0.1865,0.1435     as add3 but gate one left stock
    coh3    -K 3 -M 0 -T -1,0.1961,0.1561     coherent metric drives everything
    est3    -K 3 -M 2 -T 0.1561,0.12,0.10     coherent search, stock gate
    ref0    -K 3 -M 1 -T -1,99,99             CONTROL: the coherent metric is
                                              computed but can never open the
                                              gate, while t1 = -1 leaves gate
                                              one permanently open.  Anything
                                              this arm gains is the extra
                                              refinement, not the metric.

The thresholds are not guesses.  `tools/csyncroc.py` measures the distribution
of every metric over a false-alarm population -- 2000 single-signal files at
-45 dB, far below decodable, plus 150 busy sub-band files counting only
candidates more than 6 Hz from any transmission -- and places each arm's
threshold the same number of standard deviations above its own false-alarm mean
as 0.12 sits above the stock one.  That is 5.48 sd for `minsync2` and 3.64 sd
for the 0.10 used in pass 2.

A tail-probability match would be better and is not measurable, which is itself
worth knowing: **the `min_snr` floor in the candidate detector already keeps
noise out, so 0.12 is beyond every one of the 1774 false-alarm candidates
collected.  On an ordinary file the sync gate is not what rejects noise.**  Four
thousand pure-noise files produce essentially no candidates at all.

## 3. What the model predicted

`tools/csynctheory.py` reimplements the metric in numpy: no decoder, no wav
files.  SNR needed for a deflection of 5 noise standard deviations, per symbol,
AWGN, no residual frequency error:

    stock, noncoherent    +0.12 dB
    block-coherent L=2    -0.88 dB    1.01 dB better
    block-coherent L=3    -1.36 dB    1.48 dB
    block-coherent L=4    -1.68 dB    1.81 dB
    block-coherent L=6    -2.20 dB    2.33 dB

with two limits that decided the design.  A residual frequency error destroys
it: at 0.125 Hz, where the coarse search leaves you, L=3 keeps only 0.83 dB and
L=6 is worse than the stock metric.  Fading destroys it the same way: at a
1.0 Hz Doppler spread L=3 LOSES 0.21 dB and only L=2 still gains anything.  Both
are the statement that the channel has to hold still for L*0.68 s.

So the prediction was 1 to 1.5 dB at the gate on a steady path and nothing on a
fast one.  The gate part of that did not survive contact with the decoder; see
section 6.

## 4. The feature is off by default

`tools/csyncident.py`, 216 files across six SNRs and three channels:

    bin/wsprd-csync with no -K  vs  bin/wsprd-branch   0 mismatches, stdout
                                                       compared byte for byte
    -K 3 -M 1 -T -1,0.1865,0.1435, -P 1 vs -P 8        0 mismatches
    additive mode lost a decode the stock gate found   0 files

## 5. The measurement agrees with the one that started this

The baseline arm reproduces `results/jtw/NOTES.md` to three decimal places, on
the same wav files: `tools/csyncrun.py` seeds each trial from the SNR value with
the same `seed0` and the same generation order `tools/jtwrun.py` used.

| channel | SNR | det | syn | dec | jtw NOTES |
|---|---|---|---|---|---|
| AWGN         | -31.0 | 0.996 | 0.956 | 0.948 | 0.996 / 0.956 / 0.948 |
| AWGN         | -32.5 | 0.972 | 0.520 | 0.400 | 0.972 / 0.520 / 0.400 |
| AWGN         | -33.5 | 0.878 | 0.258 | 0.108 | 0.878 / 0.258 / 0.108 |
| Rayleigh 0.1 | -31.0 | 0.996 | 0.888 | 0.164 | 0.996 / 0.888 / 0.164 |
| Rayleigh 1.0 | -26.0 | 1.000 | 0.994 | 0.648 | 1.000 / 0.994 / 0.648 |
| Rayleigh 1.0 | -27.5 | 1.000 | 0.910 | 0.084 | 1.000 / 0.910 / 0.084 |

The 50% thresholds reproduce `results/SUMMARY.txt` as well: -32.36 against
-32.30, -29.91 against -29.90, -28.38 against -28.40, -26.30 against -26.30.

## 6. The diagnostic: the change moves what it aimed at

500 files per point, `-C 500 -o 4 -d -N 14 -P 1`.  `det` a candidate formed
within 1 Hz of the truth in some pass, `syn` such a candidate also survived the
sync gate, `dec` the message decoded, `nwat` candidates handed to the
demodulator per file.

| channel | SNR | arm | det | syn | dec | nwat |
|---|---|---|---|---|---|---|
| AWGN | -33.5 | base  | 0.878 | 0.258 | 0.108 | 0.36 |
| AWGN | -33.5 | add3c | 0.878 | 0.312 | 0.114 | 0.42 |
| AWGN | -33.5 | ref0  | 0.878 | 0.328 | 0.144 | 0.44 |
| AWGN | -33.5 | add3  | 0.878 | 0.378 | 0.168 | 0.50 |
| AWGN | -32.5 | base  | 0.972 | 0.520 | 0.400 | 0.80 |
| AWGN | -32.5 | ref0  | 0.972 | 0.630 | 0.500 | 0.93 |
| AWGN | -32.5 | add3  | 0.972 | 0.694 | 0.558 | 1.01 |
| AWGN | -31.0 | base  | 0.996 | 0.956 | 0.948 | 2.03 |
| AWGN | -31.0 | add3  | 0.996 | 0.982 | 0.974 | 2.16 |
| R 0.1 | -31.0 | base | 0.996 | 0.888 | 0.164 | 1.77 |
| R 0.1 | -31.0 | add3 | 0.996 | 0.970 | 0.172 | 1.98 |
| R 1.0 | -27.5 | base | 1.000 | 0.910 | 0.084 | 1.80 |
| R 1.0 | -27.5 | add3 | 1.000 | 0.958 | 0.086 | 1.97 |

On AWGN at -33.5 dB the gate now passes 37.8% of files instead of 25.8% and the
decode rate follows, 10.8% to 16.8%.  That is the mechanism working exactly as
intended: `det` is untouched, `syn` moves, `dec` moves with it.

On both fading channels `syn` also moves -- 0.888 to 0.970, 0.910 to 0.958 --
**and the decode rate does not follow at all**, 0.164 to 0.172 and 0.084 to
0.086.  On a fading channel the sync gate is not the binding constraint; the
demodulator is, and letting more candidates past a gate they were already
clearing 89% of the time buys nothing.  That is why the whole idea is an
AWGN-only result, and it was visible in the baseline breakdown before any of
this was built.

## 7. 50% decode thresholds

dB re 2500 Hz, callsign not in the hash table.  3000 files per arm per channel,
six SNR points on a 0.5 dB grid, 500 trials each.  The grid is split between two
boxes on the whole and half dB and merged; every arm at a given SNR sees
byte-identical wav files, so the arms are paired.  95% CI from a profile
likelihood, `tools/fit_ci.py`.

| channel | arm | SNR50 | 95% CI | gain | decodes / 3000 |
|---|---|---|---|---|---|
| AWGN         | base  | -32.36 | [-32.40,-32.31] |  --   | 1597 |
| AWGN         | add3c | -32.46 | [-32.51,-32.42] | +0.10 | 1694 |
| AWGN         | coh3  | -32.52 | [-32.57,-32.47] | +0.16 | 1744 |
| AWGN         | ref0  | -32.53 | [-32.58,-32.49] | +0.17 | 1757 |
| AWGN         | add2  | -32.60 |                 | +0.24 | 1833 |
| AWGN         | add3  | -32.65 | [-32.70,-32.60] | +0.29 | 1864 |
| AWGN         | est3  | -32.05 |                 | -0.31 | 1337 |
| Rayleigh 0.1 | base  | -29.91 | [-29.96,-29.85] |  --   | 1631 |
| Rayleigh 0.1 | add3  | -29.92 | [-29.97,-29.86] | +0.01 | 1640 |
| Rayleigh 0.1 | coh3  | -30.01 | [-30.06,-29.95] | +0.10 | 1713 |
| Rayleigh 0.3 | base  | -28.38 | [-28.43,-28.33] |  --   | 1614 |
| Rayleigh 0.3 | add3  | -28.38 | [-28.43,-28.33] |  0.00 | 1614 |
| Rayleigh 0.3 | coh3  | -28.42 | [-28.48,-28.37] | +0.04 | 1652 |
| Rayleigh 1.0 | base  | -26.30 | [-26.35,-26.26] |  --   | 1549 |
| Rayleigh 1.0 | add3  | -26.31 | [-26.35,-26.26] | +0.01 | 1553 |
| Rayleigh 1.0 | coh3  | -26.29 | [-26.34,-26.24] | -0.01 | 1537 |

`add2` and `add3c` on the fading channels are the same story as `add3`: 0.00 dB
everywhere.  At Rayleigh 0.3 Hz the additive arms return *literally the same
1614 decodes* as the baseline, which is as clean a null as this harness can
produce.

### Where the 0.29 dB actually comes from

`ref0` is the control that matters.  It computes the coherent metric and then
refuses to let it open the gate, while leaving gate one (`minsync1`) permanently
open so that every candidate gets the fine lag and frequency search.

    gate one left open, coherent gate disabled  (ref0)   +0.17 dB
    coherent gate on, gate one left stock       (add3c)  +0.10 dB
    both                                        (add3)   +0.29 dB

The two effects are independent and add.  Two thirds of the headline number has
nothing to do with block-coherent anything: it is that `minsync1` is currently
throwing away candidates before they are properly refined, and refining them
raises their stock sync value enough to clear `minsync2`.

The coherent metric's own contribution is 0.10 dB on AWGN and zero on all three
fading channels.

### The coherent metric is a better scorer and a worse search objective

`est3` lets the coherent metric steer the drift refinement and the fine lag and
frequency search while the stock metric still decides the gate.  It **loses
0.31 dB**.  `coh3`, which does both, gains only 0.16 dB where the gate alone
would have given more.

The reason showed up in the calibration: maximised over its own search, the
coherent metric's false-alarm distribution has twice the spread of the stock
metric's (sd 0.0218 against 0.0109), because it has a quarter as many effective
degrees of freedom, so the search finds bigger noise peaks and drags the
frequency and lag estimate onto them.  Scoring at a point somebody else chose is
where it is worth something; choosing the point is where it is not.  That is
also why `-M 1` is both the cheap mode and the best one.

## 8. False decodes, at volume

WSPR has no CRC, so a false decode becomes a bad public spot.  Scored by
`tools/harness.py` unmodified: the transmitted messages are known, so anything
else the decoder prints is counted.

| corpus | arm | files | transmissions | decoded | false |
|---|---|---|---|---|---|
| busy sub-band, 25 signals | base | 1200 | 30000 | 27876 | 1 |
| busy sub-band, 25 signals | ref0 | 1200 | 30000 | 27999 | 1 |
| busy sub-band, 25 signals | add3 | 1200 | 30000 | 28053 | 1 |
| pure noise, generated at -60 dB | base | 4000 | 4000 | 0 | 0 |
| pure noise, generated at -60 dB | ref0 | 4000 | 4000 | 0 | 0 |
| pure noise, generated at -60 dB | add3 | 4000 | 4000 | 0 | 0 |

**The branch baseline is not zero on the busy sub-band at this volume**: one
false decode in 1200 files, `TU6DDR/7 13`.  Every arm produces the same single
false decode, from the same file, with the same text.  Neither `add3` nor `ref0`
introduces one, and `add3` recovers 177 more genuine spots on the same corpus.

Across the sensitivity campaign, 12000 files per arm over four channels, the
false-decode counts are also identical to the baseline: 0 on AWGN, 1 at
Rayleigh 0.1 Hz, 0 at 0.3 Hz, 2 at 1.0 Hz, for both `base` and `add3`, and the
same files in each case.  `coh3` produced two fewer, one at Rayleigh 0.1 Hz and
one at 1.0 Hz, and `est3` produced one that the baseline did not, at Rayleigh
0.3 Hz.

`add3c`, `coh3` and `add2` were queued for the busy sub-band test and cut short
to free a quiet machine for the CPU measurement.  `add3` differs from `add3c`
only by leaving gate one open, which refines more candidates and can only raise
their sync value, so in practice it admits everything `add3c` does and more; its
clean result therefore covers `add3c`, and `ref0` covers the gate-one change on
its own.  `coh3` and `add2` are untested on the busy sub-band and that is a gap.

## 9. CPU per file

hyperv-gha, Intel i7-13700 VM, 6 physical / 12 logical cores, one decoder
process at a time, `-P 1`, no `-Y`.  The machine was genuinely quiet: load
average 1.74 before, 1.17 after, all of it this measurement.

| corpus | base | add3c | ref0 | add3 | coh3 |
|---|---|---|---|---|---|
| single signal, -28 dB, decodes on pass 0 | 0.105 s | 0.107 | 0.108 | 0.108 | 0.087 |
| single signal, -32 dB, marginal          | 0.267 s | 0.270 | 0.326 | 0.340 | 0.233 |
| single signal, -36 dB, nothing there     | 0.021 s | 0.021 | 0.023 | 0.023 | 0.023 |
| busy sub-band, 25 signals per file       | 5.092 s | 5.212 | 5.514 | 5.607 | 4.100 |

`add3` costs +10% on the busy sub-band, which is the case that matters, and +27%
on a marginal single signal, which is where leaving gate one open makes the
decoder refine every candidate instead of a tenth of them.  `add3c`, which
leaves gate one alone, is free (+2%).  Most of `add3`'s cost is `ref0`'s cost,
i.e. the refinement, not the metric.

`coh3` is 19% CHEAPER than the baseline on a busy band.  Its gate is a better
discriminator in both directions: on the single-signal diagnostic it passes the
true signal slightly more often than the baseline does and still hands the
demodulator fewer candidates overall (0.30 against 0.36 per file on AWGN at
-33.5 dB, 1.13 against 1.80 at Rayleigh 1.0 Hz and -27.5 dB).

## 10. Verdict

The idea is sound, it does what the physics says it should, and it is not worth
shipping as a sensitivity change.

  * Block-coherent scoring is a genuinely better sync metric.  The gate passes
    47% more files at -33.5 dB on AWGN and the decode rate rises with it.
  * It buys 0.10 dB on a steady path and nothing measurable on Rayleigh at 0.1,
    0.3 or 1.0 Hz.  On a fading channel the sync gate is not the constraint, and
    no improvement to it can be.
  * The model predicted 1.5 dB at the gate.  The gap is not the metric failing;
    it is that only about half a dB of decodable signal is sitting behind the
    gate in the first place, and the demodulator cannot decode most of what is
    let through.  This is the same ceiling that limited the threshold-loosening
    experiment to 0.55 dB.
  * It is safe: no new false decodes in 1200 busy sub-band files, 4000
    pure-noise files and 12000 campaign files, and the calibration is at matched
    false-alarm rate rather than at a looser threshold.
  * The one finding worth carrying forward is the control.  Leaving `minsync1`
    open so every candidate gets the fine lag and frequency search is worth
    0.17 dB on AWGN by itself, needs no new metric, and produced no false
    decodes here.  It costs 8% CPU on a busy band.  It deserves its own
    experiment with its own safety testing rather than a footnote in this one --
    in particular on the fading and busy-band cases, where it also did nothing
    for sensitivity, and at higher busy-band volume than 1200 files.

Negative results, stated plainly:

  * Normalising by the ruled-out tones instead of by all four is algebraically
    the stock metric.  Not measured, because it is provable.
  * Using the coherent metric as the SEARCH objective costs 0.31 dB.
  * L=2 and L=3 are indistinguishable in practice (0.24 against 0.29 dB on
    AWGN, both nothing elsewhere), so the block length is not a knob worth
    turning.
  * Per-symbol weighting, the third idea in the brief, was not reached.  The
    diagnostic says why it would not have mattered: on the fading channels where
    weighting is supposed to help, the sync gate is already passing 89 to 100%
    of files and the loss is downstream.

## 11. Reproduction

    cd /home/tf/wspr-lab
    python3 tools/mkcsync.py            # regenerates src/wsprd-csync
    ./build.sh csync                    # -> bin/wsprd-csync

    python3 tools/csynctheory.py        # the mechanism in numpy, no decoder
    python3 tools/csyncident.py 5       # -K absent is byte-identical to the
                                        # branch; -P 1 == -P 8; additive mode
                                        # never loses a stock decode

    # the false-alarm distributions and the calibrated thresholds
    python3 tools/csyncroc.py --Ls 2,3 --noise 2000 --band 150 --single 250 \
        --jobs 5 --out results/csync/roc.json

    # the sensitivity campaign, split across two boxes on the 0.5 dB grid
    ARMS='add3:-K 3 -M 1 -T -1,0.1865,0.1435
          add2:-K 2 -M 1 -T -1,0.1702,0.1326
          add3c:-K 3 -M 1 -T 99,0.1865,0.1435
          coh3:-K 3 -M 0 -T -1,0.1961,0.1561
          est3:-K 3 -M 2 -T 0.1561,0.12,0.10
          ref0:-K 3 -M 1 -T -1,99,99'
    bash tools/csynccampaign.sh A 5 base: $ARMS      # whole dB,  this box
    ssh hyperv-gha 'bash tools/csynccampaign.sh B 5 base: '"$ARMS"   # half dB

    # the det / syn / dec breakdown, per arm
    bash tools/csyncdiag.sh base 5
    bash tools/csyncdiag.sh add3 5 -K 3 -M 1 -T -1,0.1865,0.1435

    # the acceptance test
    bash tools/csyncsafety.sh 5 "add3:-K 3 -M 1 -T -1,0.1865,0.1435"

    # CPU, on whichever box is quiet
    bash tools/csynccpu.sh "" "base:" "add3:-K 3 -M 1 -T -1,0.1865,0.1435"

    python3 tools/csynctable.py         # every table above
    python3 tools/fit_ci.py <json...>   # thresholds with confidence intervals

`-march=native` differs between the two boxes, so the source has to be synced
and rebuilt on the remote rather than the binary copied; a copied binary dies
with SIGILL.  Every run used `-P 1` inside the harness and at most 5 jobs per
box.

### Files

    tools/mkcsync.py       generates src/wsprd-csync from src/wsprd-stage
    tools/csynctheory.py   the metric in numpy: sections A to D
    tools/csyncident.py    off-by-default, determinism and additivity checks
    tools/csyncroc.py      false-alarm distributions and threshold calibration
    tools/csyncrun.py      sensitivity harness with the det/syn/dec counters
    tools/csynccampaign.sh the SNR campaign, one half of the grid per box
    tools/csyncdiag.sh     the det/syn/dec breakdown
    tools/csyncsafety.sh   false decodes at volume
    tools/csynccpu.sh      CPU per file, one process at a time
    tools/csynctable.py    merges the halves and prints every table
    tools/csynccal.py      an earlier, superseded calibration probe

    results/csync/*.json   every run, in tools/harness.py's format
    results/csync/TABLES.txt  the output of tools/csynctable.py
    results/csync/CI.txt      thresholds with 95% confidence intervals
    results/csync/theory.txt  the output of tools/csynctheory.py
    results/csync/roc.log     the calibration, including the false-alarm means
                              and standard deviations the thresholds came from
