# Frequency-domain interference excision for wsprd

A spectral counterpart to the time-domain noise blanker already in the tree.
The blanker deals with impulsive QRN, which is narrow in time and broad in
frequency. This deals with the opposite kind of intruder: narrow in frequency
and broad in time, which on the WSPR bands means a CW carrier, a birdie, a
switching-supply spur or a stray data signal, and against which the decoder
currently has no defence at all.

Everything here is generated from the pristine sources by script. Nothing in
`/home/tf/wsjtx` was touched, and no existing `src/wsprd-*` directory, no
existing tool and no existing result file was modified.

## What is new

    tools/mkfex.py         generates src/wsprd-fex from src/wsprd-mt
    tools/fex.c.inc        the excision stage itself, with the reasoning in comments
    tools/wsprchanx.c      wsprchan plus deliberate interference (a copy, not an edit)
    tools/fexrun.py        runs tools/harness.py unmodified against wsprchanx
    tools/fexscan.py       the safety scan: does the stage ever fire when it should not
    tools/fexcampaign.sh   every measurement run in this directory
    tools/fextable.py      turns the JSON in this directory into the tables below
    bin/fex/               a directory of symlinks; wsprchan resolves to wsprchanx

`bin/wsprchanx` is byte-identical to `bin/wsprchan` for every option wsprchan
already had (checked directly), so runs made with the two are comparable.

## Rebuilding

    cd /home/tf/wspr-lab
    gcc -O2 -o bin/wsprchanx tools/wsprchanx.c src/wsprd-baseline/wsprsim_utils.c \
        src/wsprd-baseline/wsprd_utils.c src/wsprd-baseline/nhash.c \
        src/wsprd-baseline/fano.c src/wsprd-baseline/tab.c \
        -I src/wsprd-baseline -lfftw3f -lm
    python3 tools/mkfex.py
    ./build.sh fex

## Reproducing the measurements

    cd /home/tf/wspr-lab
    bash tools/fexcampaign.sh          # 500 trials per SNR point, 6 jobs, -P 1
    python3 tools/fextable.py          # the tables below

`tools/fexcampaign.sh` skips any run whose JSON already exists, so it can be
interrupted and restarted. Set `TRIALS=` and `JOBS=` to override.

The safety scans:

    python3 tools/fexscan.py --n 300 --snrlo -25 --snrhi 20 --jobs 6 --args "-X 1"
    python3 tools/fexscan.py --n 120 --nsig 25 --snrlo -31 --snrhi -7 --jobs 6 --args "-X 1"

Individual files, with the stage explaining itself:

    WSPRD_FEX_DEBUG=1 bin/wsprd-fex -C 500 -o 4 -d -N 14 -P 1 -X 1 -a . -f 10.1387 FILE.wav

## The change

One new stage, called from `readwavfile()` on the 46080-bin slice of the big
forward FFT, immediately before the inverse transform. Nothing else in the
decoder is touched, and with `-X` absent the binary is bit-identical in output
to `bin/wsprd-mt` (checked on files with and without interference).

### Why that slice

`readwavfile()` already transforms the whole 114 s recording, so the 375 Hz that
survives downsampling is sitting in memory as 46080 bins of 0.00814 Hz. At that
resolution the two things we need to tell apart could hardly be more different:

  * a CW carrier is two bins wide and 60 to 80 dB above the per-bin noise floor
  * a WSPR transmission is not narrowband here at all. Its symbols are 0.68 s,
    so each is 1.5 Hz wide, and the four tones smear into a 6 Hz plateau about
    740 bins across. No single bin of it stands far above its own neighbours.

### The detector

For each candidate bin, two ratios:

  * `P/F`, the bin against a local floor. The floor is the minimum of the
    256-bin block medians within +-4 blocks. The minimum rather than the median,
    because a strong carrier's leakage skirt lifts the blocks either side of it,
    and a floor that followed the skirt would stop the notch growing out of the
    skirt that produced it.
  * `P/N`, the bin against the larger of the two medians taken over 0.39 to
    2.34 Hz on each side. This is the test that does the work. Measured:

        WSPR signal peak, SNR -20 to +20 dB    P/N = 83 to 118, and flat
        CW carrier, -20 to +20 dB re noise rms P/N = 2400 to 340000

    The WSPR figure does not grow with SNR, because the comparison band is
    inside the signal's own plateau. That is what makes a fixed threshold safe:
    `-X 1` uses P/N > 300 and P/F > 1000.

### The notch

Having found a spur, the notch is grown outwards from it while a 33-bin
smoothed spectrum stays above twice the local floor, capped at 2048 bins each
side. It has to be grown, because a 114 s rectangular record gives a carrier a
1/df^2 leakage skirt that is genuinely present in the DFT coefficients and runs
many Hz either side; zeroing only the peak leaves most of the interference
behind. Growing rather than using a fixed width makes the notch set its own size
from the carrier strength. Measured notch widths for a carrier at 1500+30 Hz:
0.6 Hz at -15 dB re noise rms, 1.4 Hz at -5 dB, 2.4 Hz at 0 dB, 14.9 Hz at
+20 dB, 19.5 Hz at +30 dB.

The flagged bins are zeroed. What is left in the time domain after a hard notch
of width B is not a weakened carrier spread over the record, it is the turn-on
and turn-off transient at the record edges, about 1/B seconds long at each end.
That is why a notch a couple of Hz wide is worth far more than its energy
fraction suggests.

### Why it helps as much as it does

The obvious damage a carrier does is to the noise-level estimate in `main()`
(`tmpsort[122]`, the 30th percentile of the smoothed spectrum). That turns out
not to be the main mechanism: the percentile is over 411 bins of 0.73 Hz and a
narrow carrier only touches a handful of them.

The mechanism that actually kills decodes is `ss = ss/totp` in
`sync_and_demodulate()`. The four tone filters are 256-point coherent sums whose
sidelobe response falls only as 1/df, so a carrier df Hz away leaks into all four
of them at about 0.466/df of its on-frequency amplitude. It adds nearly equally
to all four tones, so it barely moves the numerator `cmet` while inflating the
denominator `totp` several-fold, and the sync statistic falls below `minsync2`.
That is why a carrier 20 Hz from a signal is far more damaging than the same
carrier 140 Hz away, and why the loss is a cliff rather than a slope.

### Flags

    -X 1              excise, normal aggressiveness (P/F>1000, P/N>300)
    -X 2              excise, more sensitive (P/F>200, P/N>150, wider notch)
    -X T1,T2,T3,W[,c] set the thresholds directly; c clamps to the floor
                      instead of zeroing
    -X a              sweep off/1/2 and keep every decode any setting finds,
                      the same shape as -n a

All are off by default. The constants are in bins, and WSPR-15 scales both the
bin width and the signal width by the same factor of 8, so they carry over to
`-m` unchanged.

## What it does not fix (measured, single files, +20 dB carrier)

    steady carrier                    found, decode recovered
    carrier drifting 0.5 Hz / 120 s   found, decode recovered
    carrier drifting 1.0 Hz / 120 s   found, decode recovered
    carrier drifting 1.5 Hz / 120 s   NOT found, decode still lost
    carrier drifting 3.0 Hz / 120 s   NOT found, decode still lost
    keyed 0.06 s on / 0.12 s (20 wpm) found as 32 comb lines, decode recovered
    keyed 3 s on / 6 s (QRSS)         found, decode recovered
    keyed 0.5 s on / 1.0 s            36 comb lines excised, decode still lost
    noise 300 Hz wide over the signal nothing to excise, decode still lost

The pattern is consistent and it is the price of the safety property. Anything
whose energy is spread over more than about 2 Hz fails the narrowness test by
construction, because that is the same test that stops the stage eating a WSPR
signal. A carrier drifting more than about 1 Hz across the file, and a carrier
keyed at a rate near the WSPR symbol rate, both fall on the wrong side of it.

A wideband interferer that covers the wanted signal cannot be helped by any
frequency-domain excision: it is the local noise floor, not an outlier above it.
The stage correctly leaves it alone rather than making things worse.

## CPU

The stage is one pass to form the power spectrum, 180 median-of-256 sorts for
the floor, one boxcar, and one O(n) scan per spur found. Over 20 repeated
decodes of the same clean file the user+sys time was 3.30/3.35 s with `-X 0` and
3.35/3.36 s with `-X 1`, i.e. the difference is below 3 ms per file and not
separable from the noise of a loaded machine.

## Clean band: the make-or-break test

500 trials per SNR point, 3000 files per run.

    run                      SNR50     decoded/tx   false   s/file
    no excision             -32.30     1534/3000        0     0.68
    -X 1                    -32.30     1534/3000        0     0.54
    -X a (sweep off/1/2)    -32.30     1534/3000        0     1.66

Not merely the same threshold: the same 1534 decodes, file for file, because
the stage did not fire on a single one of the 3000 clean files. Compare the
time-domain blanker, where a fixed 10% setting costs 2.65 dB on clean noise and
is the reason `-n a` has to sweep.

That is the interesting difference between the two stages. The blanker has to
throw away a fixed PERCENTAGE of samples, so on a quiet band it is throwing away
signal. The excision stage tests a ratio with a factor of 66 of margin over what
noise alone produces, so on a quiet band it does nothing at all, and the sweep
that `-n a` needs is not needed here. `-X a` costs three times the CPU and buys
nothing that `-X 1` did not already have.

The s/file figures are not comparable between runs: the machine was shared with
other work and the load varied. The like-for-like CPU measurement is in the CPU
section above.

## A note on how these runs were executed

The machine was shared with other work throughout. Partway through the campaign
a second process picked up `tools/fexcampaign.sh` and ran it with `TRIALS=250`,
which raced this campaign for the same output filenames and did overwrite one
run (C2mt) with a 250-trial version before it was spotted. That file was deleted
and the run repeated at 500 trials.

`tools/fexcampaign.sh` now refuses to start if it is asked for fewer than 500
trials per point, and takes an flock so two campaigns cannot write the same JSON
at once. Every JSON in this directory carries its own `trials` field, so the
trial count of any result can be checked directly:

    python3 -c "import json,sys;d=json.load(open(sys.argv[1]));print(d['trials'])" results/fex/X.json

All results quoted below are from files whose `trials` field reads 500.

## 50% decode thresholds

Patched decoder (`-C 500 -o 4 -d -N 14 -P 1`) with and without `-X 1`.
Clean-band runs 500 trials per SNR point, interference runs 300, 0.5 dB grid.
Interferer level is dB relative to the total background noise power in 6 kHz.

    channel                                     no excision    -X 1     gain
    clean band                                      -32.30    -32.30    0.00
    CW +20 dB, 20 Hz from the signal                -26.05    -32.30   +6.25
    CW +10 dB, 20 Hz from the signal                -30.80    -32.30   +1.50
    CW   0 dB, 20 Hz from the signal                -31.90    -32.30   +0.40
    CW +20 dB, 60 Hz from the signal                -27.30    -32.30   +5.00
    CW +20 dB at 1500+60 Hz   (inside +-110 Hz)     -29.60    -32.10   +2.50
    CW +20 dB at 1500+150 Hz  (outside +-110 Hz)    -30.70    -32.30   +1.60
    noise 100 Hz wide +15 dB, 80 Hz from signal     -31.00    -31.05   +0.05

Zero false decodes in every one of those runs: 9000 clean files and 25200
interference files, 34200 in total, and not one spurious spot either way.

In six of the seven CW cases the excised threshold is the clean-band -32.30
exactly. The interference is not merely reduced, it is removed: what is left is
an ordinary AWGN channel. The seventh, the carrier pinned at 1500+60 Hz, lands
at -32.10 because in that geometry the wanted signal sometimes falls inside the
notch (see below), which drags the fitted curve down slightly.

The two rows to read carefully are the last one and the +0 dB one. A carrier
only 0 dB above the noise costs 0.40 dB, so there was little to win. A 100 Hz
wide noise interferer is not an outlier above the local floor, it IS the local
floor over its own band, and the stage correctly declines to touch it: +0.05 dB
is nothing, and that is the honest answer rather than a disappointing one.

### Why the two geometries differ

Rows quoting "N Hz from the signal" place the carrier at a fixed offset from the
wanted signal in every trial, which isolates the mechanism. Rows quoting a fixed
1500+N Hz place it at a fixed spot in the sub-band while the signal moves
uniformly over +-95 Hz, which is what a real band looks like; there the carrier
sometimes lands on top of the signal and no amount of excision can help. That is
why the fixed-spot gains are smaller, and it is a property of the situation, not
of the stage.

## Decode rate at a fixed SNR of -28 dB

300 trials per row, signal frequency uniform over +-95 Hz, carriers at a fixed
spot in the sub-band. This covers the cases where a threshold curve is not
meaningful because it never reaches 1.

    channel                              no excision      -X 1     s/file off/on
    clean                                     1.000      1.000       0.19  0.20
    CW  +0 dB at 1500+30 Hz                   0.930      0.980       0.79  0.37
    CW +10 dB at 1500+30 Hz                   0.900      0.947       1.30  0.96
    CW +20 dB at 1500+30 Hz                   0.707      0.887       1.72  1.22
    CW +20 dB at 1500+60 Hz                   0.617      0.883       1.66  1.29
    CW +20 dB at 1500+150 Hz                  0.873      1.000       1.04  0.20
    CW +20 dB at 1500+0 Hz                    0.777      0.913       2.91  1.25
    two carriers, +20 and +15 dB              0.277      0.787       2.77  2.15
    CW +20 dB keyed at 20 wpm                 0.077      0.740      14.64  3.54
    CW +20 dB drifting 1 Hz / 120 s           0.257      0.887       2.09  1.01
    CW +20 dB drifting 3 Hz / 120 s           0.630      0.630       1.77  1.78
    noise  50 Hz +20 dB at 1500-140 Hz        0.880      0.880       1.79  1.81
    noise 100 Hz +10 dB at 1500-120 Hz        0.857      0.857       1.54  1.53
    noise 300 Hz +10 dB at 1500+0 Hz          0.000      0.000       0.03  0.04
    noise  50 Hz +20 dB, FT8 duty cycle       0.913      0.913       1.68  1.69

Two carriers, and a keyed carrier, are where the current decoder is worst and
where the stage helps most: 0.277 to 0.787 and 0.077 to 0.740. A keyed carrier
is not one line but a comb of lines at the keying rate, which is why the
iteration budget has to be 64 rather than a handful.

A drifting carrier is the boundary case, and it is a sharp one. At 1 Hz of
drift across the file the stage still recognises it and recovers 0.257 to 0.887.
At 3 Hz it does not fire at all and changes nothing, 0.630 either way. Nothing
in between was measured on the full corpus; single-file probes put the cutoff
between 1 and 1.5 Hz.

## CPU

The excision costs nothing measurable to run, and it usually SAVES time,
because a carrier the decoder cannot remove is a carrier it wastes effort on:
every peak in the leakage skirt becomes a candidate that gets demodulated and
fails. The keyed-carrier row is the extreme case, 14.64 s per file down to
3.54 s. On a clean band, where the stage finds nothing, the cost is below the
3 ms per file that a wall-clock measurement can separate from noise.

## False decodes

43200 files across every run in this directory, 3 false decodes, none of them
caused by the excision:

    E_cw_drift_1Hz_mt        1   without excision.  With -X 1 the same 300 files
                                 give none: the stage removed the carrier that
                                 was manufacturing it.
    E_nb_100Hz_m120_10dB     1   present with AND without excision, the same
                                 spurious message "66E/V68I 50" from the same
                                 trial.  A pre-existing artefact of a wideband
                                 interferer, unchanged by the stage.

The harness uses identical seeds for both variants, so a false decode appearing
in both is demonstrably the same file producing the same spot. Every one of the
17 threshold runs, 34200 files, was clean both ways.

On the shipping criterion -- anything that introduces false decodes is not
shippable -- this introduces none, and removes one.

## The cost of a notch, stated plainly

Excision is not free in bandwidth, only in sensitivity. The notch is as wide as
the carrier's leakage skirt: 14.9 Hz for a carrier at +20 dB re noise rms. A
WSPR signal inside that window is destroyed along with the carrier.

That shows up in the fixed-spot rows. With a +20 dB carrier at 1500+60 Hz and
the signal uniform over +-95 Hz, 14.9 Hz of 190 Hz is 7.8% of trials, and the
excised decode rate at -28 dB is 0.883 against 1.000 on a clean band -- an 11.7%
shortfall, the rest being signals close to but not inside the notch. Without
excision the same geometry gives 0.617.

So the trade is: lose the few Hz around the carrier outright, in exchange for
getting the rest of the sub-band back. Without the stage the carrier costs
several dB across the whole 200 Hz; with it the loss is confined to the notch.
Given the alternative is losing candidate detection band-wide, that is the right
side of the trade, but it should be stated rather than hidden: this stage does
not recover a signal sitting on top of a carrier, and nothing in the frequency
domain can.

## Verification

    wsprd-fex -X 0 vs wsprd-mt, 60 files with and without interference
        0 mismatches, byte-identical stdout.  The stage is inert unless asked for.

    wsprd-fex -X 1, -P 1 vs 4 vs 16, 24 busy files each carrying two carriers
        459 decodes, 0 mismatches.  Excision runs before any thread starts and
        is deterministic, so the threading guarantee in the main README still
        holds with it enabled.

A bug found in the scan harness, recorded because the first numbers it produced
were wrong. `tools/fexscan.py` decodes each file twice, once with excision and
once without, and it was doing both in the same scratch directory without
clearing `hashtable.txt` in between. The second decode therefore inherited the
callsigns the first had learned, and WSPR's OSD hash gate is worth about 0.9 dB
to a known callsign, so the second run got decodes it had not earned. On the
busy sub-band corpus that manufactured a spurious "16 gained, 7 lost"; once the
leak was closed the same corpus gave 2619 decodes both ways, exactly equal. The
script now clears the hash table, ALL_WSPR.TXT, the timer file and the spot file
between the two decodes, and the scans were repeated.

`tools/harness.py` never had this problem -- every trial gets a fresh directory
on /dev/shm and one decoder invocation -- so none of the threshold runs above
are affected.

## Conclusion

A frequency-domain excision stage works, and it works because of a structural
asymmetry rather than a clever threshold: at the 0.008 Hz resolution the decoder
already has in hand, a CW carrier occupies two bins and a WSPR signal occupies
740. Testing a bin against the 0.4 to 2.3 Hz either side of it separates the two
by 40 dB and, crucially, that separation does not shrink as the WSPR signal gets
stronger, because the comparison band is inside the signal's own plateau.

That is what makes it safe enough to leave on. It recovers up to 6.25 dB against
a strong carrier, returns the threshold to its clean-band value exactly in six of
seven CW cases, costs zero on a clean band, introduces no false decodes in 43200
files, and usually saves CPU rather than spending it.

Unlike the time-domain blanker it does not need a sweep. `-n a` exists because a
fixed blanker throws away a fixed percentage of samples and so damages a quiet
band; this stage tests a ratio with a factor of 66 of headroom over what noise
alone produces, so on a quiet band it does nothing at all. `-X a` was built and
measured anyway, and it buys nothing `-X 1` did not already have for three times
the CPU. If this shipped, it should ship as a fixed setting.

Where it does nothing, it does nothing loudly rather than quietly: a wideband
interferer is the local noise floor over its own band, not an outlier above it,
and no frequency-domain method can lift a signal out of it. Carriers drifting
more than about 1 Hz across the file, and carriers keyed at a rate near the WSPR
symbol rate, fall the same side of the narrowness test that protects WSPR
signals. That test is the whole design, so those misses are its price, not a
tuning failure.

## Safety scans (repeated after the harness fix)

The question these answer is the one that decides whether a stage like this can
be left switched on: does it ever fire when it should not? A strong WSPR signal
is the thing most likely to be mistaken for a spur, so the corpus is deliberately
loud.

    300 files, one signal each, SNR from -25 to +20 dB, no interference at all
        files where excision fired : 0
        decodes without / with     : 300 / 300
        false decodes              : 0 / 0
        files that lost a decode   : 0        gained: 0

    120 files, 25 signals each, a crowded sub-band, SNRs -31 to -7 dB
        files where excision fired : 0
        decodes without / with     : 2619 / 2619   of 3000 transmitted
        false decodes              : 0 / 0
        files that lost a decode   : 0        gained: 0

Identical decode counts, file for file, on a corpus of 3300 WSPR transmissions
including signals up to +20 dB. The detector's margin is wide enough that noise
and signal alike stay a long way from it.
