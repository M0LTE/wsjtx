# Adaptive coherent demodulation: measuring the channel instead of guessing

The coherent demodulator in `tools/coherent.c.inc` tracks the channel by smoothing decision-directed complex amplitudes over `nsmooth` symbols, leaving each symbol out of its own estimate. The right window is set by the channel's coherence time, which nobody knows in advance, so the shipping decoder tries eight windows (162, 81, 41, 27, 15, 9, 5, 3 symbols, `-N 7` through `-N 14`) and keeps the first that decodes. This build measures the coherence time and uses it.

Everything here was built and measured on `hyperv-gha`. Binaries are not portable between the two machines (`build.sh` compiles `-march=native`), so rebuild before running.

## What was built

    tools/mkdop.py       generates src/wsprd-dop from src/wsprd-mt
    ./build.sh dop       -> bin/wsprd-dop
    tools/dopcheck.py    estimator vs the simulator's known Doppler spread
    tools/rundop.sh      the threshold sweep
    tools/rundopfix.sh   the same, one fixed window at a time
    tools/doptable.py    prints the threshold table from the JSON

New decoder option:

    -A m   adaptive coherent demodulation, replacing the eight fixed windows
           1  one trial, rectangular window of the measured matched length
           2  one trial, Wiener weights matched to the measured spectrum
           3  three trials: Wiener matched, then 4x and 1/4 the measured spread
           4  three trials: rectangular matched, then 4x and 1/4
           7  one trial, Wiener weights on a decision-free pilot channel estimate

`-A 0` (the default) is the unmodified decoder; `wsprd-dop` without `-A` produces
byte-identical output to `wsprd-mt`.

Environment variables, for measurement only:

    WSPRD_DOPLOG=1   log every coherence estimate to stderr as a DOPEST line
    WSPRD_DOPLOG=2   also dump the raw autocorrelation as a DOPR line
    WSPRD_COHONLY=1  drop the six noncoherent trials, leaving only the adaptive one
    WSPRD_FIXW=n     replace the coherent trials with one fixed window of n symbols

## How the estimate is made

### A pilot that costs nothing

WSPR's sync vector fixes which pair of tones a symbol can use: `pr3[i]` selects
tones {0,2} or {1,3}, and the data bit picks one of the two. Those two tones are
exactly orthogonal over one symbol, and the symbol's starting phase is
data-independent (that is the property the whole coherent demodulator rests on).
So the sum of the two candidate tone correlations,

    u_i = z_t0[i] + z_t1[i]

carries the signal at full amplitude with its true phase, plus the noise from
both tones and nothing else. No decision is involved, so no decision error can
bias it. It is 3 dB noisier than a correct decision-directed estimate and
infinitely better than a wrong one.

Its autocorrelation is

    R(0) = A^2 + 2N        R(m) = A^2 * rho(m),  m >= 1

with the noise independent from symbol to symbol because the matched filter is
exactly one symbol long.

### The model rho is fitted to

rho is not the physical channel's autocorrelation: it is the channel as the
demodulator sees it, which is the physical channel averaged over one symbol by
that same matched filter. For a Gaussian Doppler spectrum of standard deviation
sigma, writing beta = 2 pi^2 sigma^2 T^2 with T the symbol period,

    rho_beta(m) = INT_-1^1 (1-|x|) exp(-beta (m+x)^2) dx,  normalised to rho(0)=1

the triangular kernel being the autocorrelation of the one-symbol boxcar. That
correction is not cosmetic. At a 1 Hz spread the physical channel has fallen to
0.10 by the next symbol but the observed one is still at 0.28, and a model that
ignores the averaging reads every fast channel as far faster than it is.

rho_beta does not depend on the data, so the whole family is tabulated once at
start-up over 80 values of beta and the fit is a scan: for each beta, fit the
scale A^2 by weighted least squares, keep the beta with the smallest residual.
R(0) minus the fitted A^2 gives the noise and hence the per-symbol SNR. The ITU
2-sigma spread the simulator's `-R` takes is 2*sigma.

### Two things had to be taken out before the fit worked

**A spurious coherent floor.** The demodulator fits its own residual carrier and
drift by maximising |sum w_i exp(-j psi_i)| over a grid of some 10^5 hypotheses
and de-rotates by the winner. Taking the maximum of that many projections leaves
the sequence with a coherent component aligned with the noise rather than the
channel, and in the autocorrelation it appears as a constant floor at every lag.
Measured on the channel where the true correlation past lag 2 is known to be zero
(1 Hz spread), the floor is 7/162 of the total power. A fit without that term has
only one knob and is forced to explain the floor by declaring the channel stable,
which on a fast path is exactly backwards. Subtracting `DOP_KAPPA*R(0)/162`
before fitting fixes it. A path with a genuine specular component would leave a
floor of its own and read as slower than its diffuse part really is; the
simulator's fading is pure Rayleigh, so that case is untested.

**A tail that is noise.** 110 seconds holds only about ten independent fades even
at 0.1 Hz, so past the coherence time the sample autocorrelation of one
realisation is a random walk about zero with excursions as large as a third of the
peak. Fitting it with equal weight is fitting noise, and it was what made the
first working version collapse to "stable" on every fading channel. The fit range
now follows the decay: find the lag where the correlation has fallen to a fifth,
and fit out to twice that. On a stable path nothing falls and all 60 lags are used.

### From rho to the smoother

The demodulator normalises its channel estimate to unit modulus and gets its
maximal-ratio weighting from the data term instead, so what matters is the
*direction* of the estimate, not its scale. The figure of merit is therefore the
correlation Cov(hhat, h)/sqrt(Var(hhat)), not the mean-square error.

  * `-A 1` maximises that over rectangular leave-one-out windows and returns
    the best `nsmooth`, not restricted to the eight-value ladder.
  * `-A 2` maximises it over all weightings, which is the leave-one-out Wiener
    solution: solve (rho + I/gamma) g = rho over the taps either side of the
    symbol, excluding the symbol itself. The taps run out to where rho has
    fallen to exp(-9). One Cholesky of at most 162x162 per call, which is about
    2% of what the demodulator already spends on its carrier-fit FFTs.


## 1. Does the estimator work?

Against the simulator's `-R`, which is the ITU 2-sigma Doppler spread, 150 files
per cell, single signal pinned at 0 Hz offset and 0 DT, taking the estimate
belonging to the candidate closest to 0 Hz. `tools/dopcheck.py`, raw data in
`results/dop/est_final.json`.

    true spread    reported fd, p25 / p50 / p75, at three SNRs in 2500 Hz
      (Hz)            -20 dB               -26 dB               -29 dB
      0.00      0.000 0.000 0.001    0.000 0.000 0.002    0.000 0.000 0.002
      0.03      0.007 0.013 0.027    0.006 0.010 0.022    0.007 0.014 0.034
      0.05      0.008 0.022 0.052    0.010 0.034 0.058    0.013 0.037 0.065
      0.10      0.052 0.090 0.125    0.052 0.100 0.139    0.058 0.112 0.139
      0.20      0.193 0.215 0.268    0.155 0.215 0.268    0.173 0.240 0.333
      0.30      0.268 0.333 0.371    0.268 0.371 0.462    0.268 0.333 0.462
      0.50      0.414 0.515 0.641    0.414 0.575 1.234    0.298 0.575 3.297
      1.00      0.298 0.575 3.297    0.030 0.414 3.297    0.000 0.333 3.297

Yes, over the range that matters, and the answer barely moves between -20 dB and
-29 dB, which is the property that makes it usable: it is measuring the channel,
not the SNR.

  * **0.1 to 0.5 Hz: accurate to about 20%.** Median 0.10 at 0.10, 0.22 at 0.20,
    0.35 at 0.30, 0.55 at 0.50.
  * **A stable path reads as stable.** Median exactly 0.000, and the window it
    picks is the full 163 symbols, which is the right answer.
  * **Below 0.05 Hz it reads low**, by about a factor of two at 0.05 and three at
    0.03. Not a bug in the fit: 110 seconds at 0.03 Hz contains about three
    independent fades, so the realisation genuinely does not decorrelate within
    the transmission. It errs towards more smoothing, which is the safe direction.
  * **At 1 Hz it falls apart.** The observed channel is down to 0.28 by lag 1 and
    to nothing by lag 2, so exactly one lag carries information and its
    measurement SNR is about 2. The quartiles run from 0.03 to the top of the
    grid. It still picks a 3-symbol window most of the time, which is the right
    behaviour, but the number it reports is not trustworthy.
  * **Above 2 Hz there is nothing to measure**, because the sync search stops
    finding the signal at all: 12 detections out of 150 at -20 dB.

The per-symbol SNR the fit returns is biased low by roughly a factor of two
(8.6 measured against 17 true at -20 dB on a stable path). That is a sensitive
ratio near the top -- 9% of error in the signal-power split halves it -- and its
only effect is to bias the chosen window slightly longer.

## 2. Thresholds, adaptive against brute force

500 trials per SNR point, 0.5 dB grid, 6 or 7 points per channel, fresh random
callsign every trial so the hash table is never primed. `tools/rundop.sh`,
`tools/doptable.py`.

    channel                variant                    SNR50   vs bf   files  false  s/file
    AWGN, stable path      stock -N 4                -31.10   +1.20    3000      0   0.225
    AWGN, stable path      brute force -N 14         -32.30   +0.00    3000      0   0.384
    AWGN, stable path      adaptive rect -A 1        -32.25   +0.05    3000      0   0.308
    AWGN, stable path      adaptive Wiener -A 2      -32.25   +0.05    3000      0   0.310
    AWGN, stable path      adaptive Wiener x3 -A 3   -32.25   +0.05    3000      0   0.341

    Rayleigh 0.05 Hz       stock -N 4                -29.50   +1.05    3500      0   0.421
    Rayleigh 0.05 Hz       brute force -N 14         -30.55   +0.00    3500      1   0.931
    Rayleigh 0.05 Hz       adaptive rect -A 1        -30.15   +0.40    3500      1   0.478
    Rayleigh 0.05 Hz       adaptive Wiener -A 2      -30.25   +0.30    3500      0   0.480
    Rayleigh 0.05 Hz       adaptive Wiener x3 -A 3   -30.40   +0.15    3500      0   0.610

    Rayleigh 0.10 Hz       stock -N 4                -29.20   +0.70    3500      0   0.573
    Rayleigh 0.10 Hz       brute force -N 14         -29.90   +0.00    3500      0   1.847
    Rayleigh 0.10 Hz       adaptive rect -A 1        -29.70   +0.20    3500      0   0.876
    Rayleigh 0.10 Hz       adaptive Wiener -A 2      -29.75   +0.15    3500      0   0.885
    Rayleigh 0.10 Hz       adaptive Wiener x3 -A 3   -29.90   +0.00    3500      0   1.124

    Rayleigh 0.30 Hz       stock -N 4                -28.30   +0.10    3500      1   0.619
    Rayleigh 0.30 Hz       brute force -N 14         -28.40   +0.00    3500      1   1.786
    Rayleigh 0.30 Hz       adaptive rect -A 1        -28.40   +0.00    3500      1   0.859
    Rayleigh 0.30 Hz       adaptive Wiener -A 2      -28.40   +0.00    3500      1   0.864
    Rayleigh 0.30 Hz       adaptive Wiener x3 -A 3   -28.40   +0.00    3500      1   1.193

    Rayleigh 1.00 Hz       stock -N 4                -26.35   +0.00    3000      0   0.399
    Rayleigh 1.00 Hz       brute force -N 14         -26.35   +0.00    3000      0   1.182
    Rayleigh 1.00 Hz       adaptive rect -A 1        -26.35   +0.00    3000      0   0.585
    Rayleigh 1.00 Hz       adaptive Wiener -A 2      -26.35   +0.00    3000      0   0.592
    Rayleigh 1.00 Hz       adaptive Wiener x3 -A 3   -26.35   +0.00    3000      0   0.752

The AWGN and 0.1 Hz brute-force figures reproduce SUMMARY.txt (-32.30 and -29.85,
here -29.90) on an independent corpus, so the two sets of numbers are comparable.

**The brute force wins, by 0.05 dB on AWGN and 0.15 to 0.30 dB on the slow
fading channels where there is anything to win.** At 0.3 Hz and above nothing
separates any of them, and at 0.3 Hz coherent demodulation of any kind is worth
only 0.1 dB over stock.

## 3. Is it the estimator that is short, or the extra rolls of the dice?

The brute force gets eight attempts at the codeword and the adaptive version
gets one. Each attempt is another draw against the Fano and OSD decoders, so
some of the gap is nothing to do with picking the right window. To separate the
two, every fixed window was measured on its own, with the same trial budget the
adaptive version spends: six noncoherent trials plus one coherent window.
`tools/rundopfix.sh`, `tools/dopfixtable.py`, 500 trials per point.

    channel                 w162     w41     w27     w15      w9      w5      w3   best fixed
    AWGN, stable          -32.25  -32.25  -32.20  -32.15  -32.10  -31.80  -31.50   w162 -32.25
    Rayleigh 0.05 Hz      -29.85  -29.85  -29.95  -30.30  -30.35  -30.15  -29.85     w9 -30.35
    Rayleigh 0.10 Hz      -29.40  -29.40  -29.40  -29.45  -29.75  -29.75  -29.50     w9 -29.75
    Rayleigh 0.30 Hz      -28.25  -28.25  -28.25  -28.25  -28.25  -28.30  -28.35     w3 -28.35
    Rayleigh 1.00 Hz      -26.35  -26.35  -26.35  -26.35  -26.35  -26.35  -26.35    any -26.35

Putting the three together:

    channel        best single window   adaptive single (-A 2)   brute force (8 windows)
    AWGN                -32.25                 -32.25                    -32.30
    0.05 Hz             -30.35                 -30.25                    -30.55
    0.10 Hz             -29.75                 -29.75                    -29.90
    0.30 Hz             -28.35                 -28.40                    -28.40
    1.00 Hz             -26.35                 -26.35                    -26.35

**The estimator does its job.** On four of the five channels one adaptive
demodulation is worth as much as the best fixed window chosen with hindsight,
and on 0.05 Hz it is 0.10 dB short. The rest of the brute force's margin --
0.05 dB on AWGN, 0.20 at 0.05 Hz, 0.15 at 0.10 Hz -- is not window selection at
all. It is having eight independent attempts at a code that is being run at its
limit, and no amount of estimating buys that back. Trying three windows around
the estimate (`-A 3`) buys most of it back for half the extra work: it ties the
brute force at 0.10 Hz and beats the best hindsight-chosen fixed window at
0.05 Hz.

## 4. CPU per file

Wall seconds per file at `-P 1`, averaged over each channel's SNR grid, from the
same runs as the threshold table.

    channel        stock   brute force   -A 1   -A 2   -A 3      -A 2 vs bf
    AWGN           0.225      0.384      0.308  0.310  0.341        0.81x
    0.05 Hz        0.421      0.931      0.478  0.480  0.610        0.52x
    0.10 Hz        0.573      1.847      0.876  0.885  1.124        0.48x
    0.30 Hz        0.619      1.786      0.859  0.864  1.193        0.48x
    1.00 Hz        0.399      1.182      0.585  0.592  0.752        0.50x
    all            0.456      1.253      0.632  0.637  0.820        0.51x

One adaptive demodulation instead of eight fixed ones halves the decoder's time
on any channel where the coherent trials actually run to completion. On AWGN the
saving is only 19%, because there the brute force decodes on its first coherent
trial and never reaches the other seven. The three-window version costs 0.65x.

The estimator itself is free. It is one 60-lag autocorrelation and a scan over
80 tabulated shapes, a few thousand operations, against the roughly 30 Mflops
the demodulator already spends on the carrier-fit FFTs. The Wiener design adds
one Cholesky of at most 162x162, about 2% on top, which is why `-A 1` and `-A 2`
have the same cost to three digits.

## 5. False decodes

16500 files per variant.

    stock -N 4                1
    brute force -N 14         2
    adaptive rect -A 1        2
    adaptive Wiener -A 2      1
    adaptive Wiener x3 -A 3   1

Every one of them is the same message in the same channel across every variant
that saw it -- `92TNP DQ97 13` at 0.3 Hz appears in all five including stock
`-N 4`, and `0QB/CZ2DCB 60` at 0.05 Hz appears in both the brute force and
`-A 1`. They come from the noncoherent trials that all the variants share.
**No false decode in this corpus is attributable to the coherent demodulator,
adaptive or brute force.** The rate is about one per 8000 files either way.

## 6. Correctness checks

    wsprd-dop with no -A, 20 files, vs wsprd-mt -N 14      identical output
    -A 2 and -A 3, 30 files, -P 1 vs -P 8                  identical output

The adaptive path is a pure function of the sample data, so the speculative
parallel trial grid stays deterministic.

## Verdict

The coherence time is measurable and the estimator measures it: from 0.1 to
0.5 Hz it lands within about 20% of the simulator's ground truth, it reads a
stable path as stable, and its answer is the same at -29 dB as at -20 dB. Fed
into a matched smoother, one demodulation is worth as much as the best fixed
window picked with hindsight on four of five channels.

**It still loses to the brute force**, by 0.05 dB on AWGN and 0.15 to 0.30 dB on
the slow fading channels. The reason is not the estimate. It is that eight
demodulations give the Fano and OSD decoders eight attempts at a 50-bit message
with no CRC, and at the code's limit those extra attempts are worth 0.05 to
0.20 dB on their own. That is the same finding as the GLRT experiment in
SUMMARY.txt from the other direction: near the limit, extra search helps and
extra cleverness does not.

**So the value here is CPU, not sensitivity.** `-A 2` costs half the brute
force's time for 0.15 to 0.30 dB. `-A 3` costs 0.65x for 0.00 to 0.15 dB, ties
the brute force on AWGN to within 0.05 dB and at 0.10 Hz exactly, and would be
the setting to ship if the eight-window search were ever the thing that had to
give. On the target the task set -- the 0.44 dB of headroom FST4W-120 has on the
0.1 Hz fading channel -- this closes none of it. The remaining headroom is not
in the smoothing window.

Two things that did not pan out, measured and rejected:

  * **A decision-free channel estimate** (`-A 7`): use the two-tone pilot sum as
    the channel input as well as for the estimator, so a wrong decision cannot
    poison the channel estimate at all. It is 3 dB noisier per symbol but immune
    to error propagation. Screened at 0.1 Hz on its own 750-file corpus
    (`results/dop/scr_f01_A*.json`, 250 trials at each of -29, -30, -31) it
    measured -29.65 against -29.70 for the decision-directed Wiener on the same
    corpus, so the 3 dB of extra noise costs slightly more than the decision
    errors do.
  * **Wiener weighting over a rectangular window** is real but small: 0.05 to
    0.10 dB at 0.05 and 0.10 Hz, nothing elsewhere. Worth having since it is
    free, not worth having on its own.

## Reproducing

On `hyperv-gha`, from `~/wspr-lab`:

    python3 tools/mkdop.py          # regenerate src/wsprd-dop from src/wsprd-mt
    ./build.sh dop                  # -> bin/wsprd-dop

Estimator against the simulator's known spread (about 4 minutes):

    python3 tools/dopcheck.py --trials 80 --jobs 10 --out results/dop/est_v3.json

Threshold sweep, five channels x five decoder settings, 500 trials per SNR
point on a 0.5 dB grid (about 90 minutes on 12 cores):

    bash tools/rundop.sh            # -> results/dop/<chan>_<variant>.json
    python3 tools/doptable.py       # the threshold table in section 2

Fixed-window control, six noncoherent trials plus one coherent window
(about 30 minutes):

    bash tools/rundopfix.sh         # -> results/dop/<chan>_w<N>.json
    python3 tools/dopfixtable.py    # the fixed-window table in section 3

Both sweeps skip any output file that already exists, so they can be
interrupted and restarted. `--jobs 10` leaves two of the twelve cores free;
`-P 1` in the decoder arguments is essential, or every parallel decoder spawns
a thread per core and the machine thrashes.

To watch the estimator on one file:

    WSPRD_DOPLOG=1 bin/wsprd-dop -C 500 -o 4 -d -N 14 -A 1 -P 1 -a . -f 10.1387 FILE.wav
