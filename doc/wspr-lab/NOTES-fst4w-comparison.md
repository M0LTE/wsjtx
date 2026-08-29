# FST4W-120 against WSPR-2: how these numbers were made

Everything here is reproducible from this checkout plus a WSJT-X build.
Raw per-SNR data: `curves.tsv`, `thresholds.json`, and one `*.json` per run.
Headline table: `SUMMARY.txt`.

## Answer

50% decode threshold, SNR in dB relative to 2500 Hz, callsign NOT previously
heard, one signal per file, frequency uniform in +/-95 Hz of 1500 and DT uniform
in -1.0..+1.5 s. 600 trials per SNR point on a 0.5 dB grid.

| channel | FST4W-120 | WSPR stock | WSPR patched | FST4W vs stock | FST4W vs patched |
|---|---|---|---|---|---|
| AWGN, stable path | **-32.22** | -30.94 | -32.31 | **-1.29** +/-0.04 | **+0.09** +/-0.04 |
| Rayleigh, 0.1 Hz spread | **-30.40** | -29.18 | -29.95 | **-1.22** +/-0.05 | **-0.44** +/-0.07 |
| Rayleigh, 1.0 Hz spread | **-26.16** | -26.25 | -26.27 | **+0.10** +/-0.05 | **+0.11** +/-0.06 |

Negative means FST4W-120 needs less signal, i.e. is better. Numbers are the
model-free 50% crossing; the maximum-likelihood logistic fit agrees to within
0.06 dB on the first two rows. On the 1 Hz row the logistic does NOT fit FST4W
(chi-square 257 on 14 dof) because that curve has a long shoulder rather than a
clean S, so the model-free number is the one to use there.

The WSPR columns were re-measured here from scratch and reproduce
`results/SUMMARY.txt` to within 0.07 dB on all six values (-30.94 vs -30.90,
-32.31 vs -32.30, -29.18 vs -29.20, -29.95 vs -29.85, -26.25 vs -26.30,
-26.27 vs -26.30), which is the end-to-end check that the two harnesses are
measuring the same thing.

### Side conditions, AWGN

| condition | SNR50 | vs the -d 3 baseline |
|---|---|---|
| `-d 3`, new callsign (the number above) | -32.22 | - |
| `-d 3`, call already in `fst4w_calls.txt` | -32.77 | 0.55 dB better |
| `-d 2` | -32.20 | 0.02 dB |
| `-d 1` | -32.08 | 0.14 dB worse |
| `-d 3` with F Tol 20 Hz instead of 100 Hz | -32.24 | 0.02 dB |

Two things fall out of that. FST4W has the same shape of already-known-station
bonus that WSPR does (0.55 dB, against WSPR's 0.90 dB), for the same reason: at
`-d 3` a decode that fails the CRC path is retried with the CRC thrown away
(Keff=50) and accepted only if the call/grid is already in `fst4w_calls.txt`. So
comparing FST4W against WSPR on the new-callsign case is the right pairing, and
both modes hold something back for stations you have heard before. And the
200 Hz search window costs FST4W essentially nothing - narrowing it to 40 Hz
buys 0.02 dB - so the wide-window choice is not tilting the comparison.

### Above threshold, under fast fading

The 50% threshold is not the whole story at 1 Hz Doppler spread. FST4W-120 has
a shoulder there that WSPR does not:

| SNR | FST4W-120 | WSPR stock | WSPR patched |
|---|---|---|---|
| -24.0 | 543/600 (90.5%) | 600/600 | 600/600 |
| -25.0 | 502/600 (83.7%) | 549/600 | 552/600 |

At -24 dB, 6 dB above threshold, FST4W is still losing one transmission in ten
while WSPR loses none. Pushing further up: 590/600 at -21, 596/600 at -20 and
-19, 598/600 at -18. FST4W's sync is a coherent correlation over 8-symbol
(5.5 s) blocks, and 1 Hz of Doppler spread decorrelates the channel in about
0.27 s, so the sync statistic is the plausible culprit. WSPR's threshold is
0.10 dB better than FST4W's at 1 Hz spread and its high-SNR behaviour is much
better.

### CPU per file

Child user+sys CPU, one process at a time on an otherwise idle box, 40 files per
cell (`cpu_per_file.json`):

| channel | SNR | wsprd stock | wsprd patched | jt9 -W FST4W-120 |
|---|---|---|---|---|
| AWGN | -24 | 0.113 s | 0.813 s | 0.169 s |
| AWGN | -32 | 0.144 s | 1.111 s | 0.439 s |
| Rayleigh 0.1 Hz | -32 | 0.122 s | 1.322 s | 0.537 s |
| Rayleigh 1.0 Hz | -32 | 0.050 s | 0.559 s | 0.412 s |

FST4W-120 reaches the patched WSPR decoder's AWGN sensitivity for about 40% of
its CPU, and the stock decoder's for about 3x the CPU. Note that the patched
`wsprd` defaults to one worker thread per logical processor - these are its
single-thread costs (`-P 1`).

### False decodes

| decoder | false decodes | files |
|---|---|---|
| FST4W-120 `jt9 -W` | 2 | 49950 |
| WSPR stock `wsprd` | 0 | 24000 |
| WSPR patched `wsprd -N 14` | 2 | 24000 |

Nothing separates the three at this sample size. FST4W's two both came through
the 24-bit CRC path (`8W9FWJ BF93 27` at -33.5 dB, `8E6KUS RA79 10` at -33.0 dB
in the already-known run); the patched WSPR's two were compound-callsign
nonsense at -26 and -27 dB under 1 Hz fading.

## What is being compared

| | WSPR-2 | FST4W-120 |
|---|---|---|
| period | 120 s | 120 s |
| transmission | 110.6 s | 109.3 s |
| symbols | 162 | 160 |
| baud | 1.4648 | 1.4634 |
| occupied bw | 5.9 Hz | 5.9 Hz |
| payload | 50 bits | 50 bits |
| error control | K=32 r=1/2 convolutional, Fano/OSD, **no CRC** | (240,74) LDPC carrying 50 payload bits plus a **24-bit CRC** |
| decoder | `wsprd` | `jt9 -W` |

Same period, near-identical transmission length, same occupied bandwidth, same
payload. That is what makes FST4W-120 the fair comparison for WSPR.

## Binaries

    /tmp/verify-tree                 WSJT-X source with the WSPR patch applied
    /tmp/verify-tree/build/jt9       FST4W decoder
    /tmp/verify-tree/build/wsprd     patched wsprd      (bin/wsprd-verify -> this)
    bin/wsprd-baseline               pristine wsprd
    bin/fst4simx                     signal generator, see below
    bin/wsprchan                     WSPR channel simulator

`bin/fst4simx` is built from `tools/fst4simx.f90`, which is
`lib/fst4/fst4sim.f90` with exactly two additions: a 10th argument that seeds
libc `rand()`, and an optional 11th argument giving the output filename. Stock
`fst4sim` never calls `sgran()`, so every invocation would otherwise produce the
identical noise realisation, which is useless for statistics. Seeded with 1 it
is byte-identical to stock `fst4sim`, checked for both AWGN and fading:

    /tmp/verify-tree/build/fst4sim "M0LTE IO91 30" 120 1500 0.0 0.0 0.0 1 -28 T
    bin/fst4simx                    "M0LTE IO91 30" 120 1500 0.0 0.0 0.0 1 -28 T 1 a1.wav
    md5sum 000000_0001.wav a1.wav      # identical

Build it with:

    cd /tmp/verify-tree/build
    gfortran -c -O3 -funroll-all-loops -fno-f2c -fno-second-underscore \
      -I/tmp/verify-tree/build -I/tmp/verify-tree -I/tmp/verify-tree/lib/fst4 \
      /home/tf/wspr-lab/tools/fst4simx.f90 -o /tmp/fst4simx.o
    c++ -O3 -fopenmp -pthread /tmp/fst4simx.o -o /home/tf/wspr-lab/bin/fst4simx \
      libwsjt_fort.a libwsjt_cxx.a -lfftw3f_threads -lfftw3f -lm \
      -lboost_log_setup -lboost_log -lboost_filesystem -lboost_regex \
      -lboost_serialization -lboost_thread -lboost_atomic -lboost_chrono \
      -lboost_container -lboost_date_time -lgfortran -lquadmath

(the exact link line is `build/CMakeFiles/fst4sim.dir/link.txt`).

## Driving jt9 for FST4W-120

    jt9 -W -p 120 -f 1500 -F 100 -d 3 -a <dir> -t <dir> <NNNN_NNNN.wav>

  * `-W` selects FST4W (mode 241). `-Y` is the same plus hash22 printing.
  * `-p 120` is the T/R period, which is what picks nsps=8200 / NN=160.
  * `-f 1500` is the Rx frequency, `-F 100` the tolerance; jt9 caps the
    tolerance at 100 Hz for FST4W, giving a 200 Hz search window - the same
    width as the WSPR sub-band `wsprd` searches. Verified by decoding signals
    placed at 1405, 1450, 1500, 1550, 1595 and 1601 Hz: all six decode. The
    source puts the window edge at nfqso + 1.5*baud + ntol = 1602.2 Hz.
  * `-d 3` is Deepest. `mainwindow.cpp` drives BOTH decoders from the one
    Decode-menu setting: Deepest sends `wsprd -C 500 -o 4 -d [-N 14]` and
    `jt9 -d 3`. So `-d 3` against `wsprd -d` is the like-for-like pair.
  * `-a` and `-t` must point at a scratch directory: jt9 writes
    `fst4w_calls.txt`, `decoded.txt`, `timer.out` and `jt9_wisdom.dat` there.
  * the filename must be long enough for jt9's UTC parser - it reads 4 digits
    before `.wav` if the character 5 back is `_`, else 6 digits, and indexes
    off the end of a short name without checking. `d.wav` crashes it;
    `0001_0000.wav` is fine.
  * DT search covers -1.0 s to +2.0 s of nominal. Verified: generated DT of
    -1.0/-0.5/0.0/0.7/1.5/2.0 all decode and are reported back exactly; 2.5
    does not. (`emedelay` is never initialised for FST4W in `jt9.f90`; it lands
    on zeroed heap in practice, which is the branch that gives 0..3 s.)

## FST4W's equivalent of the WSPR hash table

`fst4_decode.f90` reads `<data_dir>/fst4w_calls.txt` at startup. At `-d 3` only,
a decode that fails the 24-bit CRC path is retried with Keff=50 - no CRC at all
- and is accepted only if the call/grid it produces is already in that file.
That is precisely WSPR's hash-table gate in a different costume, so the harness
treats it the same way: a fresh scratch directory and a fresh random callsign
per trial, so "this station has not been heard before" is the default case.
`--wcallseed` pre-loads the file with the transmitted call/grid to measure the
already-known case, mirroring the WSPR `--hashseed` oracle test.

Nothing else in the FST4W path is a-priori: for 50-bit messages the code sets
`ntmax=nblock` and skips all the AP passes (`iwspr.eq.1` branch), so `-c`/`-x`
callsign hints cannot help it.

## SNR calibration - the part that mattered most

Both simulators state SNR in dB relative to a 2500 Hz noise bandwidth, and both
scale the signal amplitude by exactly `10**(snr/20)`, so only the offset needed
checking. Three independent checks:

1. **From the source.** `wsprchan.c`: noise rms 1000 over 0..6000 Hz, signal
   amplitude `1000*sqrt(lin*4*2500/12000)`, so signal power / noise-in-2500 Hz
   `= lin`. `fst4sim.f90`: `bandwidth_ratio=2500/6000`,
   `sig=sqrt(2*bandwidth_ratio)*10**(0.05*snrdb)`, gain 100 on unit-variance
   noise, giving the same `lin`. Identical conventions.

2. **From the samples** (`tools/snrcal.py`, output in
   `snr_calibration_absolute.txt`). Both simulators are deterministic given the
   seed and generate their noise independently of the requested SNR, so
   `file(snr) - file(-100 dB)` isolates the signal exactly. Measuring signal and
   noise power that way, and correcting for the two int16 quantisers the
   difference carries, the two agree to **0.014 dB or better on AWGN** across
   0/-10/-20/-28 dB, and to 0.06 dB on the fading channels where the WSPR-side
   sampling error is itself +/-0.10 dB.

3. **From the decoder.** Generating FST4W at -15/-20/-25/-28 makes jt9 report
   -15.00/-20.00/-24.97/-27.98, the same test that made `wsprd` report
   -15/-20/-25/-28. Below about -30 the reported value flattens out (-32.92 at a
   generated -34.0) because only the favourable noise realisations decode.
   Under 1 Hz Doppler spread jt9's FST4W SNR estimate runs about 1.4 dB low
   (generated -24.0 reported -25.43); that is an estimator bias in the report
   line only, and does not touch the generated SNR, which check 2 pins directly.

## Fading conventions

`wsprchan -R W` and `fst4sim`'s `fdop` argument are the same quantity: the
2-sigma width of a Gaussian Doppler power spectrum. From the source, wsprchan
shapes a complex Gaussian by `exp(-f^2/(2*sigma^2))` in power with
`sigma = W/2`; `watterson.f90` shapes by `exp(-(f/fspread)^2)` in amplitude,
i.e. `exp(-2f^2/fspread^2)` in power, also `sigma = fspread/2`.

Measured as well (`tools/fadecheck.py`, output in
`fading_convention_check.txt`): the half-width of the autocorrelation of the
fading power envelope matches the analytic `0.2650/W` for both simulators at
W = 0.1, 0.3, 1.0 and 3.0 Hz, and matches between them.

Not used here: `fst4sim` also accepts a NEGATIVE fdop, which switches to
`lorentzian_fading.f90`, a Lorentzian rather than Gaussian Doppler spectrum with
its own width convention (`b=6`, so the parameter is roughly the -20 dB width).
Everything below uses the positive/Gaussian path, which is the one that matches
the WSPR lab.

## The harness

`tools/harness_fst4w.py` is the FST4W twin of `tools/harness.py` and keeps the
same discipline:

  * one scratch directory per trial on /dev/shm, deleted afterwards
  * a fresh random callsign per trial
  * frequency offset uniform in +/-95 Hz of 1500, DT uniform in -1.0..+1.5 s,
    the same ranges the WSPR harness uses
  * the simulator's own `Message:` line is checked against what was requested,
    so a message the packer silently mangled cannot be scored as a failure
  * scoring is exact on (call, grid, power) after normalising the power to an
    integer; anything else the decoder prints is a false decode

One difference from the WSPR harness, and it does not matter: `packjt77`
refuses `Q` as a callsign prefix (ITU reserves Q for Q-codes), so
`harness_fst4w.py` drops Q from the prefix alphabet. The old WSPR packer accepts
it. The callsign is 28 bits of payload in both modes and the decoder never sees
the letters, so this cannot shift a threshold.

## Reproducing

    cd /home/tf/wspr-lab
    python3 tools/snrcal.py 16                    # SNR calibration
    python3 tools/fadecheck.py 10                 # fading convention
    bash    tools/runall2b.sh                     # every curve (about 90 min on 10 cores)
    python3 tools/fst4w_summary.py                # -> SUMMARY.txt, curves.tsv, thresholds.json
    python3 tools/cpubench.py 40                  # CPU per file, needs an idle box

`tools/fit_ci.py` fits one curve: maximum-likelihood logistic with a
profile-likelihood 95% confidence interval on the 50% threshold, a Pearson
chi-square goodness of fit, and a model-free 50% crossing by linear
interpolation between the bracketing grid points as a cross-check. It will fit
a ceiling below 1 if a likelihood-ratio test demands one, which is why every
curve here is measured up to a decode rate of at least 0.99 - otherwise a grid
that simply stops before the curve flattens gets mistaken for a real outage
floor.

## Caveats

  * `wsprd` in the patched tree defaults to one worker thread per logical
    processor. Inside a parallel harness it must be given `-P 1`, or the box is
    oversubscribed 16x. Decoder output is identical at any `-P`
    (`tools/determinism.py`), only the timing changes.
  * The wall-clock `s/file` column in the run logs was taken while another
    session was using about 6 of this box's 16 cores, so treat it as throughput
    under contention, not as CPU cost. `cpu_per_file.json` is the clean number,
    measured as child user+sys CPU one process at a time.

## What this does NOT establish

  * **Crowded sub-bands.** Every trial here has one signal in the file. WSPR's
    decoder does multi-pass subtraction and the lab has separate busy-band
    numbers for it; `jt9 -W` does no subtraction pass at all, so its behaviour
    with 25 signals in 200 Hz is unmeasured and is very likely worse than these
    single-signal numbers suggest.
  * **Impulsive noise.** The blanker sweep that is worth 7.4 dB to the patched
    WSPR decoder under heavy QRN has no counterpart here. FST4W has a blanker
    (`nexp_decode/256 - 3`, default fixed 0%), and driving it from jt9 needs
    `-X`, which was not exercised.
  * **Drift.** Not tested for either mode.
  * **Other Doppler models.** `fst4sim`'s negative-fdop Lorentzian path is not
    used; only the Gaussian/Watterson path that matches `wsprchan -R`.
  * **A real band.** Everything is simulated. Nothing here has been checked
    against off-air recordings.
  * **jt9's reported SNR under fading.** It runs about 1.4 dB low at 1 Hz
    Doppler spread (generated -24.0 reported -25.4). That is a defect in the
    report line, not in the measurement, because the generated SNR was pinned
    directly from the samples - but it means jt9's own SNR column cannot be used
    to cross-calibrate a fading channel.
