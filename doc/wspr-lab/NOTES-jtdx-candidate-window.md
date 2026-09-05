# JTDX's WSPR candidate-detection windows, ported and measured

Date: 2026-08-26. Lead from `results/forks/NOTES.md` section 2a.

**Verdict: it does nothing worth having. Do not ship it.**

## 0. Results

Every arm is the same binary, `bin/wsprd-jtw`, differing only in `-W`. Every arm
at a given SNR sees byte-identical wav files, so the arms are paired. 500 trials
per SNR point, six points on a 0.5 dB grid, 3000 files per arm per channel.
Arguments `-C 500 -o 4 -d -N 14 -P 1`, i.e. what would ship.

### 50% decode threshold, dB re 2500 Hz, callsign not in the hash table

| channel | baseline | JTDX (+w2 +w3) | change | decodes |
|---|---|---|---|---|
| AWGN, stable          | -32.35 | -32.35 | 0.00 dB | 1597 -> 1607 / 3000 |
| Rayleigh, 0.1 Hz      | -29.90 | -29.90 | 0.00 dB | 1631 -> 1640 / 3000 |
| Rayleigh, 1.0 Hz      | -26.30 | -26.30 | 0.00 dB | 1549 -> 1552 / 3000 |

The baselines reproduce `results/SUMMARY.txt` to 0.05 dB (-32.30, -29.85,
-26.30), which is the harness agreeing with itself on a different seed stream.

Twenty-two extra decodes in nine thousand files, which at the local slope of
about 0.45 per dB is under 0.01 dB. It is not a gain, it is noise with a sign.

### Which window, and is either a better detector at all

Appending a pass can only add decodes, so it cannot tell you whether a window is
good. `-W 4` and `-W 5` replace the half-sine in all three stock passes and
answer that directly. AWGN:

| arm | SNR50 | change | decodes /3000 |
|---|---|---|---|
| baseline, half-sine only        | -32.35 |  0.00 | 1597 |
| + one w2 pass                   | -32.35 |  0.00 | 1604 |
| + one w3 pass                   | -32.35 |  0.00 | 1600 |
| + both, i.e. JTDX               | -32.35 |  0.00 | 1607 |
| w2 REPLACING the half-sine      | -32.25 | +0.10 | 1502 |
| w3 REPLACING the half-sine      | -32.20 | +0.15 | 1455 |

Both JTDX windows are worse detection windows than the one already there. The
one-symbol "matched filter" costs 0.10 dB and 95 decodes; the two-symbol boxcar
costs 0.15 dB and 142 decodes. Bolting them on as extra passes hides that,
because the stock passes still run and still find everything they found before.

### The windows measured directly, pass-0 probe

All three windows applied to the same unmodified file, so no pass-order or
subtraction confound. "peak" is the strongest normalised bin within 1 Hz of the
truth before the min_snr floor, i.e. the raw detector output.

| channel | peak w1 - w2 | peak w1 - w3 |
|---|---|---|
| AWGN, -31.0 to -33.5    | w2 is 0.42 to 0.45 dB down | w3 is 0.12 to 0.13 dB down |
| Rayleigh 0.1, -28.5 to -31.0 | 0.43 to 0.44 dB down | 0.12 to 0.14 dB down |
| Rayleigh 1.0, -25.0 to -27.5 | 0.39 to 0.40 dB down | 0.15 dB down |

Flat across 13 dB of SNR and three channels, which is what a fixed processing
loss looks like. The candidate-formation rate follows: at -33.5 dB on AWGN, a
candidate lands within 1 Hz of the truth in 87.8% of files with the stock
window, 86.2% with w3 and 81.8% with w2.

### Candidate detection is not the bottleneck, and this is the useful finding

  `det` a candidate formed within 1 Hz of the truth in some pass;
  `syn` such a candidate also survived sync refinement and the minsync2 gate;
  `dec` the message decoded.  Baseline arm.

| channel | SNR | det | syn | dec |
|---|---|---|---|---|
| AWGN         | -31.0 | 0.996 | 0.956 | 0.948 |
| AWGN         | -32.5 | 0.972 | 0.520 | 0.400 |
| AWGN         | -33.5 | 0.878 | 0.258 | 0.108 |
| Rayleigh 0.1 | -29.0 | 1.000 | 1.000 | 0.810 |
| Rayleigh 0.1 | -31.0 | 0.996 | 0.888 | 0.164 |
| Rayleigh 1.0 | -26.0 | 1.000 | 0.994 | 0.648 |
| Rayleigh 1.0 | -27.5 | 1.000 | 0.910 | 0.084 |

**The premise this experiment was launched on is wrong.** With `-d` in effect,
candidate detection does not give out before the decoder does. On both fading
channels it is at or within a fraction of a percent of 1.000 while the decode
rate falls to 8%, so acquisition is not costing anything at all there. On AWGN
detection is still 88% where decoding is 11%.

Where the loss does sit, on AWGN at the very bottom, is between detection and
the demodulator: 87.8% detected, 25.8% surviving sync refinement and the
`minsync2` gate, 10.8% decoded. That is the sync stage, not the candidate list,
and nothing in the JTDX change touches it. It is a real lead and nobody has
looked at it.

### False decodes

Zero on AWGN in 3000 files either way. On the fading channels the two arms
produce the *same* false decodes, from the same files: one at Rayleigh 0.1 Hz
(-30.5 dB) and two at Rayleigh 1.0 Hz. Those are baseline false decodes, not
introduced by the change. The change introduces none. It also removes none.

### CPU per file

Measured on hyperv-gha (Intel i7-13700 VM, 6 physical / 12 logical cores) with
one decoder process at a time, `-P 1`, no `-Y`. The box was NOT quiet: load
average 15.2 before and 16.4 after, all of it other people's work. Both arms ran
back to back under the same conditions, so the ratio is meaningful even though
the absolute figures are inflated.

| corpus | baseline | JTDX | change |
|---|---|---|---|
| single signal, -28 dB, decodes on pass 0 | 0.185 s | 0.193 s | +4% |
| single signal, -32 dB, marginal          | 0.462 s | 0.692 s | +50% |
| single signal, -36 dB, nothing there     | 0.039 s | 0.040 s | +3% |
| busy sub-band, 25 signals per file       | 8.806 s | 10.930 s | +24% |

The busy sub-band is the case that matters, and it costs 24% more CPU for one
extra decode in 250. The marginal single-signal case costs 50% more, because
that is exactly where the two extra passes do their full work and find nothing.
The executed pass count goes from 2 to 4: pass 1 is skipped when pass 0 decodes
nothing, in both arms.

### Verdict

Do not ship it. It buys nothing on any channel, it costs a quarter of the
decoder's CPU on a busy band, and the window it is built around is measurably a
worse detector than the one already in the tree. The lead was worth chasing and
the answer is no.

## 1. What JTDX actually does

Source: `jtdx-project/jtdx`, `lib/wsprd/wsprd.c`, cloned to `/tmp/jtdx-clone`
(HEAD `2a0e2bea8`, 2022-03-01). The file's first line is
`//last time modified by Igor UA3DJY on 20191203`, which is the only attribution
that exists: the GitHub repository's history begins with a squashed "Initial
commit" on 2020-02-21, so the change predates any visible history and there is
no changelog entry, user-guide note or mailing-list post about it anywhere in
the tree.

`lib/wsprd/wsprd.c:949-956` defines two extra analysis windows next to the stock
half-sine:

    float w[512]; float w2[512]; float w3[512];
    for(i=0; i<512; i++) { w[i]=sin(0.006147931*i); }
    for(i=0; i<128; i++) { w2[i]=0.; }
    w2[128]=1.9; w2[383]=1.9;
    for(i=129; i<383; i++) { w2[i]=1.; }
    for(i=0; i<128; i++) { w2[i+384]=0.; }
    w3[0]=1.9; w3[511]=1.9;
    for(i=1; i<511; i++) { w3[i]=1.0; }

`wsprd.c:778-779` raises the pass count with the stock value commented out on
the line above:

    //    int npasses=3;
        int npasses=4;

and `wsprd.c:992` and `wsprd.c:1004` swap the window in, destructively, at the
top of the two later passes:

    if(ipass == 2 ) { for(i=0; i<512; i++) { w[i]=w2[i]; } ... }
    if(ipass == 3 ) { for(i=0; i<512; i++) { w[i]=w3[i]; } ... }

### Where the earlier survey was right, and where it was not

Right: `w2` is rectangular over samples 128..383, which at `nspersym=256` is
exactly one symbol, zero-padded to 512. That is the matched filter for a single
WSPR tone. `npasses` does go 3 to 4 with the old value commented out, and the
extra passes do redo the entire candidate search, coarse sync search and decode
attempt against the new spectrum.

Wrong, and it matters: **`w3` is not "the same idea across the full window"**.
It is a boxcar over all 512 samples, i.e. still *two* symbol times. It is not a
one-symbol matched filter and has nothing to do with the matched-filter
argument; it only replaces the half-sine taper with a rectangle. Whatever
`w3` is for, it is not what the lead said it was.

Also not in the survey, and worth knowing:

* Both windows set their two end samples to 1.9 rather than 1.0. There is no
  comment and no obvious reason. It is a hand tweak.
* JTDX's base decoder is a WSJT-X ~2.1 snapshot whose pass structure is not
  ours. Its passes 1, 2 and 3 all carry *identical* demodulator settings
  (`nblocksize=3`, `maxdrift=0`, `minsync2=0.10`); the only thing that
  distinguishes them after the change is the window. So in JTDX, pass 2 was
  already a straight repeat of pass 1 and UA3DJY repurposed it. Upstream's
  current three passes are not like that: passes 0 and 1 are single-symbol with
  drift search, pass 2 is the block/coherent pass. We have no redundant pass to
  repurpose.
* JTDX's `-d` is the older form: it steps every second bin and *replaces* the
  local-maximum search, where ours steps every third bin and *adds* to it.
* JTDX's `min_snr` is the stock -8 dB, with an older -7 commented out. This is
  not a `min_snr` change.
* The extra windows affect more than candidate detection. `ps[][]` is also what
  the coarse DT/frequency/drift search reads, so changing the window changes
  acquisition as well as detection. See section 3.

## 2. What was ported

`tools/mkjtw.py` generates `src/wsprd-jtw` from `src/wsprd-stage`, adding one
runtime switch so that both arms of every comparison are the same binary:

    -W 0   stock: 3 passes, half-sine only          (default)
    -W 1   JTDX:  5 passes, + one with w2, + one with w3
    -W 2   4 passes, the extra one with w2 only
    -W 3   4 passes, the extra one with w3 only
    -W 4   3 passes, w2 REPLACING the half-sine everywhere
    -W 5   3 passes, w3 REPLACING the half-sine everywhere

`-W 0` is byte-identical to the branch decoder: `tools/jtwident.py` runs
`bin/wsprd-stage` and `bin/wsprd-jtw -W 0` over the same files with the same
arguments and compares stdout exactly. 216 files across six SNRs and three
channels, zero mismatches.

Our pass structure is passes 0 and 1 single-symbol with `maxdrift=4`, pass 2
block/coherent with `nblocksize=nbtrials`, `maxdrift=0`, `minsync2=0.10`. The
faithful port of "add extra passes that differ only in the window" is therefore
to append passes that repeat pass 2's settings, which is what `-W 1` does. That
is also the generous reading: every extra pass is an extra decode attempt on top
of an unmodified stock search, so it can only ever add decodes. Repurposing our
pass 2 the way JTDX repurposed its own would have thrown away the coherent
demodulator, which would not have been a fair test.

Rather than clobbering `w[]` in place, the port selects a window per pass through
a pointer. That is equivalent and non-destructive.

`-W 4` and `-W 5` are not what JTDX does. They exist because the "append a pass"
experiment cannot tell you whether the window itself is better or worse, only
whether an extra attempt finds anything. Replacing the window in all three stock
passes answers the actual signal-processing question.

### Instrumentation

`-Y <hz>` gives the decoder the true signal offset from 1500 Hz and makes it
report on stderr, read-only, for every pass:

    #CANDDET pass N win N npk N mindf X snr X    candidate list, per pass
    #CANDSYN pass N nwat N mindf X sync X        what survived the minsync2 gate
    #PROBE   win N npk N mindf X csnr X peaksnr X above N

`#PROBE` runs on pass 0 only and builds the detection spectrum three times, once
per window, on the *same unmodified file*. That is the comparison with no pass
order or subtraction confound: it answers "given this file, which window puts a
candidate nearer the signal, and how strong is the bin it puts there".

## 3. Why it cannot work: the mechanism

`tools/jtwtheory.py` reimplements the detector in numpy and needs no decoder and
no wav files. Run it to see all of this.

**Per-frame processing gain for a stationary tone**, `(sum w)^2 / sum(w^2)`:

    w1 half-sine, 2 sym    +26.17 dB
    w2 boxcar,    1 sym    +24.06 dB
    w3 boxcar,    2 sym    +27.08 dB

The one-symbol window integrates for half as long, so it starts about 2 dB
behind the stock window and 3 dB behind the two-symbol boxcar. The matched-filter
argument has to buy back at least that much before it breaks even.

**It cannot, because of the 7-point smoothing.** The candidate spectrum is not
one FFT bin, it is `smspec[i] = sum of 7 adjacent psavg bins`, which spans
7 x 0.7324 = 5.13 Hz. The four WSPR tones span 3 x 1.4648 = 4.39 Hz. Every tone
is inside the smoothing window. A two-symbol window does smear each tone across
bins, exactly as the lead said, but it smears it into bins that are summed
anyway, so nothing is lost. Removing the smoothing changes the ordering; leaving
it in, the one-symbol window is simply 0.4 to 0.5 dB down.

**The one place the argument is sound is acquisition, and there is no headroom
there.** `ps[][]` is read directly, unsmoothed, by the coarse DT/frequency/drift
search, which correlates the sync vector one symbol at a time. A two-symbol
window puts two symbols' sync bits into every value that search reads, and
adjacent sync bits agree only half the time, so half the correlation really is
thrown away. Noiseless coarse sync at the true offset:

    w1 half-sine, 2 sym    0.607
    w2 boxcar,    1 sym    0.973      <- the matched filter, and it shows
    w3 boxcar,    2 sym    0.492

But the coarse search does not need to be good, it only needs to be right, and
the stock window already is. Fraction of trials where the coarse search's argmax
lands on the exactly correct time offset and frequency bin, 60 trials per point:

    gen SNR      w1        w2        w3
      -26.0    1.000     1.000     0.300
      -29.0    1.000     1.000     0.383
      -32.0    1.000     1.000     0.367

The half-sine acquires perfectly at every SNR at which a decode is even
possible, so a better acquisition metric buys nothing. `w3`, meanwhile, is
actively bad at this: its boxcar leakage lands the coarse frequency estimate on
the wrong bin two thirds of the time, even with no noise at all.

## 4. Reproduction

    cd /home/tf/wspr-lab
    python3 tools/mkjtw.py          # regenerates src/wsprd-jtw from src/wsprd-stage
    ./build.sh jtw                  # -> bin/wsprd-jtw

    python3 tools/jtwident.py 5     # proves -W 0 is byte-identical to the branch
    python3 tools/jtwtheory.py      # the mechanism, in numpy, no decoder involved

    ./tools/jtwcampaign.sh  A 5     # whole-dB half of the SNR grid, this box
    ./tools/jtwcampaign2.sh A 5     # the -W 4 / -W 5 follow-up
    ssh hyperv-gha './tools/jtwcampaign.sh B 5; ./tools/jtwcampaign2.sh B 5'
    python3 tools/jtwtable.py       # merges both halves, prints every table

    ./tools/jtwcpu.sh               # CPU per file, one decoder process at a time

The SNR grid is split between the two boxes so the halves interleave on 0.5 dB:
this box takes the whole dB, hyperv-gha the half dB. `tools/jtwrun.py` seeds each
trial from the SNR *value* rather than from its position in the list, so every
arm at a given SNR sees byte-identical wav files and the arms are paired. The
channel simulator is bit-reproducible across the two machines; this was checked.

Note that `-march=native` differs between the boxes, so the source must be synced
and rebuilt on the remote, not the binary copied.

Files: `results/jtw/{chan}_{arm}_{A,B}.json`, in exactly `tools/harness.py`'s
format so `tools/fit.py` and `tools/mergejson.py` work on them unchanged, plus
the candidate-detection counters. Scoring is imported from `tools/harness.py`
rather than reimplemented, so these thresholds are directly comparable with
everything in `results/SUMMARY.txt`.
