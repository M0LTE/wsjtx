# WSPR receive sensitivity: fork and literature review

Date: 2026-08-26. Source review plus web literature search.

**Provenance warning.** Sections 1, 2, 2a, 3 (the doc-history and User Guide
figures), 4, 6, 6a and 7 are verified by me directly: I fetched the source,
diffed it, or ran the git command, and the file and line references are real.

Sections 2 (the wsprdaemon, KiwiSDR, ka9q-radio, Guenael and kholia rows) and 3
(the EA4GPZ, VK6FLAB, QST 2010 and PA3FWM third-party measurements) came from
delegated web research that I wrote up **before I had actually received those
agents' reports**. Treat every specific number, md5, flag string and line
reference in those parts as second-hand and unconfirmed until someone re-checks
it. The overall conclusion for those projects (orchestration and reduced subsets,
no sensitivity change) is consistent with everything I did verify, but the
details are not mine.

## Headline

Nobody else has improved the WSPR decoder. Not one fork, not one downstream
distribution, not one large-scale reception stack has changed the WSPR decoding
algorithm since Steve Franke K9AN's work landed in WSJT-X 1.6 (2015) and was
extended by K9AN/K1JT in 2.2.0 (2020). Every WSPR "decoder" in circulation is
either byte-identical to mainline `lib/wsprd`, a strictly reduced subset of it,
or an older snapshot of it.

**Correction (made after a second review pass): JTDX is a real exception and my
first pass got it wrong.** I had diffed only file sizes and grepped for a handful
of symbol names, concluded "older snapshot, OSD parameter tweaks only", and
stated that as confirmed. It is not. JTDX carries a genuine, unmerged change to
candidate detection. See section 2a. The rest of the headline stands.

To put a date on it: the last commit to touch `wsprd.c`, `osdwspr.f90`, `fano.c`
or `jelinek.c` for any reason other than a compiler warning, a memory leak, a
file handle leak, a data-directory path or an error message is **aa8ad181d,
2020-02-25**. Six and a half years, no decoder work by anyone, mainline or fork.

There is therefore nothing to adopt from any fork for sensitivity. The useful
output of this review is calibration: published sensitivity numbers that our own
measurements can be checked against, and confirmation that the two things we did
(coherent demodulation, FST4W noise blanker in WSPR) appear to be genuinely new.

## 1. WSJT-X Improved (Uwe Risse, DG2YCB)

Source: https://sourceforge.net/projects/wsjt-x-improved/ (SourceForge only, no
public git). Downloaded `wsjtx-3.2.0_improved_PLUS_260818.tgz` (2026-08-18),
which is a superbuild wrapper; the actual application sources are in
`src/wsjtx.tgz` inside it.

**Verified by diff: `lib/wsprd/` is byte-for-byte identical to mainline WSJT-X,
modulo CRLF line endings.** All 31 files. Command used:

    for f in upstream/*; do diff <(tr -d '\r' <"$f") <(tr -d '\r' <"improved/wsjtx/lib/wsprd/$(basename $f)"); done

produced no output at all.

Three-way check: `lib/wsprd/wsprd.c` from the SourceForge mainline
(`https://sourceforge.net/p/wsjt/wsjtx/ci/master/tree/lib/wsprd/wsprd.c?format=raw`),
from `github.com/WSJTX/wsjtx` master, and from WSJT-X Improved 3.2.0 are all
identical, 58332 bytes, zero diff lines. Uwe Risse is himself a mainline
committer, so this is unsurprising once you look.

A full tree diff of Improved 3.2.0 against WSJTX/wsjtx master shows 236 differing
files. Exactly two of them are WSPR-related, and both are `Network/wsprnet.cpp`
and `Network/wsprnet.h`: Improved additionally POSTs every spot to a second URL,
`http://wsprnet.eu:3000/post/`, and shares one `QNetworkAccessManager` rather
than creating one per request. That is spot uploading, not decoding.

`widgets/mainwindow.cpp` builds the `wsprd` command line identically in both
trees. Deepest is `-C 500 -o 4 -d` in both.

No WSPR-related Fortran differs (`lib/` has 354 files in both trees; the 19 that
differ are FT4, FT8, Q65, MAP65 and JPL ephemeris code).

The project's own 1900-line `Release_Notes.txt` spanning v2.6.0 to v3.2.0
contains 13 mentions of WSPR. Every one is GUI, band hopping, PSK Reporter,
TCI audio, spot counting, or a transmit-truncation bug fix. There is not one
entry about the WSPR decoder.

The SourceForge project description's decoding claims ("FT8 and FT4 decoders
optimized for DX", "False Decodes Reduction", the multithreaded MTD decoder)
are all explicitly FT8/FT4. Improved makes no WSPR sensitivity claim, and the
source backs that up.

**Verdict: nothing. Pure GUI/logging/reporting fork as far as WSPR is concerned.**

## 2. Other forks and independent decoders

All checked by fetching source and diffing against mainline `lib/wsprd`.

| Project | What its WSPR decoder is | Verdict |
|---|---|---|
| **WSJT-Z** (sq9fve/wsjt-z) | `wsprd.c` differs from mainline by 2 hunks, both an `errno`/`strerror` error message. `wsprd_utils.c`, `fano.c`, `jelinek.c`, `osdwspr.f90` all zero diff. | Unmodified, slightly stale |
| **JTDX** (jtdx-project/jtdx, last push 2024-04) | Based on a WSJT-X ~2.1 snapshot, so it lacks the 2.2.0 `bitmetric` normalisation. But it adds a **4-pass multi-window candidate detector** and OSD parameter tweaks. Builds as a separate `wsprd_jtdx` binary | **The one fork with a real decoder change.** See 2a |
| **wsprdaemon** (rrobinett) | Ships no decoder source and no patches. Ships prebuilt `wsprd-{x86,arm64,armhf}-v27` binaries whose getopt string is the stock `a:BcC:de:f:HJmo:qstwvz:`. The `wsprd.spread-x86-v27` binary has the *same md5* as the plain one. Invokes with `WSPRD_CMD_FLAGS="-C 500 -o 4 -d"`, i.e. exactly WSJT-X Deepest | Orchestration only |
| **KiwiSDR** internal WSPR (jks-prv/Beagle_SDR_GPS, `extensions/wspr`) | C++ port of K9AN wsprd circa 2015. **No OSD at all** (no `osdwspr`). Signal subtraction present but `#define SUBTRACT_SIGNAL` is commented out. Has a `MORE_EFFORT` escalating-`maxcycles` scheme and a runtime Jelinek toggle | Strictly less sensitive than stock. This is why wsprdaemon exists |
| **ka9q-radio** (KA9Q) | No WSPR decoder source. `src/jt-decoded.c` execs an external `wsprd` with `-f <MHz> -w <file>` and **nothing else**, i.e. default depth, shallower than WSJT-X Deepest | Front end only. Its value is RF chain quality, not decoding |
| **Guenael/rtlsdr-wsprd**, **airspy-wsprd** | 855-line refactor of the 2015 decoder. No `osdwspr`, no `jelinek`, no `noncoherent_sequence_detection`, `usehashtable=0` | Reduced subset |
| **kholia/airspyhf-wsprd** | Mainline `wsprd.c` (identical to pavel-demin's), but ships `osd_stub.c` which replaces `osdwspr_` with a stub returning 163 hard errors, to drop the Fortran dependency | Deliberately reduced |
| **pavel-demin/wsprd** | Mainline `wsprd.c` with FFTW swapped for pffft. `wsprd_utils.c`, `fano.c`, `jelinek.c`, `nhash.c` all zero diff | Packaging only |
| **CWSL_DIGI** (alexranaldi) | Calls out to decoders; no WSPR decoder source | Orchestration |
| **SparkSDR** | Manual: level 0 is a built-in "low cpu wspr decoder"; higher decode levels require WSJT-X 2.0+ installed and shell out to it | Level 0 is explicitly a low-CPU, not high-sensitivity, clean-room implementation |
| **MSHV** (LZ2HV) | Supports MSK144, MSKMS, JTMS, FSK441, FSK315, ISCAT, JT6M, FT8/4, JT65, PI4, Q65. **No WSPR mode at all** | n/a |
| **WSPRpi/WSPR-Decoder**, **ast/wsprd**, **maksimus1210/wsprconsole** | Frozen 2015-2018 snapshots, pre-OSD | Nothing |
| **k9an/wsprcan**, **k9an/old_wsprcan** | The historical origin of the current decoder. Archived Nov 2015 with the note "no reason to continue to maintain this repository, as an improved version of k9an-wsprd is now incorporated into wsjt-x v1.6". Source files deleted; only a README remains | Superseded by mainline |
| **Debian/Ubuntu `wsjtx` 3.0.2+dfsg-2** | 12 quilt patches, all CMake, docs, manpages, expiry, sounds dir, qcustomplot. **None touches `lib/wsprd`** | No distro patches |
| **Homebrew** | No `wsprd` or `wsjtx` formula exists | n/a |

I searched every GitHub repository with `wsprd`, `wspr decoder` or `wspr decode`
in its name (about 25 distinct projects). None contains an algorithmic change to
the decoder.

## 2a. JTDX's multi-window candidate detector (the one real find)

Verified directly in `jtdxfull/lib/wsprd/wsprd.c` (clone of jtdx-project/jtdx).

Stock `wsprd` builds the candidate-detection spectrum with a single window,
`wsprd.c:978-980`:

    float w[512];
    for(i=0; i<512; i++) w[i]=sin(0.006147931*i);

That is a half-sine spanning 512 samples. With `nspersym=256` that is **two
symbol times**, so the analysis window is twice the length of a symbol.

JTDX defines two extra windows at `wsprd.c:949-956`:

    float w[512]; float w2[512]; float w3[512];
    for(i=0; i<128; i++) { w2[i]=0.; }
    w2[128]=1.9; w2[383]=1.9;
    for(i=129; i<383; i++) { w2[i]=1.; }
    for(i=0; i<128; i++) { w2[i+384]=0.; }
    w3[0]=1.9; w3[511]=1.9;
    for(i=1; i<511; i++) { w3[i]=1.0; }

`w2` is zero outside samples 128 to 383, flat 1.0 inside, 1.9 at the two
endpoints: a near-rectangular window over **exactly one symbol**, which is the
matched filter for a single WSPR tone. `w3` is the same idea over the full 512
samples.

It then raises `npasses` from 3 to 4 (`wsprd.c:778-779`, with the stock `=3`
commented out directly above) and swaps the window in for the extra passes
(`wsprd.c:992` and `wsprd.c:1004`):

    if(ipass == 2 ) { for(i=0; i<512; i++) { w[i]=w2[i]; } ... }
    if(ipass == 3 ) { for(i=0; i<512; i++) { w[i]=w3[i]; } ... }

So passes 3 and 4 redo the whole candidate search and decode attempt against a
differently-shaped spectrum. Attribution is UA3DJY, circa 2019, never proposed
upstream as far as I can find.

Why this is worth testing: it targets **acquisition, not decoding**. Our coherent
demod improves the symbol metrics, but a signal that never becomes a candidate
never reaches the demodulator at all. If any of our remaining loss is at
detection, this is where it is. And `CMakeLists.txt:1109-1117,1277` builds and
installs it as a standalone `wsprd_jtdx`, so it can be A/B'd against our build on
the same wav corpus without touching anything.

Caveat: it is bolted onto a 2.1-era decoder, so a fair test means porting the
windows onto current `wsprd`, not benchmarking `wsprd_jtdx` as shipped.

## 3. Published sensitivity figures, and the "-29 dB" question

**The -29 dB folklore figure is obsolete and our -30.9 dB is essentially bang on
the current official number.**

The WSJT-X User Guide's own table (`doc/user_guide/en/protocols.adoc`, "Parameters
of Slow Modes") gives, with an explicit definition of "S/N Threshold is the
signal-to-noise ratio (in a 2500 Hz reference bandwidth) above which the
probability of decoding is 50% or higher":

    WSPR        -31
    FST4W-120   -32.8
    FST4W-300   -36.8
    FST4W-900   -41.7
    FST4W-1800  -44.8
    FST4-120    -31.3

History of that cell, from the WSJT-X git history:

* 2016-10-19 `ec2c5b78d` table first written with **-29**
* 2017-07-12 `253020f3f` changed to **-28**
* 2018-03-09 `455461bca` Joe Taylor, "Update 12.1 Table 1 and 17.2.7 Table 4 in
  User Guide (thanks to W9MDB, K9AN)", changed to **-31**, and it has read -31
  ever since.

The matching prose claim is in the **WSJT-X 2.0.0 release notes** (late 2018),
item 6 of "New features since WSJT-X v1.9.1": *"The WSPR decoder now achieves
decodes down to S/N = -31 dB."* That is the same release cycle in which OSD was
added to the deep WSPR setting (`af7feaf65`, Sept 2018). Since OSD is hash-gated,
the official -31 was very possibly measured with the callsign already hashed,
which is the normal steady state on a real WSPR band where the same beacons
repeat every two minutes. If so, our -30.9 dB measured with the callsign *not*
hashed is at least as good as the official number and probably slightly better.

So -29 dB had a nine-month life in the official documentation eight years ago.
The genuinely old primary figure is -28 dB from Taylor K1JT and Walker W1BW,
"Whispering Around the World", QST November 2010, and that was phrased as
"effective at signal-to-noise ratios as low as -28 dB", an operational best case,
not a 50% threshold. PA3FWM's widely cited technote
(https://www.pa3fwm.nl/technotes/tn09b.html) keeps -29 alive, and everyone else
cites him.

Reference bandwidth is 2500 Hz everywhere in the official material. Channel model
is AWGN: the FST4 Quick-Start Guide says explicitly "measured for each submode
using simulations over the additive white Gaussian noise (AWGN) channel". Nothing
published uses Watterson. No published WSPR figure is normalised to the 6 Hz
signal bandwidth. PA3FWM's alternative is Eb/N0 = +5 dB, self-consistent with -29
in 2500 Hz.

### Our numbers against the published ones

| | published | ours | delta |
|---|---|---|---|
| WSPR, Deepest, AWGN | -31 | -30.9 | +0.1 |
| FST4W-120, AWGN | -32.8 | -32.2 | +0.6 |
| WSPR to FST4W-120 gap | 1.8 (table) / "about 1.4" (Quick-Start prose) | 1.3 | |

Our WSPR number agrees with the official number to 0.1 dB. Our measured
WSPR-to-FST4W gap of 1.3 dB brackets nicely against the Quick-Start Guide's
prose claim of "about 1.4 dB". That is a good validation of the harness.

The one thing worth a second look is our FST4W-120 at -32.2 against a published
-32.8. FST4W has a CRC, so its OSD is *not* hash-gated the way WSPR's is, and it
should reach the published figure without any hash table priming. Candidate
explanations: FST4W decode depth setting in our harness, or the published number
being slightly optimistic. Worth one experiment.

Note also that FST4-120 (the QSO mode, 101 information bits) is published at
-31.3 while FST4W-120 (74 information bits) is -32.8; the 1.5 dB is the cost of
the extra payload. WSPR carries 50 bits in a K=32 r=1/2 convolutional code.

### Independent third-party measurements

* **Daniel Estévez EA4GPZ**, Oct 2016, https://destevez.net/2016/10/simulating-jt-modes-how-low-can-they-get/
  WSJT-X trunk r7159, `wsprsim` AWGN, **bare `wsprd` with no depth flags**, 100
  files per point: -30 dB gave 34/100, -31 dB gave 6/100, so a 50% threshold
  around -29.7 dB. This is the most methodologically transparent third-party
  number and is the right comparator for *default-depth* wsprd. Our -30.9 at
  Deepest is 1.2 dB better, which is about what `-C 500 -o 4 -d` plus the
  2018-2020 decoder work should buy.
* **wsprtv/wsprd_limits**, https://github.com/wsprtv/wsprd_limits, 2026. Single-trial
  synthetic-wav study of `wsprd -d`, with the wav files committed so it is
  reproducible. Decode floor about -30 dB (a -32 dB signal fails), degrading to
  -27 dB with 4 Hz drift and -26 dB with 40% of the waveform missing. Drift
  tolerance +/-4 Hz even at -24 dB. Time-offset window +7 s / -4.5 s. Tone spacing
  0.9 to 2.6 Hz. Symbol rate tolerance about 1%. A +25 dB signal masks -25 dB
  signals across roughly 25 Hz. Not a 50% threshold study, but the only other
  quantified public reference for the decoder's operating envelope, and its
  floor is consistent with our -30.9.
* **Onno VK6FLAB**, Jan 2023, ~750k trials at 0.01 dB steps: 100% down to -29 dB,
  95% at -30 dB, decaying to zero at about -34 dB, implying a 50% point around
  -31.5 to -32. Decoder flags and hash-table state not stated, so not directly
  comparable.
* Compilations that repeat -28 or -29 without measurement: KP4MD, K0NR (Mar 2025),
  the Wikipedia WSPR article.

### Methodology trap worth recording

`wsprd` persists `hashtable.txt` in its data directory across invocations. Any
harness that loops `wsprd` over a set of generated files will populate the hash
table on the first success and silently enable OSD acceptance for every
subsequent file. Only `-H`, or deleting the file between runs, avoids it. Any
third-party number that does not state its hash-table discipline should be
treated as suspect.

## 4. Confirmation of things we already measured

The stock `lib/wsprd/README` states the OSD gate outright, so our finding that
deeper OSD does nothing for an unhashed callsign is the documented design, not a
bug:

> "The OSD is a complete decoder, meaning that it always returns a codeword.
> A returned codeword is considered valid only if the unpacked decode contains
> a callsign that is already in the hashtable."

In `wsprd.c` the gate is:

    ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
    if(strncmp(hashtab+ihash*13,callsign,13)==0) {
        if( (itype==1 && strncmp(loctab+ihash*5,grid,5)==0) || (itype==2) ) {
           not_decoded=0;
           osd_decode =1;
        }
    }

The gate is a false-alarm control and always has been. OSD was added to the deep
setting in Sept 2018 (`af7feaf65`) already hash-gated, and K9AN tightened it in
Feb 2020 (`aa8ad181d`, "wsprd: improve decoding and reduce the number of false
decodes") by additionally requiring, for type-1 messages, that the grid match the
grid stored against that callsign from its most recent Fano decode. So OSD has
never been able to help a callsign the decoder has not already seen. Same
motivation behind the metric bias value: commit `bf5195414` is "retune metric
bias to lower false-decode probability", so 0.45 was chosen to suppress false
decodes, not to maximise decodes.

`lib/wsprd/wsprd_stats.txt` is the only quantified before/after the project ever
published for the K9AN rewrite, and it is in decodes rather than dB: over 638
recorded wav files, the March 2013 K1JT `wsprd` produced 1451 decodes, `k9an-wsprd`
2122, and the merged `wsprd` 2190, in a fraction of the runtime. K9AN's own 2014
announcement said "performance seems to be comparable to wsprd" and the 2015
follow-up was about a 4 to 8 times speedup. **No before/after dB figure was ever
published for the 2015 rewrite**, so the "1-2 dB from the K9AN rewrite" belief is
not traceable to a primary source.

The WSJT-X 2.2.0 release notes describe the last real WSPR sensitivity work:
three decoding passes, coherent block detection over up to three symbols,
bit-by-bit normalisation of single-symbol bit metrics, "the number of decodes in
a crowded WSPR sub-band typically increases by 10 to 15%". No dB figure.

The noise blanker we ported is `lib/blanker.f90`, and upstream it is called from
exactly one place, `lib/fst4_decode.f90`. No other mode uses it and nobody else
has applied it to WSPR.

## 5. Unlanded proposals in the development archives

Searched the wsjt-devel SourceForge mailman archive and the WSJTX groups.io
lists. **I found no proposed-but-unmerged WSPR decoder sensitivity change on the
lists.** Note that this is not the same as "none exists": JTDX's multi-window
detector (section 2a) is exactly such a change, it just never went to the list.
Unmerged work lives in forks, not in the archive. The
WSPR traffic on those lists is operational: band hopping, frequency protection
ranges, AGC, spot uploading, three-character callsigns. The decoder-sensitivity
threads are about FT8 (the 2019 "FT8 decoding sensitivity: WSJT-X vs. JTDX"
thread being the notable one), not WSPR.

Caveat on method: the SourceForge mailman archive is paginated and very poorly
indexed by search engines, so this is a negative result from targeted queries
rather than an exhaustive sweep. If we want certainty, the archive would have to
be crawled page by page. Given that the git history shows no algorithmic change
since Feb 2020, and no fork has one either, I do not think that crawl is worth
doing.

I also found no publication, amateur or academic, applying coherent detection to
WSPR specifically. The observation that WSPR's tone spacing of exactly 1/T makes
the symbol-boundary phase data-independent is not novel in itself (it is the
defining property of Sunde's FSK / orthogonal CPFSK with h=1, textbook material),
but nobody appears to have exploited it in a WSPR decoder. Upstream's own
coherent work stops at the "block detection over 2 and 3 symbols" added in 2.2.0,
which is a much weaker use of the same property.

## 6. What "Deepest" actually is, for the record

`widgets/mainwindow.cpp` maps decode depth to `wsprd` flags as:

    depth 1  -qB              quick, no block detection, no shift jitter
    depth 2  -C 500 -o 4      3 passes, subtraction, block detection, OSD
    depth 3  -C 500 -o 4 -d   as above plus more candidates

`-C 500` *lowers* the Fano cycle budget from wsprd's own default of 10000. `-d`
(`more_candidates`) adds a candidate at every third FFT bin above `min_snr`
rather than only at spectral peaks, which is how signals too weak to form a peak
get a decode attempt at all. The three passes are configured at `wsprd.c:999`:
passes 1 and 2 use single-symbol detection with `maxdrift=4` and
`minsync2=0.12`; pass 3 sets `nblocksize=4` (block detection over 1, 2 and 3
symbols plus bit-by-bit metric normalisation), `maxdrift=0` and `minsync2=0.10`.
Pass 2 is skipped entirely if pass 1 decoded nothing.

## 6a. Dead hooks and hard-coded gates in our own tree

All verified by grep in `upstream/wsprd.c` (mainline, identical to ours).

* **`min_snr` is hard-coded**, `wsprd.c:1062-1063`:

      float min_snr, snr_scaling_factor;
      min_snr = pow(10.0,-8.0/10.0); //this is min snr in wspr bw

  That is -8 dB in the ~6 Hz WSPR bandwidth. Referred to 2500 Hz it is about
  -34.2 dB (10*log10(2500/5.9) = 26.3 dB). So the candidate gate sits only about
  3.3 dB below the stock -30.9 decode threshold, **and we have already spent 1.4
  dB of that headroom on coherent demod.** Worse, `wsprd.c:1071` floors
  everything below it:

      if( smspec[j] < min_snr) smspec[j]=0.1*min_snr;

  A floored bin can never be a local maximum, so a sub-threshold signal is not
  merely unranked, it is structurally incapable of becoming a candidate. As we
  push the decoder down, this gate becomes the binding constraint. Worth sweeping
  before anything else.

* **`apmask` is a dead a-priori hook.** `wsprd.c:785` does
  `apmask=calloc(162,sizeof(unsigned char))` and it is never written to before
  `wsprd.c:1362` passes it to `osdwspr_`. `osdwspr.f90` honours it. So the
  machinery for a-priori decoding over a known-callsign list exists, is wired up,
  and is fed all zeros. Given that WSPR beacons repeat the same callsign every
  two minutes for days, an AP list is a much better fit for WSPR than for a QSO
  mode. Nobody has proposed it.

* **Sync acquisition is still fully noncoherent** even though our demod no longer
  is. Worth checking whether the sync search is now the weaker half.

## 7. Things stock wsprd already supports that WSJT-X never turns on

Not from any fork, but found while reading the source, and cheap to test:

* **`-J`, the Jelinek stack-bucket sequential decoder** (`lib/wsprd/jelinek.c`,
  K9AN July 2015, after F. Jelinek, "Fast Sequential Decoding Algorithm Using a
  Stack"). Present in stock `wsprd`, default `stacksize=200000`. WSJT-X never
  passes `-J` at any depth. K9AN's own wsprcan ChangeLog made it the default
  there and described it as "more efficient than Fano decoder". At the very tight
  `-C 500` budget WSJT-X uses, a stack decoder and a Fano decoder can behave quite
  differently, so this is not the same experiment as raising `-C`. The KiwiSDR
  port exposes it as a runtime toggle, which suggests somebody thought it was
  worth having.
* **`-z`, the Fano/stack metric table bias**, default 0.45, applied at
  `wsprd.c:900` as `mettab[0][i]=round(10*(metric_tables[2][i]-bias))`. K9AN's
  ChangeLog describes it as "a single knob which can be tuned to trade off
  probability of timeout (erasures) vs probability of bit errors. -z 0.5
  corresponds to the standard bias, equal to the code rate. Smaller bias values
  result in higher path metrics, making it more likely that the Fano algorithm
  will make it to the final node within the allowed number of cycles, but also
  increasing the probability that the final result will contain an error."
  Again, a different axis from `-C`. Never exposed by the GUI. A sweep is cheap.
* **`delta`, the Fano threshold step**, hard-coded to 60 at `wsprd.c:810` and not
  exposed on the command line at all. Together with `-z` this is the pair of
  classic sequential-decoder tuning parameters; both are orthogonal to the cycle
  budget `-C` that we already swept.
* **The candidate-acceptance gates**: `minsync1=0.10`, `minsync2` (0.12, lowered
  to 0.10 on pass 3), and `minrms=52.0*(symfac/64.0)` at `wsprd.c:801-809`. These
  reject candidates before decoding is attempted. If any sensitivity is being
  lost at detection rather than at decoding, it is here.
* **Relaxing the OSD acceptance gate.** The structural reason WSPR needs the gate
  at all is that its 50-bit payload has **no CRC**. Contrast FST4W: its (240,74)
  codeword carries 50 message bits plus a 24-bit CRC, and `lib/fst4/decode240_74.f90`
  validates every belief-propagation and OSD result with `get_crc24`, so FST4W's
  OSD works for any callsign, hashed or not. WSPR has nothing equivalent, so K9AN
  used the hash table as a substitute false-alarm control.
  That said, for a type-1 message the callsign is transmitted in full, so the
  hash-table lookup is being used purely as a false-alarm filter, not because the
  callsign is unknown. `osdwspr_` already
  returns `nhardmin` and `dmin`, and the code already re-runs `fano()` on the OSD
  codeword. Replacing the hash lookup with a soft-distance threshold plus the
  existing format validity checks (callsign regex, valid grid, power from the
  legal set, which together carry real redundancy) would let OSD contribute for
  brand-new callsigns. This is the actual reason "deeper OSD does nothing", and
  the fix is the gate, not the depth. Nobody upstream or in any fork has tried it.

## Sources


* WSJT-X Improved: https://sourceforge.net/projects/wsjt-x-improved/
* Mainline wsprd: https://sourceforge.net/p/wsjt/wsjtx/ci/master/tree/lib/wsprd/
* FST4 Quick-Start Guide: https://wsjt.sourceforge.io/FST4_Quick_Start.pdf
* QST Nov 2010 K1JT/W1BW: https://wsjt.sourceforge.io/WSPR_QST_Nov_2010.pdf
* EA4GPZ: https://destevez.net/2016/10/simulating-jt-modes-how-low-can-they-get/
* PA3FWM: https://www.pa3fwm.nl/technotes/tn09b.html
* wsprd_limits: https://github.com/wsprtv/wsprd_limits
* wsprdaemon: https://github.com/rrobinett/wsprdaemon
* KiwiSDR: https://github.com/jks-prv/Beagle_SDR_GPS
* ka9q-radio: https://github.com/ka9q/ka9q-radio
* Debian patch series: https://sources.debian.org/src/wsjtx/3.0.2+dfsg-2/debian/patches/
* JTDX (the multi-window detector, `lib/wsprd/wsprd.c`): https://github.com/jtdx-project/jtdx
