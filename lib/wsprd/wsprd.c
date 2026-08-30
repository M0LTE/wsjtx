/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.
 
 File name: wsprd.c
 
 Copyright 2001-2018, Joe Taylor, K1JT
 
 Much of the present code is based on work by Steven Franke, K9AN,
 which in turn was based on earlier work by K1JT.
 
 Copyright 2014-2018, Steven Franke, K9AN
 
 License: GNU GPL v3
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>
#include <errno.h>

#include "fano.h"
#include "jelinek.h"
#include "nhash.h"
#include "wsprd_utils.h"
#include "wsprsim_utils.h"

#define max(x,y) ((x) > (y) ? (x) : (y))

extern void osdwspr_ (float [], unsigned char [], int *, unsigned char [], int *, float *);

// Possible PATIENCE options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT,
// FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
#define PATIENCE FFTW_ESTIMATE
fftwf_plan PLAN1=NULL,PLAN2=NULL,PLAN3=NULL;

unsigned char pr3[162]=
{1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0};

int printdata=0;

int g_nbpct  = 0;     /* noise blanker: blank the strongest g_nbpct% of samples */
int g_ndrop  = 1;     /* ...and this many samples after each hit                */
int g_nbauto = 0;     /* sweep the blanker and keep whatever any setting finds  */
int g_nthreads = 0;   /* 0 = one worker per logical processor                   */
short *g_rawbuf=NULL;   /* raw samples, stashed at the sweep's first read     */
size_t g_rawn=0;
float *g_oobenv=NULL;   /* magnitude of the out-of-band detection stream      */
float  g_oobsig=0.0f;   /* its robust sigma                                   */
int    g_oobmode=0;     /* 1 = sweep active: read once, reprocess thereafter  */
float  g_nbk=0.0f;      /* blanking threshold in sigmas, 0 = no blanking      */
int g_recall = 0;     /* -r: remember recent decodes, search there again     */

/* minutes on a linear scale good enough to age priors within a year */
static long rv_minutes(int yy,int mo,int dd,int hh,int mn)
{
    static const int cum[12]={0,31,59,90,120,151,181,212,243,273,304,334};
    long days = yy*365L + yy/4 + cum[mo-1] + dd;
    if( (yy%4)==0 && mo<=2 ) days--;
    return days*1440L + hh*60L + mn;
}

int g_guard  = 0;     /* -G: gate bare type-2 Fano decodes on symbol agreement */
int g_keepdt = 0;     /* 1 = an unpackable OSD result fails just that attempt,
                         instead of abandoning the remaining DT offsets       */
float g_minsync1 = -1.0f; /* -S overrides for the sync gates; negative values */
float g_minsync2 = -1.0f; /* leave the stock thresholds alone                 */

/* Per-symbol signal and noise, taken from the four tone energies.  For each
   symbol the sync vector allows two tones and rules out the other two; the
   ruled-out pair is pure noise, and the allowed pair carries the signal plus
   the same noise again.  With e_used and e_unused summed per symbol:
       mean(e_unused) = 4 sigma^2                (two noise-only tones)
       mean(e_used)   = A_i^2 + 4 sigma^2
   sigma is estimated once across the whole transmission, because flat fading
   moves the signal and leaves the noise where it is; A_i is estimated per
   symbol and smoothed over a few symbols to keep it from being dominated by
   its own two-sample noise.  Scaling soft symbol i by (A_i/sigma)^wexp is
   maximal-ratio combining at wexp=1; wexp=0.5 measured best and is what the
   weighted trials use. */
static void persym_weights(const double *e_used, const double *e_unused,
                           double *w, int nsym, int win, double wexp)
{
    int i,k;
    double s2=0.0;
    for(i=0;i<nsym;i++) s2 += e_unused[i];
    s2 /= (4.0*nsym);
    if(s2<=0.0){ for(i=0;i<nsym;i++) w[i]=1.0; return; }
    int half=win/2;
    for(i=0;i<nsym;i++){
        double acc=0.0; int n=0;
        for(k=i-half;k<=i+half;k++){
            if(k<0||k>=nsym) continue;
            acc += e_used[k]-4.0*s2;
            n++;
        }
        double a2 = n? acc/n : 0.0;
        if(a2<0.0) a2=0.0;
        w[i] = pow(sqrt(a2/s2), wexp);
    }
    /* keep the average weight at unity so the global scaling downstream, and
       the metric table it feeds, see the same overall level as before */
    double m=0.0;
    for(i=0;i<nsym;i++) m+=w[i];
    m/=nsym;
    if(m>0.0) for(i=0;i<nsym;i++) w[i]/=m;
}

/* wall clock, because CPU-seconds stop meaning much once we thread */
static double wsecs(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

//***************************************************************************
unsigned long readc2file(char *ptr_to_infile, float *idat, float *qdat,
                         double *freq, int *wspr_type)
{
    float *buffer;
    double dfreq;
    int i,ntrmin;
    char c2file[15];
    size_t nr;
    FILE* fp;
    
    fp = fopen(ptr_to_infile,"rb");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open data file '%s': %s\n", ptr_to_infile, strerror(errno));
        return 1;
    }
    nr=fread(c2file,sizeof(char),14,fp);
    nr=fread(&ntrmin,sizeof(int),1,fp);
    nr=fread(&dfreq,sizeof(double),1,fp);
    *freq=dfreq;
    
    buffer=calloc(2*65536,sizeof(float));
    nr=fread(buffer,sizeof(float),2*45000,fp);
    fclose(fp);
    
    *wspr_type=ntrmin;
    
    for(i=0; i<45000; i++) {
        idat[i]=buffer[2*i];
        qdat[i]=-buffer[2*i+1];
    }
    free(buffer);
    
    if( nr == 2*45000 ) {
        return (unsigned long) nr/2;
    } else {
        return 1;
    }
}

/* -------------------------------------------------------------------------
   Frequency-domain interference excision.

   readwavfile() already takes one big forward FFT of the whole recording, so
   the 375 Hz slice it keeps around 1500 Hz is sitting in memory as 46080 bins
   of 0.00814 Hz each.  That is the natural place to deal with a narrowband
   intruder: a CW carrier, a birdie, a switching-supply spur.

   Two facts make it work at that resolution:

     - a carrier is two bins wide, and 60..80 dB above the per-bin noise floor
     - a WSPR transmission is not narrowband here at all.  Its 162 symbols are
       0.68 s each, so every symbol is 1.5 Hz wide, and the four tones smear
       into a 6 Hz plateau roughly 740 bins across.  No single bin of it ever
       stands far above its own neighbours.

   So the test is not "is this bin big" -- a loud WSPR signal is also big -- but
   "is this bin big COMPARED WITH the 0.4 to 2.3 Hz either side of it".  That
   ratio is ~55 dB for a carrier and ~15 dB for a WSPR signal however strong the
   WSPR signal is, because the comparison band is inside its own plateau.

   Having found a spur, the notch has to be grown outwards, because a 114 s
   rectangular record gives a carrier a 1/df^2 leakage skirt that is genuinely
   present in the DFT coefficients and that runs many Hz either side.  Zeroing
   only the peak leaves most of the interference behind.  Growth stops where the
   smoothed spectrum falls back to the local floor, which makes the notch width
   set itself from the carrier strength instead of from a constant.

   What is left after a hard notch of width B is not a weak carrier spread over
   the record; it is the turn-on and turn-off transient of the record edges,
   about 1/B seconds long at each end.  That is why a notch a couple of Hz wide
   is worth far more than its energy fraction suggests.

   The local floor is the MINIMUM of the 256-bin block medians within +-4
   blocks, not the median of them.  A floor that followed the leakage skirt
   would stop the notch growing out of the skirt that produced it.
   ------------------------------------------------------------------------- */

int   g_fex      = 0;       /* off unless asked for: it can remove a signal  */
int   g_fexauto  = 0;       /* sweep the aggressiveness, keep every decode  */
float g_fexT1    = 1000.0f; /* peak / local floor, to be a spur at all      */
float g_fexT2    = 300.0f;  /* peak / max(side medians), the narrowness test*/
float g_fexT3    = 2.0f;    /* smoothed / floor, where notch growth stops   */
int   g_fexW     = 2048;    /* hard cap on notch half width, in bins        */
int   g_fexclamp = 0;       /* 1 = clamp to the floor, 0 = zero the bins    */
long  g_fexbins  = 0;       /* diagnostics for the last readwavfile()       */
int   g_fexhits  = 0;
float g_fexf[16];           /* frequencies of the spurs found, Hz re 1500   */

#define FEX_BS   256        /* floor block, 2.08 Hz                         */
#define FEX_MB   4          /* floor = min of block medians over +-MB       */
#define FEX_SM   33         /* boxcar width for the growth test             */
#define FEX_G1   48         /* narrowness test window, inner edge (0.39 Hz) */
#define FEX_G2   288        /* ... and outer edge (2.34 Hz)                 */
#define FEX_ITER 64         /* give up after this many spurs.  A keyed CW
                               signal is a comb of lines at the keying rate,
                               not one line, so the budget has to be generous */

static int fexcmp(const void *a, const void *b)
{
    float x=*(const float*)a, y=*(const float*)b;
    return (x<y) ? -1 : (x>y);
}
static inline int fexwrap(int i, int n){ i%=n; return (i<0)? i+n : i; }

static float fexmed(const float *p, int n, int lo, int len, float *scratch)
{
    for(int i=0;i<len;i++) scratch[i]=p[fexwrap(lo+i,n)];
    qsort(scratch,len,sizeof(float),fexcmp);
    return scratch[len/2];
}

static void excise_spurs(fftwf_complex *X, int n, double binhz)
{
    if( !g_fex || n < 8*FEX_BS ) return;
    float *P  = (float*)malloc((size_t)n*sizeof(float));
    float *Ps = (float*)malloc((size_t)n*sizeof(float));
    float *F  = (float*)malloc((size_t)n*sizeof(float));
    unsigned char *mk = (unsigned char*)calloc((size_t)n,1);
    float *sc = (float*)malloc((size_t)(FEX_G2+FEX_BS+8)*sizeof(float));
    int nb = n/FEX_BS;
    float *bm = (float*)malloc((size_t)nb*sizeof(float));
    g_fexbins=0; g_fexhits=0;
    if( !P || !Ps || !F || !mk || !sc || !bm ) goto done;

    for(int i=0;i<n;i++) P[i]=X[i][0]*X[i][0]+X[i][1]*X[i][1];

    for(int b=0;b<nb;b++) bm[b]=fexmed(P,n,b*FEX_BS,FEX_BS,sc);
    for(int b=0;b<nb;b++){
        float m=1e30f;
        for(int d=-FEX_MB;d<=FEX_MB;d++){
            float v=bm[fexwrap(b+d,nb)];
            if(v<m) m=v;
        }
        m /= 0.6931f;                 /* median of an exponential -> its mean */
        if(m<=0.0f) m=1e-30f;
        for(int i=0;i<FEX_BS;i++) F[b*FEX_BS+i]=m;
    }

    {   double s=0.0;
        for(int d=-(FEX_SM/2); d<=FEX_SM/2; d++) s+=P[fexwrap(d,n)];
        for(int i=0;i<n;i++){
            Ps[i]=(float)(s/FEX_SM);
            s += P[fexwrap(i+1+FEX_SM/2,n)] - P[fexwrap(i-FEX_SM/2,n)];
        }
    }

    for(int it=0; it<FEX_ITER; it++){
        int best=-1; float bestr=0.0f;
        for(int i=0;i<n;i++){
            if(mk[i]) continue;
            float r=P[i]/F[i];
            if(r>bestr){ bestr=r; best=i; }
        }
        if( best<0 || bestr < g_fexT1 ) break;

        float L=fexmed(P,n,best-FEX_G2,FEX_G2-FEX_G1,sc);
        float R=fexmed(P,n,best+FEX_G1,FEX_G2-FEX_G1,sc);
        float N=(L>R)?L:R;
        if(N < F[best]) N=F[best];
        if( getenv("WSPRD_FEX_DEBUG") )
            fprintf(stderr,"fex:   cand %+8.3f Hz  P/F=%10.3g  P/N=%10.3g  %s\n",
                    (best<=n/2? best : best-n)*binhz, bestr, P[best]/N,
                    (P[best] >= g_fexT2*N)?"EXCISE":"broad, left alone");
        if( P[best] < g_fexT2*N ){
            /* Broad.  A WSPR signal, or an interferer wide enough that there is
               nothing to excise.  Retire the neighbourhood and look elsewhere. */
            for(int d=-FEX_G2; d<=FEX_G2; d++) mk[fexwrap(best+d,n)] |= 2;
            continue;
        }

        int wlo=0, whi=0;
        while(wlo<g_fexW){ int j=fexwrap(best-wlo-1,n); if(Ps[j] > g_fexT3*F[j]) wlo++; else break; }
        while(whi<g_fexW){ int j=fexwrap(best+whi+1,n); if(Ps[j] > g_fexT3*F[j]) whi++; else break; }
        for(int d=-wlo; d<=whi; d++){
            int j=fexwrap(best+d,n);
            if(!(mk[j]&1)){ mk[j]|=1; g_fexbins++; }
            mk[j]|=2;
        }
        if(g_fexhits<16) g_fexf[g_fexhits]=(float)((best<=n/2? best : best-n)*binhz);
        g_fexhits++;
    }

    for(int i=0;i<n;i++){
        if( !(mk[i]&1) ) continue;
        if( g_fexclamp ){
            if(P[i]>F[i]){ float g=sqrtf(F[i]/P[i]); X[i][0]*=g; X[i][1]*=g; }
        } else {
            X[i][0]=0.0f; X[i][1]=0.0f;
        }
    }

    if( getenv("WSPRD_FEX_DEBUG") ){
        fprintf(stderr,"fex: T1=%.0f T2=%.0f T3=%.2f W=%d clamp=%d  spurs=%d bins=%ld (%.2f%%)",
                g_fexT1,g_fexT2,g_fexT3,g_fexW,g_fexclamp,g_fexhits,g_fexbins,100.0*g_fexbins/n);
        for(int i=0;i<g_fexhits && i<16;i++) fprintf(stderr," %+.2fHz",g_fexf[i]);
        fprintf(stderr,"\n");
    }
done:
    free(P); free(Ps); free(F); free(mk); free(sc); free(bm);
}

/* Hard noise blanker, the same rule as lib/blanker.f90 which FST4 and FST4W
   already use.  Find the amplitude exceeded by npct/ndropmax of the samples,
   zero those and the ndropmax samples that follow. */
static void blank_impulses(short *x, size_t nz, int npct, int ndropmax)
{
    if(npct<=0) return;
    static int hist[32769];
    memset(hist,0,sizeof hist);
    for(size_t i=0;i<nz;i++){
        if(x[i]==-32768) x[i]=-32767;
        hist[abs(x[i])]++;
    }
    long n=0, want=(long)(nz*0.01*npct/ndropmax+0.5), nthresh=0;
    for(int i=32768;i>=0;i--){ n+=hist[i]; if(n>=want){ nthresh=i; break; } }
    int ndrop=0;
    for(size_t i=0;i<nz;i++){
        if(ndrop>0){ x[i]=0; ndrop--; continue; }
        if(abs(x[i])>nthresh){ x[i]=0; ndrop=ndropmax; }
    }
}

//***************************************************************************
unsigned long readwavfile(char *ptr_to_infile, int ntrmin, float *idat, float *qdat )
{
    size_t i, j, npoints, nr;
    int nfft1, nfft2, nh2, i0;
    double df;
    
    nfft2=46080; //this is the number of downsampled points that will be returned
    nh2=nfft2/2;
    
    if( ntrmin == 2 ) {
        nfft1=nfft2*32;      //need to downsample by a factor of 32
        df=12000.0/nfft1;
        i0=1500.0/df+0.5;
        npoints=114*12000;
    } else if ( ntrmin == 15 ) {
        nfft1=nfft2*8*32;
        df=12000.0/nfft1;
        i0=(1500.0+112.5)/df+0.5;
        npoints=8*114*12000;
    } else {
        fprintf(stderr,"This should not happen\n");
        return 1;
    }
    
    float *realin;
    fftwf_complex *fftin, *fftout;
    
    FILE *fp;
    short int *buf2;
    
    fp = fopen(ptr_to_infile,"rb");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open data file '%s'\n", ptr_to_infile);
        return 1;
    }
    
    buf2 = calloc(npoints,sizeof(short int));
    if( g_oobmode && g_rawbuf != NULL ) {
        /* the sweep reprocesses the same recording; no need to read it again */
        memcpy(buf2,g_rawbuf,g_rawn*sizeof(short int));
        nr=g_rawn;
        fclose(fp);
    } else {
        nr=fread(buf2,2,22,fp);      //Read and ignore header
        nr=fread(buf2,2,npoints,fp); //Read raw data
        fclose(fp);
        if( g_oobmode && nr>0 ) {
            g_rawbuf=malloc(npoints*sizeof(short int));
            memcpy(g_rawbuf,buf2,npoints*sizeof(short int));
            g_rawn=nr;
        }
    }
    if( g_oobmode ) {
        /* Blank where the OUT-OF-BAND detection stream is impulsive.  The
           stream cannot contain the WSPR signals, so this can never mutilate
           them, which is what made the raw percentile blanker fabricate
           callsigns on a clean crowded band. */
        if( g_nbk>0.0f && g_oobenv!=NULL && nr>0 ) {
            float th=g_nbk*g_oobsig;
            for(i=0; i<npoints; i++){
                if( g_oobenv[i]>th ){
                    buf2[i]=0;
                    if(i>0) buf2[i-1]=0;
                    if(i+1<npoints) buf2[i+1]=0;
                }
            }
        }
    } else {
        if(nr>0) blank_impulses(buf2,nr,g_nbpct,g_ndrop);
    }
    if( nr == 0 ) {
        free(buf2);
        fprintf(stderr, "No data in file '%s'\n", ptr_to_infile);
        return 1;
    }	
    
    /* re-entrant: the blanker sweep reads the file more than once */
    if(PLAN1){ fftwf_destroy_plan(PLAN1); PLAN1=NULL; }
    if(PLAN2){ fftwf_destroy_plan(PLAN2); PLAN2=NULL; }
    realin=(float*) fftwf_malloc(sizeof(float)*nfft1);
    fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*(nfft1/2+1));
    PLAN1 = fftwf_plan_dft_r2c_1d(nfft1, realin, fftout, PATIENCE);
    
    for (i=0; i<npoints; i++) {
        realin[i]=buf2[i]/32768.0;
    }
    
    for (i=npoints; i<(size_t)nfft1; i++) {
        realin[i]=0.0;
    }
    free(buf2);
    
    fftwf_execute(PLAN1);
    fftwf_free(realin);

    if( g_oobmode && g_oobenv==NULL ) {
        /* Build the detection stream once, from the unblanked spectrum:
           300..5500 Hz minus a generous 1200..1800 Hz around the WSPR band. */
        fftwf_complex *oc=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*(nfft1/2+1));
        memcpy(oc,fftout,sizeof(fftwf_complex)*(nfft1/2+1));
        int jlo=(int)(300.0/df), jhi=(int)(5500.0/df);
        int klo=(int)(1200.0/df), khi=(int)(1800.0/df);
        int jj;
        for(jj=0;jj<=nfft1/2;jj++){
            if( jj<jlo || jj>jhi || (jj>=klo && jj<=khi) ){ oc[jj][0]=0.0f; oc[jj][1]=0.0f; }
        }
        float *ox=(float*)fftwf_malloc(sizeof(float)*nfft1);
        fftwf_plan op=fftwf_plan_dft_c2r_1d(nfft1,oc,ox,FFTW_ESTIMATE);
        fftwf_execute(op);
        fftwf_destroy_plan(op);
        fftwf_free(oc);
        g_oobenv=malloc(npoints*sizeof(float));
        for(i=0;i<npoints;i++) g_oobenv[i]=fabsf(ox[i]);
        fftwf_free(ox);
        long ns=0;
        float *sm=malloc((npoints/8+1)*sizeof(float));
        for(i=0;i<npoints;i+=8) sm[ns++]=g_oobenv[i];
        qsort(sm,ns,sizeof(float),floatcomp);
        g_oobsig=1.4826f*sm[ns/2];
        free(sm);
        if(g_oobsig<=0.0f) g_oobsig=1.0f;
    }
    
    fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nfft2);
    
    for (i=0; i<(size_t)nfft2; i++) {
        j=i0+i;
        if( i>(size_t)nh2 ) j=j-nfft2;
        fftin[i][0]=fftout[j][0];
        fftin[i][1]=fftout[j][1];
    }
    
    fftwf_free(fftout);

    /* Excise narrowband interference from the slice before it is transformed
       back.  Everything downstream, the noise-level percentile, the candidate
       list, and above all the totp normalisation inside sync_and_demodulate,
       then sees a band without the intruder in it. */
    excise_spurs(fftin, nfft2, df);

    fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nfft2);
    PLAN2 = fftwf_plan_dft_1d(nfft2, fftin, fftout, FFTW_BACKWARD, PATIENCE);
    fftwf_execute(PLAN2);
    
    for (i=0; i<(size_t)nfft2; i++) {
        idat[i]=fftout[i][0]/1000.0;
        qdat[i]=fftout[i][1]/1000.0;
    }
    
    fftwf_free(fftin);
    fftwf_free(fftout);
    return nfft2;
}

//***************************************************************************
void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode)
{
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/
    
    float fplast=-10000.0;   /* was static: shared mutable state, not thread safe */
    static const float dt=1.0/375.0, df=375.0/256.0;
    static const float pi=3.14159265358979323846;
    float twopidt, df15=df*1.5, df05=df*0.5;
    
    int i, j, k, lag;
    float i0[162],q0[162],i1[162],q1[162],i2[162],q2[162],i3[162],q3[162];
    float p0,p1,p2,p3,cmet,syncmax,fac;
    float c0[256],s0[256],c1[256],s1[256],c2[256],s2[256],c3[256],s3[256];
    float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float f0=0.0, fp, ss, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[162];
    int best_shift = 0, ifreq;
    
    syncmax=-1e30;
    if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0; f0=*f1;}
    if( mode == 1 ) {lagmin=*shift1;lagmax=*shift1;f0=*f1;}
    if( mode == 2 ) {lagmin=*shift1;lagmax=*shift1;ifmin=0;ifmax=0;f0=*f1;}
    
    twopidt=2*pi*dt;
    /* The tone tables below depend on the frequency, the drift and the symbol
       index but not on the lag, so the symbol loop goes outside the lag loop and
       each table is built once and used for every lag. */
#define MAXLAGS 512
    int nlags, ilag;
    float sslag[MAXLAGS], totplag[MAXLAGS];
    for(ifreq=ifmin; ifreq<=ifmax; ifreq++) {
        f0=*f1+ifreq*fstep;
        nlags=0;
        for(lag=lagmin; lag<=lagmax && nlags<MAXLAGS; lag=lag+lagstep) nlags++;
        for(ilag=0; ilag<nlags; ilag++){ sslag[ilag]=0.0; totplag[ilag]=0.0; }
        {
            for (i=0; i<162; i++) {
                fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
                if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
                    dphi0=twopidt*(fp-df15);
                    cdphi0=cos(dphi0);
                    sdphi0=sin(dphi0);
                    
                    dphi1=twopidt*(fp-df05);
                    cdphi1=cos(dphi1);
                    sdphi1=sin(dphi1);
                    
                    dphi2=twopidt*(fp+df05);
                    cdphi2=cos(dphi2);
                    sdphi2=sin(dphi2);
                    
                    dphi3=twopidt*(fp+df15);
                    cdphi3=cos(dphi3);
                    sdphi3=sin(dphi3);
                    
                    c0[0]=1; s0[0]=0;
                    c1[0]=1; s1[0]=0;
                    c2[0]=1; s2[0]=0;
                    c3[0]=1; s3[0]=0;
                    
                    for (j=1; j<256; j++) {
                        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
                    }
                    fplast = fp;
                }
                
                for(ilag=0, lag=lagmin; ilag<nlags; ilag++, lag=lag+lagstep) {
                i0[i]=0.0; q0[i]=0.0;
                i1[i]=0.0; q1[i]=0.0;
                i2[i]=0.0; q2[i]=0.0;
                i3[i]=0.0; q3[i]=0.0;
                
                for (j=0; j<256; j++) {
                    k=lag+i*256+j;
                    if( (k>0) && (k<np) ) {
                        i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                        q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                        i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                        q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                        i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                        q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                        i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                        q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
                    }
                }
                p0=i0[i]*i0[i] + q0[i]*q0[i];
                p1=i1[i]*i1[i] + q1[i]*q1[i];
                p2=i2[i]*i2[i] + q2[i]*q2[i];
                p3=i3[i]*i3[i] + q3[i]*q3[i];
                
                p0=sqrt(p0);
                p1=sqrt(p1);
                p2=sqrt(p2);
                p3=sqrt(p3);
                
                totplag[ilag]=totplag[ilag]+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                sslag[ilag] = (pr3[i] == 1) ? sslag[ilag]+cmet : sslag[ilag]-cmet;
                if( mode == 2) {                 //Compute soft symbols
                    if(pr3[i]==1) {
                        fsymb[i]=p3-p1;
                    } else {
                        fsymb[i]=p2-p0;
                    }
                }
                } // lag loop
            }
        }
        for(ilag=0, lag=lagmin; ilag<nlags; ilag++, lag=lag+lagstep) {
            ss=sslag[ilag]/totplag[ilag];
            if( ss > syncmax ) {          //Save best parameters
                syncmax=ss;
                best_shift=lag;
                fbest=f0;
            }
        }
    } //freq loop
    
    if( mode <=1 ) {                       //Send best params back to caller
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }
    
    if( mode == 2 ) {
        *sync=syncmax;
        for (i=0; i<162; i++) {              //Normalize the soft symbols
            fsum=fsum+fsymb[i]/162.0;
            f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
        }
        fac=sqrt(f2sum-fsum*fsum);
        for (i=0; i<162; i++) {
            fsymb[i]=symfac*fsymb[i]/fac;
            if( fsymb[i] > 127) fsymb[i]=127.0;
            if( fsymb[i] < -128 ) fsymb[i]=-128.0;
            symbols[i]=fsymb[i] + 128;
        }
        return;
    }
    return;
}

void noncoherent_sequence_detection(float *id, float *qd, long np,
                                    unsigned char *symbols, float *f1,  int *shift1,
                                    float *drift1, int symfac, int *nblocksize, int *bitmetric,
                                    float wexp)
{
    /************************************************************************
     *  Noncoherent sequence detection for wspr.                            *
     *  Allowed block lengths are nblock=1,2,3,6, or 9 symbols.             *
     *  Longer block lengths require longer channel coherence time.         *
     *  The whole block is estimated at once.                               *
     *  nblock=1 corresponds to noncoherent detection of individual symbols *
     *     like the original wsprd symbol demodulator.                      *
     ************************************************************************/
    float fplast=-10000.0;   /* was static: shared mutable state, not thread safe */
    static const float dt=1.0/375.0, df=375.0/256.0;
    static const float pi=3.14159265358979323846;
    float twopidt, df15=df*1.5, df05=df*0.5;
    
    int i, j, k, lag, itone, ib, b, nblock, nseq, imask;
    float xi[512],xq[512];
    float is[4][162],qs[4][162],cf[4][162],sf[4][162],cm,sm,cmp,smp;
    float p[512],fac,xm1,xm0;
    float c0[257],s0[257],c1[257],s1[257],c2[257],s2[257],c3[257],s3[257];
    float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float f0, fp, fsum=0.0, f2sum=0.0, fsymb[162];
    
    twopidt=2*pi*dt;
    f0=*f1;
    lag=*shift1;
    nblock=*nblocksize;
    nseq=1<<nblock;
    int bitbybit=*bitmetric;
    
    for (i=0; i<162; i++) {
        fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
        if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
            dphi0=twopidt*(fp-df15);
            cdphi0=cos(dphi0);
            sdphi0=sin(dphi0);
            
            dphi1=twopidt*(fp-df05);
            cdphi1=cos(dphi1);
            sdphi1=sin(dphi1);
            
            dphi2=twopidt*(fp+df05);
            cdphi2=cos(dphi2);
            sdphi2=sin(dphi2);
            
            dphi3=twopidt*(fp+df15);
            cdphi3=cos(dphi3);
            sdphi3=sin(dphi3);
            
            c0[0]=1; s0[0]=0;
            c1[0]=1; s1[0]=0;
            c2[0]=1; s2[0]=0;
            c3[0]=1; s3[0]=0;
            
            for (j=1; j<257; j++) {
                c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
            }
            
            fplast = fp;
        }
        
        cf[0][i]=c0[256]; sf[0][i]=s0[256];
        cf[1][i]=c1[256]; sf[1][i]=s1[256];
        cf[2][i]=c2[256]; sf[2][i]=s2[256];
        cf[3][i]=c3[256]; sf[3][i]=s3[256];
        
        is[0][i]=0.0; qs[0][i]=0.0;
        is[1][i]=0.0; qs[1][i]=0.0;
        is[2][i]=0.0; qs[2][i]=0.0;
        is[3][i]=0.0; qs[3][i]=0.0;
        
        for (j=0; j<256; j++) {
            k=lag+i*256+j;
            if( (k>0) && (k<np) ) {
                is[0][i]=is[0][i] + id[k]*c0[j] + qd[k]*s0[j];
                qs[0][i]=qs[0][i] - id[k]*s0[j] + qd[k]*c0[j];
                is[1][i]=is[1][i] + id[k]*c1[j] + qd[k]*s1[j];
                qs[1][i]=qs[1][i] - id[k]*s1[j] + qd[k]*c1[j];
                is[2][i]=is[2][i] + id[k]*c2[j] + qd[k]*s2[j];
                qs[2][i]=qs[2][i] - id[k]*s2[j] + qd[k]*c2[j];
                is[3][i]=is[3][i] + id[k]*c3[j] + qd[k]*s3[j];
                qs[3][i]=qs[3][i] - id[k]*s3[j] + qd[k]*c3[j];
            }
        }
    }
    
    for (i=0; i<162; i=i+nblock) {
        for (j=0;j<nseq;j++) {
            xi[j]=0.0; xq[j]=0.0;
            cm=1; sm=0;
            for (ib=0; ib<nblock; ib++) {
                b=(j&(1<<(nblock-1-ib)))>>(nblock-1-ib);
                itone=pr3[i+ib]+2*b;
                xi[j]=xi[j]+is[itone][i+ib]*cm + qs[itone][i+ib]*sm;
                xq[j]=xq[j]+qs[itone][i+ib]*cm - is[itone][i+ib]*sm;
                cmp=cf[itone][i+ib]*cm - sf[itone][i+ib]*sm;
                smp=sf[itone][i+ib]*cm + cf[itone][i+ib]*sm;
                cm=cmp; sm=smp;
            }
            p[j]=xi[j]*xi[j]+xq[j]*xq[j];
            p[j]=sqrt(p[j]);
        }
        for (ib=0; ib<nblock; ib++) {
            imask=1<<(nblock-1-ib);
            xm1=0.0; xm0=0.0;
            for (j=0; j<nseq; j++) {
                if((j & imask)!=0) {
                    if(p[j] > xm1) xm1=p[j];
                }
                if((j & imask)==0) {
                    if(p[j]>xm0) xm0=p[j];
                }
            }
            fsymb[i+ib]=xm1-xm0;
            if( bitbybit == 1 ) {
                fsymb[i+ib]=fsymb[i+ib]/(xm1 > xm0 ? xm1 : xm0);
            }
        }
    }
    if( wexp != 0.0f ) {
        double eu[162], ev[162], wt[162];
        for(i=0;i<162;i++){
            int t0=pr3[i], t1=pr3[i]+2, u0=1-pr3[i], u1=3-pr3[i];
            eu[i] = (double)is[t0][i]*is[t0][i]+(double)qs[t0][i]*qs[t0][i]
                  + (double)is[t1][i]*is[t1][i]+(double)qs[t1][i]*qs[t1][i];
            ev[i] = (double)is[u0][i]*is[u0][i]+(double)qs[u0][i]*qs[u0][i]
                  + (double)is[u1][i]*is[u1][i]+(double)qs[u1][i]*qs[u1][i];
        }
        persym_weights(eu,ev,wt,162,5,wexp);
        for(i=0;i<162;i++) fsymb[i]=(float)(fsymb[i]*wt[i]);
    }
    for (i=0; i<162; i++) {              //Normalize the soft symbols
        fsum=fsum+fsymb[i]/162.0;
        f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
    }
    fac=sqrt(f2sum-fsum*fsum);
    for (i=0; i<162; i++) {
        fsymb[i]=symfac*fsymb[i]/fac;
        if( fsymb[i] > 127) fsymb[i]=127.0;
        if( fsymb[i] < -128 ) fsymb[i]=-128.0;
        symbols[i]=fsymb[i] + 128;
    }
    return;
}


/* Small zero-padded DFT helper used by the coherent demodulator to find the
   residual carrier frequency left over after the noncoherent sync search. */
#define VNFFT 2048
static fftwf_plan VPLAN=NULL;
static fftwf_complex *vin=NULL,*vout=NULL;

static void vensure(void)
{
  if(!VPLAN){   /* call once before any threads start; see main() */
    vin =(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*VNFFT);
    vout=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*VNFFT);
    VPLAN=fftwf_plan_dft_1d(VNFFT,vin,vout,FFTW_FORWARD,FFTW_ESTIMATE);
  }
}

static double vpeak(const double *zr,const double *zi,int i0,int L,
                    const double *psi,double *ubest,
                    fftwf_complex *in,fftwf_complex *out)
{
  for(int m=0;m<VNFFT;m++){ in[m][0]=0.f; in[m][1]=0.f; }
  for(int i=0;i<L;i++){
    double c=cos(-psi[i0+i]), s=sin(-psi[i0+i]);
    in[i][0]=(float)(zr[i0+i]*c-zi[i0+i]*s);
    in[i][1]=(float)(zr[i0+i]*s+zi[i0+i]*c);
  }
  /* fftwf_execute() reuses the buffers the plan was built with, so it cannot be
     called from two threads at once.  fftwf_execute_dft() runs the same plan on
     buffers the caller supplies, which can. */
  fftwf_execute_dft(VPLAN,in,out);
  double best=0; int mb=0;
  for(int m=0;m<VNFFT;m++){
    double v=(double)out[m][0]*out[m][0]+(double)out[m][1]*out[m][1];
    if(v>best){best=v;mb=m;}
  }
  if(ubest) *ubest = (mb<VNFFT/2? (double)mb : (double)mb-VNFFT)/VNFFT;
  return best;
}

/*****************************************************************************
  Decision-directed COHERENT demodulator.

  WSPR is continuous-phase 4-FSK with tone spacing exactly 1/T.  The phase a
  symbol adds is 2*pi*f0*T + pi*(2s-3), and 2s-3 is odd for every one of the
  four tones, so every tone advances the carrier phase by the same amount
  modulo 2*pi.  The transmit phase at each symbol boundary is therefore
  independent of the data, which means a coherent phase reference can be built
  from the received signal alone -- no pilot, no known preamble.

  Stock wsprd never uses this.  Its demodulator is noncoherent: it throws the
  phase away and keeps only tone magnitudes, then optionally re-combines 2 or 3
  neighbouring symbols coherently (nblock).  Noncoherent detection of binary
  orthogonal signalling costs about 1 dB at high SNR and a great deal more at
  the error rates WSPR actually runs at.

  Here we:
    1. correlate against all four tones per symbol, as usual;
    2. rotate out the known, data-independent inter-symbol phase ramp so a
       stable channel shows a constant phase;
    3. take a tentative hard decision per symbol, fit the residual carrier
       frequency and drift to those decisions (the decoder's own estimates are
       only good to a few tenths of a hertz, which is many cycles over 110 s
       -- noncoherent detection never cared, coherent detection does);
    4. smooth the decision-directed complex amplitudes over nsmooth symbols to
       track the channel, leaving each symbol out of its own estimate so it
       cannot vote for itself;
    5. form COHERENT soft symbols, Re{z1 conj(h)} - Re{z0 conj(h)}, which is
       also maximal-ratio weighted across a fading channel for free.

  nsmooth is the one knob: long windows win on a stable path, short ones keep
  working when the path is fading.  Trying several is cheap.
 *****************************************************************************/
void coherent_sequence_detection(float *id, float *qd, long np,
                                 unsigned char *symbols, float *f1, int *shift1,
                                 float *drift1, int symfac, int nsmooth,
                                 float wexp)
{
    const double dt=1.0/375.0, df=375.0/256.0;
    const double pi=3.14159265358979323846, twopidt=2*pi*dt;
    const double df15=df*1.5, df05=df*0.5, tsym=256.0*dt;
    const int nsym=162;
    int i,j,k,t;
    float is[4][162],qs[4][162];
    double cfs[162],sfs[162];
    double c0[257],s0[257],c1[257],s1[257],c2[257],s2[257],c3[257],s3[257];
    double Zr[4][162],Zi[4][162];
    double wr[162],wi[162],psi[162];
    double f0=*f1, drift=*drift1;
    int lag=*shift1;

    for (i=0; i<nsym; i++) {
        double fp = f0 + (drift/2.0)*((double)i-81.0)/81.0;
        double d0=twopidt*(fp-df15), d1=twopidt*(fp-df05);
        double d2=twopidt*(fp+df05), d3=twopidt*(fp+df15);
        double cd0=cos(d0),sd0=sin(d0),cd1=cos(d1),sd1=sin(d1);
        double cd2=cos(d2),sd2=sin(d2),cd3=cos(d3),sd3=sin(d3);
        c0[0]=1;s0[0]=0;c1[0]=1;s1[0]=0;c2[0]=1;s2[0]=0;c3[0]=1;s3[0]=0;
        for (j=1;j<257;j++){
            c0[j]=c0[j-1]*cd0-s0[j-1]*sd0; s0[j]=c0[j-1]*sd0+s0[j-1]*cd0;
            c1[j]=c1[j-1]*cd1-s1[j-1]*sd1; s1[j]=c1[j-1]*sd1+s1[j-1]*cd1;
            c2[j]=c2[j-1]*cd2-s2[j-1]*sd2; s2[j]=c2[j-1]*sd2+s2[j-1]*cd2;
            c3[j]=c3[j-1]*cd3-s3[j-1]*sd3; s3[j]=c3[j-1]*sd3+s3[j-1]*cd3;
        }
        /* the per-symbol phase advance: identical for all four tones, which is
           the whole reason this works.  Take it from tone 0. */
        cfs[i]=c0[256]; sfs[i]=s0[256];
        for(t=0;t<4;t++){ is[t][i]=0.0f; qs[t][i]=0.0f; }
        for (j=0;j<256;j++){
            k=lag+i*256+j;
            if(k>0 && k<np){
                float xi=id[k], xq=qd[k];
                is[0][i]+= xi*c0[j]+xq*s0[j];  qs[0][i]+= xq*c0[j]-xi*s0[j];
                is[1][i]+= xi*c1[j]+xq*s1[j];  qs[1][i]+= xq*c1[j]-xi*s1[j];
                is[2][i]+= xi*c2[j]+xq*s2[j];  qs[2][i]+= xq*c2[j]-xi*s2[j];
                is[3][i]+= xi*c3[j]+xq*s3[j];  qs[3][i]+= xq*c3[j]-xi*s3[j];
            }
        }
    }

    /* rotate out the cumulative inter-symbol phase ramp */
    double pr=1.0, pq=0.0;
    for (i=0;i<nsym;i++){
        for(t=0;t<4;t++){
            Zr[t][i]= is[t][i]*pr + qs[t][i]*pq;
            Zi[t][i]= qs[t][i]*pr - is[t][i]*pq;
        }
        double a=pr*cfs[i]-pq*sfs[i], b=pr*sfs[i]+pq*cfs[i];
        double m=sqrt(a*a+b*b); if(m>0){a/=m;b/=m;}
        pr=a; pq=b;
    }

    /* tentative decisions -> decision-directed complex amplitudes */
    for (i=0;i<nsym;i++){
        int t0=pr3[i], t1=pr3[i]+2;
        double p0=Zr[t0][i]*Zr[t0][i]+Zi[t0][i]*Zi[t0][i];
        double p1=Zr[t1][i]*Zr[t1][i]+Zi[t1][i]*Zi[t1][i];
        int tb = (p1>p0)? t1 : t0;
        wr[i]=Zr[tb][i]; wi[i]=Zi[tb][i];
    }

    /* fit residual frequency and drift to those decisions */
    fftwf_complex vfi[VNFFT], vfo[VNFFT];   /* per-call FFT buffers */
    {
        double bdd=0.0,bu=0.0,bval=-1.0;
        for(int stage=0;stage<2;stage++){
            double lo,hi,st;
            if(stage==0){lo=-1.2;hi=1.2;st=0.01;} else {lo=bdd-0.012;hi=bdd+0.012;st=0.0015;}
            for(double dd=lo; dd<=hi+1e-12; dd+=st){
                double acc=0.0;
                for(i=0;i<nsym;i++){ psi[i]=acc;
                    acc += 2.0*pi*tsym*dd*((double)i-(double)(nsym-1)/2.0)/(double)nsym; }
                double u,v=vpeak(wr,wi,0,nsym,psi,&u,vfi,vfo);
                if(v>bval){ bval=v; bdd=dd; bu=u; }
            }
        }
        double acc=0.0;
        for(i=0;i<nsym;i++){ psi[i]=acc+2.0*pi*bu*i;
            acc += 2.0*pi*tsym*bdd*((double)i-(double)(nsym-1)/2.0)/(double)nsym; }
    }
    /* de-rotate everything by the fitted carrier */
    for(i=0;i<nsym;i++){
        double c=cos(-psi[i]), s=sin(-psi[i]);
        for(t=0;t<4;t++){
            double a=Zr[t][i]*c-Zi[t][i]*s, b=Zr[t][i]*s+Zi[t][i]*c;
            Zr[t][i]=a; Zi[t][i]=b;
        }
        double a=wr[i]*c-wi[i]*s, b=wr[i]*s+wi[i]*c;
        wr[i]=a; wi[i]=b;
    }

    /* Leave-one-out smoothed channel estimate, then two decision-directed
       iterations.  Symbols are weighted by how reliable their own decision
       looks, so the ones the demodulator is unsure about do not drag the
       channel estimate around; the second pass re-decides using the coherent
       metric from the first, which is a better decision to direct with. */
    double fsymb[162], csr[163], csi[163];
    int half = nsmooth/2;
    for(int iter=0; iter<2; iter++){
        if(iter>0){
            /* re-decide from the coherent metric of the previous iteration */
            for(i=0;i<nsym;i++){
                int t0=pr3[i], t1=pr3[i]+2;
                int tb = (fsymb[i]>0.0)? t1 : t0;
                double a0=sqrt(Zr[t0][i]*Zr[t0][i]+Zi[t0][i]*Zi[t0][i]);
                double a1=sqrt(Zr[t1][i]*Zr[t1][i]+Zi[t1][i]*Zi[t1][i]);
                double rel=fabs(a1-a0)/(a1+a0+1e-30);
                wr[i]=Zr[tb][i]*rel; wi[i]=Zi[tb][i]*rel;
            }
        } else {
            for(i=0;i<nsym;i++){
                int t0=pr3[i], t1=pr3[i]+2;
                double a0=sqrt(Zr[t0][i]*Zr[t0][i]+Zi[t0][i]*Zi[t0][i]);
                double a1=sqrt(Zr[t1][i]*Zr[t1][i]+Zi[t1][i]*Zi[t1][i]);
                int tb = (a1>a0)? t1 : t0;
                double rel=fabs(a1-a0)/(a1+a0+1e-30);
                wr[i]=Zr[tb][i]*rel; wi[i]=Zi[tb][i]*rel;
            }
        }
        csr[0]=0; csi[0]=0;
        for(i=0;i<nsym;i++){ csr[i+1]=csr[i]+wr[i]; csi[i+1]=csi[i]+wi[i]; }
        for(i=0;i<nsym;i++){
            int a=i-half, b=i+half;
            if(a<0){ b-=a; a=0; }
            if(b>nsym-1){ a-=(b-(nsym-1)); b=nsym-1; if(a<0)a=0; }
            double hr=csr[b+1]-csr[a]-wr[i];      /* leave this symbol out */
            double hi=csi[b+1]-csi[a]-wi[i];
            double hm=sqrt(hr*hr+hi*hi);
            if(hm<1e-20){ fsymb[i]=0.0; continue; }
            hr/=hm; hi/=hm;
            int t0=pr3[i], t1=pr3[i]+2;
            double r1=Zr[t1][i]*hr+Zi[t1][i]*hi;   /* Re{ z . conj(h) } */
            double r0=Zr[t0][i]*hr+Zi[t0][i]*hi;
            fsymb[i]=r1-r0;
        }
    }

    if( wexp != 0.0f ) {
        double eu[162], ev[162], wt[162];
        for(i=0;i<nsym;i++){
            int t0=pr3[i], t1=pr3[i]+2, u0=1-pr3[i], u1=3-pr3[i];
            eu[i] = Zr[t0][i]*Zr[t0][i]+Zi[t0][i]*Zi[t0][i]
                  + Zr[t1][i]*Zr[t1][i]+Zi[t1][i]*Zi[t1][i];
            ev[i] = Zr[u0][i]*Zr[u0][i]+Zi[u0][i]*Zi[u0][i]
                  + Zr[u1][i]*Zr[u1][i]+Zi[u1][i]*Zi[u1][i];
        }
        persym_weights(eu,ev,wt,nsym,5,wexp);
        for(i=0;i<nsym;i++) fsymb[i]*=wt[i];
    }
    double fsum=0.0,f2sum=0.0;
    for(i=0;i<nsym;i++){ fsum+=fsymb[i]/nsym; f2sum+=fsymb[i]*fsymb[i]/nsym; }
    double fac=sqrt(f2sum-fsum*fsum);
    if(fac<=0.0) fac=1.0;
    for(i=0;i<nsym;i++){
        double v=symfac*fsymb[i]/fac;
        if(v>127.0) v=127.0;
        if(v<-128.0) v=-128.0;
        symbols[i]=(unsigned char)(v+128.5);
    }
}

/***************************************************************************
 symbol-by-symbol signal subtraction
 ****************************************************************************/
void subtract_signal(float *id, float *qd, long np,
                     float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
    float dt=1.0/375.0, df=375.0/256.0;
    int i, j, k;
    float pi=4.*atan(1.0),twopidt, fp;
    
    float i0,q0;
    float c0[256],s0[256];
    float dphi, cdphi, sdphi;
    
    twopidt=2*pi*dt;
    
    for (i=0; i<162; i++) {
        fp = f0 + ((float)drift0/2.0)*((float)i-81.0)/81.0;
        
        dphi=twopidt*(fp+((float)channel_symbols[i]-1.5)*df);
        cdphi=cos(dphi);
        sdphi=sin(dphi);
        
        c0[0]=1; s0[0]=0;
        
        for (j=1; j<256; j++) {
            c0[j]=c0[j-1]*cdphi - s0[j-1]*sdphi;
            s0[j]=c0[j-1]*sdphi + s0[j-1]*cdphi;
        }
        
        i0=0.0; q0=0.0;
        
        for (j=0; j<256; j++) {
            k=shift0+i*256+j;
            if( (k>0) & (k<np) ) {
                i0=i0 + id[k]*c0[j] + qd[k]*s0[j];
                q0=q0 - id[k]*s0[j] + qd[k]*c0[j];
            }
        }
        
        
        // subtract the signal here.
        
        i0=i0/256.0; //will be wrong for partial symbols at the edges...
        q0=q0/256.0;
        
        for (j=0; j<256; j++) {
            k=shift0+i*256+j;
            if( (k>0) & (k<np) ) {
                id[k]=id[k]- (i0*c0[j] - q0*s0[j]);
                qd[k]=qd[k]- (q0*c0[j] + i0*s0[j]);
            }
        }
    }
    return;
}
/******************************************************************************
  Subtract the coherent component of a signal 
 *******************************************************************************/
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
    float dt=1.0/375.0, df=375.0/256.0;
    float pi=4.*atan(1.0), twopidt, phi=0, dphi, cs;
    int i, j, k, ii, nsym=162, nspersym=256,  nfilt=360; //nfilt must be even number.
    int nsig=nsym*nspersym;
    int nc2=45000;
    
    float *refi, *refq, *ci, *cq, *cfi, *cfq;
    
    refi=calloc(nc2,sizeof(float));
    refq=calloc(nc2,sizeof(float));
    ci=calloc(nc2,sizeof(float));
    cq=calloc(nc2,sizeof(float));
    cfi=calloc(nc2,sizeof(float));
    cfq=calloc(nc2,sizeof(float));
    
    twopidt=2.0*pi*dt;
    
    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/
    
    // create reference wspr signal vector, centered on f0.
    //
    for (i=0; i<nsym; i++) {
        
        cs=(float)channel_symbols[i];
        
        dphi=twopidt*
        (
         f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
         + (cs-1.5)*df
         );
        
        for ( j=0; j<nspersym; j++ ) {
            ii=nspersym*i+j;
            refi[ii]=cos(phi); //cannot precompute sin/cos because dphi is changing
            refq[ii]=sin(phi);
            phi=phi+dphi;
        }
    }

    float w[nfilt], norm=0, partialsum[nfilt]; 
    //lowpass filter and remove startup transient
    for (i=0; i<nfilt; i++) partialsum[i]=0.0;
    for (i=0; i<nfilt; i++) {
        w[i]=sin(pi*(float)i/(float)(nfilt-1));
        norm=norm+w[i];
    }
    for (i=0; i<nfilt; i++) {
        w[i]=w[i]/norm;
    }
    for (i=1; i<nfilt; i++) {
        partialsum[i]=partialsum[i-1]+w[i];
    }

    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (i=0; i<nsym*nspersym; i++) {
        k=shift0+i;
        if( (k>0) && (k<np) ) {
            ci[i+nfilt] = id[k]*refi[i] + qd[k]*refq[i];
            cq[i+nfilt] = qd[k]*refi[i] - id[k]*refq[i];
        }
    }

    // LPF.  Called once per decode from the serial part of the search, and a
    // 360-tap convolution over 45000 points is enough work to be worth sharing.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(j)
#endif
    for (i=nfilt/2; i<45000-nfilt/2; i++) {
        cfi[i]=0.0; cfq[i]=0.0;
        for (j=0; j<nfilt; j++) {
            cfi[i]=cfi[i]+w[j]*ci[i-nfilt/2+j];
            cfq[i]=cfq[i]+w[j]*cq[i-nfilt/2+j];
        }
    }

    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq+cq*refi)
    // beginning of first symbol in reference signal is at i=nfilt
    // beginning of first symbol in received data is at shift0.
    for (i=0; i<nsig; i++) {
        if( i<nfilt/2 ) {        // take care of the end effect (LPF step response) here
            norm=partialsum[nfilt/2+i];
        } else if( i>(nsig-1-nfilt/2) ) {
            norm=partialsum[nfilt/2+nsig-1-i];
        } else {
            norm=1.0;
        }
        k=shift0+i;
        j=i+nfilt;
        if( (k>0) && (k<np) ) {
            id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
            qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
        }
    }
    
    free(refi);
    free(refq);
    free(ci);
    free(cq);
    free(cfi);
    free(cfq);
    
    return;
}

unsigned long writec2file(char *c2filename, int trmin, double freq
                          , float *idat, float *qdat)
{
    int i;
    float *buffer;
    FILE *fp;
    
    fp = fopen(c2filename,"wb");
    if( fp == NULL ) {
        fprintf(stderr, "Could not open c2 file '%s'\n", c2filename);
        return 0;
    }
    unsigned long nwrite = fwrite(c2filename,sizeof(char),14,fp);
    nwrite = fwrite(&trmin, sizeof(int), 1, fp);
    nwrite = fwrite(&freq, sizeof(double), 1, fp);
    
    buffer=calloc(2*45000,sizeof(float));
    for(i=0; i<45000; i++) {
        buffer[2*i]=idat[i];
        buffer[2*i+1]=-qdat[i];
    }
    nwrite = fwrite(buffer, sizeof(float), 2*45000, fp);
    free(buffer);
    fclose(fp);
    if( nwrite == 2*45000 ) {
        return nwrite;
    } else {
        return 0;
    }
}

unsigned int count_hard_errors( unsigned char *symbols, unsigned char *channel_symbols)
{
    int i,is;
    unsigned char cw[162];
    unsigned int nerrors;
    for (i=0; i<162; i++) {
        cw[i] = channel_symbols[i] >=2 ? 1:0;
    }
    deinterleave(cw);
    nerrors=0;
    for (i=0; i<162; i++) {
        is = symbols[i] > 127 ? 1:0;
        nerrors = nerrors + (is == cw[i] ? 0:1);
    }
    return nerrors;
}


/*****************************************************************************
  One demodulate-and-decode attempt, made self-contained so the search grid can
  be run in parallel.

  For each candidate wsprd tries every combination of demodulator setting and
  DT jitter offset until one decodes.  Those attempts are completely independent
  of each other -- they only read the sample buffers -- so on a quiet frequency,
  where every attempt fails and all of them therefore run, the grid is
  embarrassingly parallel.  That is where nearly all of the decoder's time goes.

  Results stay deterministic because the parallel loop keeps the lowest-numbered
  attempt that succeeded, which is exactly the one the sequential search would
  have stopped at.
 *****************************************************************************/
struct trial_result {
    int decoded;
    int abandoned;    /* OSD produced something unpackable: stock gives up on
                         the rest of this setting's DT offsets here          */
    int blocksize, bitmetric, jitter, shift;
    unsigned int metric, cycles;
    int nhardmin, osd_decode;
    float dmin;
    unsigned char decdata[11];
    unsigned char symbols[162];
};

struct dec_ctx {
    float *idat, *qdat;
    long  npoints;
    int   (*mettab)[256];
    int   delta, ndepth, symfac, stackdecoder;
    unsigned int maxcycles, nbits, stacksize;
    float minrms;
    const char *hashtab, *loctab;
};

static void run_trial(const struct dec_ctx *c, float f1, int shift1, float drift1,
                      int blocksize, int bitmetric, int jitter, float wexp,
                      struct snode *stack, struct trial_result *out)
{
    unsigned char symbols[162*2], decdata[11], cw[162], apmask[162];
    signed char message[11];
    char callsign[13], grid[5];
    float fsymbs[162], dmin=0.0;
    unsigned int metric=0, cycles=0, maxnp=0;
    int nhardmin=0, i, not_decoded=1, osd_decode=0;
    int jittered_shift=shift1+jitter;

    memset(out,0,sizeof *out);
    out->blocksize=blocksize; out->bitmetric=bitmetric;
    out->jitter=jitter; out->shift=jittered_shift;
    memset(symbols,0,sizeof symbols);
    memset(apmask,0,sizeof apmask);

    if( blocksize < 0 ) {
        coherent_sequence_detection(c->idat, c->qdat, c->npoints, symbols, &f1,
                                    &jittered_shift, &drift1, c->symfac, -blocksize,
                                    wexp);
    } else {
        noncoherent_sequence_detection(c->idat, c->qdat, c->npoints, symbols, &f1,
                                       &jittered_shift, &drift1, c->symfac,
                                       &blocksize, &bitmetric, wexp);
    }

    float sq=0.0;
    for(i=0; i<162; i++){ float y=(float)symbols[i]-128.0f; sq += y*y; }
    if( sqrtf(sq/162.0f) <= c->minrms ) return;

    deinterleave(symbols);

    if( c->stackdecoder ) {
        not_decoded = jelinek(&metric,&cycles,decdata,symbols,c->nbits,
                              c->stacksize,stack,c->mettab,c->maxcycles);
    } else {
        not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,c->nbits,
                           c->mettab,c->delta,c->maxcycles);
    }

    if( (c->ndepth >= 0) && not_decoded ) {
        for(i=0; i<162; i++) fsymbs[i]=symbols[i]-128.0f;
        int nd=c->ndepth;
        osdwspr_(fsymbs,apmask,&nd,cw,&nhardmin,&dmin);

        unsigned char osdsym[162*2];
        memcpy(osdsym,symbols,sizeof osdsym);
        for(i=0; i<162; i++) osdsym[i]=255*cw[i];
        fano(&metric,&cycles,&maxnp,decdata,osdsym,c->nbits,
             c->mettab,c->delta,c->maxcycles);
        for(i=0; i<11; i++)
            message[i] = decdata[i]>127 ? (signed char)(decdata[i]-256)
                                        : (signed char)decdata[i];
        int n1,n2,n3,nadd,nu,ntype,itype,ihash;
        unpack50(message,&n1,&n2);
        /* An OSD codeword that will not unpack.  Stock wsprd breaks out of the
           whole DT loop here; the caller reproduces that from this flag. */
        if( !unpackcall(n1,callsign) ) { out->abandoned=1; return; }
        callsign[12]=0;
        if( !unpackgrid(n2, grid) ) { out->abandoned=1; return; }
        grid[4]=0;
        ntype = (n2&127) - 64;
        if( (ntype >= 0) && (ntype <= 62) ) {
            nu = ntype%10;
            itype=1;
            if( !(nu == 0 || nu == 3 || nu == 7) ) {
                nadd=nu;
                if( nu > 3 ) nadd=nu-3;
                if( nu > 7 ) nadd=nu-7;
                n3=n2/128+32768*(nadd-1);
                if( !unpackpfx(n3,callsign) ) { out->abandoned=1; return; }
                itype=2;
            }
            ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
            if(strncmp(c->hashtab+ihash*13,callsign,13)==0) {
                if( (itype==1 && strncmp(c->loctab+ihash*5,grid,5)==0) ||
                    (itype==2) ) {
                   not_decoded=0;
                   osd_decode =1;
                }
            }
        }
    }

    if( not_decoded ) return;
    out->decoded=1;
    out->metric=metric; out->cycles=cycles;
    out->nhardmin=nhardmin; out->dmin=dmin; out->osd_decode=osd_decode;
    memcpy(out->decdata,decdata,11);
    memcpy(out->symbols,symbols,162);
}

//***************************************************************************
void usage(void)
{
    printf("Usage: wsprd [options...] infile\n");
    printf("       infile must have suffix .wav or .c2\n");
    printf("\n");
    printf("Options:\n");
    printf("       -a <path> path to writeable data files, default=\".\"\n");
    printf("       -B disable block demodulation - use single-symbol noncoherent demod\n");
    printf("       -c write .c2 file at the end of the first pass\n");
    printf("       -C maximum number of decoder cycles per bit, default 10000\n");
    printf("       -d deeper search. Slower, a few more decodes\n");
    printf("       -e x (x is transceiver dial frequency error in Hz)\n");
    printf("       -f x (x is transceiver dial frequency in MHz)\n");
    printf("       -H do not use (or update) the hash table\n");
    printf("       -J use the stack decoder instead of Fano decoder\n");
    printf("       -m decode wspr-15 .wav file\n");
    printf("       -n x noise blanker: blank the strongest x%% of samples (0=off, default 0);\n");
    printf("            -n a sweeps 0,3,6,10,15,20%% and keeps every decode any of them finds\n");
    printf("       -N n demodulator trials on the last pass (default 4 = unchanged):\n");
    printf("            6 adds noncoherent block lengths 6 and 9,\n");
    printf("            7..14 add coherent demodulation, smoothing 162,81,41,27,15,9,5,3,\n");
    printf("            15..20 add fade-weighted retries of blocks 1,2,3 and\n");
    printf("            coherent 27,9,5, which track a slowly fading signal\n");
    printf("       -X n narrowband interference excision: 1 = on (default), 0 = off,\n");
    printf("            2 = aggressive, a = sweep, or T1,T2,T3,W[,c] for raw thresholds\n");
    printf("       -o n (0<=n<=5), decoding depth for OSD, default is disabled\n");
    printf("       -P n worker threads: 0 = one per logical processor (default), 1 = serial\n");
    printf("       -r remember where stations were decoded for half an hour, and at\n");
    printf("          full depth try those frequencies even without a candidate\n");
    printf("       -G reject a bare prefixed-callsign message from the Fano path when\n");
    printf("          it disagrees with more than 24 of the demodulated symbols, the\n");
    printf("          signature of a wrong codeword; hash-verified results are exempt\n");
    printf("       -A after an unpackable OSD result keep trying the remaining DT\n");
    printf("          offsets, rather than abandoning them as stock wsprd does\n");
    printf("       -S s1,s2 override the sync gates: s1 admits candidates to the fine\n");
    printf("          lag and frequency refinement (default 0.10), s2 admits them to\n");
    printf("          the demodulator (default 0.12, final pass 0.10); negative keeps\n");
    printf("          the default.  -S 0,-1 refines every candidate\n");
    printf("       -Y n decode passes (default 3); extra passes repeat the final pass\n");
    printf("          on what remains after another round of subtraction\n");
    printf("       -q quick mode - doesn't dig deep for weak signals\n");
    printf("       -s single pass mode, no subtraction (same as original wsprd)\n");
    printf("       -v verbose mode (shows dupes)\n");
    printf("       -w wideband mode - decode signals within +/- 150 Hz of center\n");
    printf("       -z x (x is fano metric table bias, default is 0.45)\n");
}

//***************************************************************************
int main(int argc, char *argv[])
{
    char cr[] = "(C) 2018, Steven Franke - K9AN";
    (void)cr;
    extern char *optarg;
    extern int optind;
    int i,j,k;
    unsigned char *symbols, *decdata, *channel_symbols;
    signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
    char *callsign, *grid,  *call_loc_pow;
    char *ptr_to_infile,*ptr_to_infile_suffix;
    char *data_dir=".";
    char wisdom_fname[200],all_fname[200],spots_fname[200];
    char timer_fname[200],hash_fname[200],rv_fname[220];
    char uttime[5],date[7];
    int c,delta,maxpts=65536,verbose=0,quickmode=0,more_candidates=0, stackdecoder=0;
    int usehashtable=1,wspr_type=2, ipass, nblocksize;
    /* Demodulator trials made on the final pass.  These are block lengths for
       noncoherent_sequence_detection, which supports 1,2,3,6,9 -- its arrays are
       sized 512 = 2^9 for exactly that -- but which wsprd has never asked for
       beyond 3.  Negative entries select the coherent demodulator with that many
       symbols of channel smoothing.  Trials run in order and the first success
       wins, so added trials can only ever add decodes. */
    static const int bstab[20]={1,2,3,1,6,9,-162,-81,-41,-27,-15,-9,-5,-3,
                               1,2,3,-27,-9,-5};   /* the last six are weighted */
    static const int bmtab[20]={0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    static const float bwtab[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0.5f,0.5f,0.5f,0.5f,0.5f,0.5f};
    int nbtrials=4;
    int nhardmin;
    int writec2=0,maxdrift;
    int shift1, not_decoded;
    unsigned int nbits=81, stacksize=200000;
    unsigned int npoints, cycles, metric;
    float df=375.0/256.0/2;
    float dt=1.0/375.0, dt_print;
    double dialfreq_cmdline=0.0, dialfreq, freq_print;
    double dialfreq_error=0.0;
    float fmin=-110, fmax=110;
    float f1, sync1, drift1;
    float psavg[512];
    float *idat, *qdat;
    clock_t t0,t00;
    float tfano=0.0,treadwav=0.0,tcandidates=0.0,tsync0=0.0;
    float tsync1=0.0,tsync2=0.0,tosd=0.0,ttotal=0.0;
    float ttrials=0.0, twall=0.0;
    double wt0, wtstart=wsecs();
    
    struct cand { float freq; float snr; int shift; float drift; float sync; int prior; };
    struct cand candidates[200];
#define MAXUNIQUE 50
    
    struct result { char date[7]; char time[5]; float sync; float snr;
        float dt; double freq; char message[23]; float drift;
        unsigned int cycles; int jitter; int blocksize; unsigned int metric;
        int nhardmin; int ipass; int decodetype;};
    struct result decodes[50];
    
    char *hashtab;
    hashtab=calloc(32768*13,sizeof(char));
    char *loctab;
    loctab=calloc(32768*5,sizeof(char));
    char *hashsnap=calloc(32768*13,sizeof(char));
    char *locsnap=calloc(32768*5,sizeof(char));
    int nh;
    symbols=calloc(nbits*2,sizeof(unsigned char));
    decdata=calloc(11,sizeof(unsigned char));
    channel_symbols=calloc(nbits*2,sizeof(unsigned char));
    callsign=calloc(13,sizeof(char));
    grid=calloc(5,sizeof(char));
    call_loc_pow=calloc(23,sizeof(char));
    float allfreqs[100];
    char allcalls[100][13];
    for (i=0; i<100; i++) allfreqs[i]=0.0;
    memset(allcalls,0,sizeof(char)*100*13);
    
    int uniques=0, noprint=0, ndecodes_pass=0;

#define RVMAX 60
#define RVAGE 30
    struct rvprior { char rdate[7]; char rtime[5]; float freq, dtv, drift; int stale; };
    struct rvprior rvp[RVMAX];
    int nrv=0;
    
    // Parameters used for performance-tuning:
    unsigned int maxcycles=10000;            //Decoder timeout limit
    float minsync1=0.10;                     //First sync limit
    float minsync2=0.12;                     //Second sync limit
    int iifac=8;                             //Step size in final DT peakup
    int symfac=50;                           //Soft-symbol normalizing factor
    int subtraction=1;
    int npasses=3;
    int ndepth=-1;                            //Depth for OSD
    
    float minrms=52.0 * (symfac/64.0);      //Final test for plausible decoding
    delta=60;                                //Fano threshold step
    float bias=0.45;                        //Fano metric bias (used for both Fano and stack algorithms)
    
    t00=clock();
    fftwf_complex *fftin, *fftout;
#include "./metric_tables.c"
    
    int mettab[2][256];
    
    idat=calloc(maxpts,sizeof(float));
    qdat=calloc(maxpts,sizeof(float));
    
    while ( (c = getopt(argc, argv, "a:ABcC:de:f:GHn:N:P:JmrS:o:qstwvX:Y:z:")) !=-1 ) {
        switch (c) {
            case 'a':
                data_dir = optarg;
                break;
            case 'B':
                npasses=2;
                break;
            case 'n':
                if(optarg[0]=='a'){ g_nbauto=1; break; }
                g_nbpct=(int)strtol(optarg,NULL,10);
                if(g_nbpct<0) g_nbpct=0;
                if(g_nbpct>50) g_nbpct=50;
                break;
            case 'A':
                g_keepdt=1;
                break;
            case 'G':
                g_guard=1;
                break;
            case 'r':
                g_recall=1;
                break;
            case 'S':
                sscanf(optarg,"%f,%f",&g_minsync1,&g_minsync2);
                break;
            case 'Y':
                npasses=(int)strtol(optarg,NULL,10);
                if(npasses<1) npasses=1;
                if(npasses>8) npasses=8;
                break;
            case 'P':
                g_nthreads=(int)strtol(optarg,NULL,10);
                break;
            case 'X': {
                if(optarg[0]=='a'){ g_fexauto=1; break; }
                float t1,t2,t3; int w=0; char cl=0;
                int nf=sscanf(optarg,"%f,%f,%f,%d,%c",&t1,&t2,&t3,&w,&cl);
                if(nf>=3){
                    g_fex=1; g_fexT1=t1; g_fexT2=t2; g_fexT3=t3;
                    g_fexW=(w>0)?w:2048; g_fexclamp=(cl=='c');
                } else {
                    int lev=(int)strtol(optarg,NULL,10);
                    g_fex=lev;
                    if(lev>=2){ g_fexT1=300.0f; g_fexT2=100.0f; }
                }
                break; }
            case 'N':
                nbtrials=(int)strtol(optarg,NULL,10);
                if(nbtrials<1) nbtrials=1;
                if(nbtrials>20) nbtrials=20;
                break;
            case 'c':
                writec2=1;
                break;
            case 'C':
                maxcycles=(unsigned int) strtoul(optarg,NULL,10);
                break;
            case 'd':
                more_candidates=1;
                break;
            case 'e':
                dialfreq_error = strtod(optarg,NULL);   // units of Hz
                // dialfreq_error = dial reading - actual, correct frequency
                break;
            case 'f':
                dialfreq_cmdline = strtod(optarg,NULL); // units of MHz
                break;
            case 'H':
                usehashtable = 0;
                break;
            case 'J': //Stack (Jelinek) decoder, Fano decoder is the default
                stackdecoder = 1;
                break;
            case 'm':  //15-minute wspr mode
                wspr_type = 15;
                break;
            case 'o':  //use ordered-statistics-decoder
                ndepth=(int) strtol(optarg,NULL,10);
                break;
            case 'q':  //no shift jittering
                quickmode = 1;
                break;
            case 's':  //single pass mode
                subtraction = 0;
                npasses = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'w':
                fmin=-150.0;
                fmax=150.0;
                break;
            case 'z':
                bias=strtod(optarg,NULL); //fano metric bias (default is 0.45)
                break;
            case '?':
                usage();
                return 1;
        }
    }
    
    if( access(data_dir, R_OK | W_OK)) {
      fprintf(stderr, "Error: inaccessible data directory: '%s'\n", data_dir);
      usage();
      return EXIT_FAILURE;
    }

    if( optind+1 > argc) {
        usage();
        return 1;
    } else {
        ptr_to_infile=argv[optind];
    }
    
#ifdef _OPENMP
    if( g_nthreads > 0 ) omp_set_num_threads(g_nthreads);
    int nworkers = omp_get_max_threads();
#else
    int nworkers = 1;
#endif
    vensure();                      /* build the FFTW plan before any threads */
    struct snode **stacks = NULL;
    if( stackdecoder ) {            /* the stack decoder needs one per worker */
        stacks = calloc(nworkers,sizeof(struct snode*));
        for(i=0;i<nworkers;i++) stacks[i]=calloc(stacksize,sizeof(struct snode));
    }

    // setup metric table
    for(i=0; i<256; i++) {
        mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
        mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
    }
    
    FILE *fp_fftwf_wisdom_file, *fall_wspr, *fwsprd, *fhash, *ftimer;
    strcpy(wisdom_fname,".");
    strcpy(all_fname,".");
    strcpy(spots_fname,".");
    strcpy(timer_fname,".");
    strcpy(hash_fname,".");
    if(data_dir != NULL) {
      strncpy(wisdom_fname,data_dir, sizeof wisdom_fname);
      strncpy(all_fname,data_dir, sizeof all_fname);
      strncpy(spots_fname,data_dir, sizeof spots_fname);
      strncpy(timer_fname,data_dir, sizeof timer_fname);
      strncpy(hash_fname,data_dir, sizeof hash_fname);
    }
    strncat(wisdom_fname,"/wspr_wisdom.dat",20);
    strncat(all_fname,"/ALL_WSPR.TXT",20);
    strncat(spots_fname,"/wspr_spots.txt",20);
    strncat(timer_fname,"/wspr_timer.out",20);
    strncat(hash_fname,"/hashtable.txt",20);
    snprintf(rv_fname,sizeof rv_fname,"%s/wspr_priors.txt",data_dir);
    if ((fp_fftwf_wisdom_file = fopen(wisdom_fname, "r"))) {  //Open FFTW wisdom
        fftwf_import_wisdom_from_file(fp_fftwf_wisdom_file);
        fclose(fp_fftwf_wisdom_file);
    }
    
    fall_wspr=fopen(all_fname,"a");
    fwsprd=fopen(spots_fname,"w");
    //  FILE *fdiag;
    //  fdiag=fopen("wsprd_diag","a");
    
    if((ftimer=fopen(timer_fname,"r"))) {
        //Accumulate timing data
        int nr=fscanf(ftimer,"%f %f %f %f %f %f %f %f",
               &treadwav,&tcandidates,&tsync0,&tsync1,&tsync2,&tfano,&tosd,&ttotal);
        fclose(ftimer);
        if(nr == 0) fprintf(stderr, "Empty timer file: '%s'\n", timer_fname);
    }
    ftimer=fopen(timer_fname,"w");
    
    int isc2 = 0;
    if( strstr(ptr_to_infile,".wav") ) {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".wav");
        
        t0 = clock();
        npoints=readwavfile(ptr_to_infile, wspr_type, idat, qdat);
        treadwav += (float)(clock()-t0)/CLOCKS_PER_SEC;
        
        if( npoints == 1 ) {
            return 1;
        }
        dialfreq=dialfreq_cmdline - (dialfreq_error*1.0e-06);
    } else if ( strstr(ptr_to_infile,".c2") !=0 )  {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".c2");
        isc2 = 1;
        npoints=readc2file(ptr_to_infile, idat, qdat, &dialfreq, &wspr_type);
        if( npoints == 1 ) {
            return 1;
        }
        dialfreq -= (dialfreq_error*1.0e-06);
    } else {
        printf("Error: Failed to open %s\n",ptr_to_infile);
        printf("WSPR file must have suffix .wav or .c2\n");
        return 1;
    }
    
    // Parse date and time from given filename
    strncpy(date,ptr_to_infile_suffix-11,6);
    strncpy(uttime,ptr_to_infile_suffix-4,4);
    date[6]='\0';
    uttime[4]='\0';
    
    // Do windowed ffts over 2 symbols, stepped by half symbols
    int nffts=4*floor(npoints/512)-1;
    fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
    fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
    PLAN3 = fftwf_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, PATIENCE);
    
    float ps[512][nffts];
    float w[512];
    for(i=0; i<512; i++) {
        w[i]=sin(0.006147931*i);
    }
    
    if( usehashtable ) {
        char line[80], hcall[13], hgrid[5];
        if( (fhash=fopen(hash_fname,"r+")) ) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
                hgrid[0]='\0';
                sscanf(line,"%d %s %s",&nh,hcall,hgrid);
                strcpy(hashtab+nh*13,hcall);
                if(strlen(hgrid)>0) strcpy(loctab+nh*5,hgrid);
            }
        } else {
            fhash=fopen(hash_fname,"w+");
        }
        fclose(fhash);
    }
    
    if( g_recall && wspr_type==2 ) {
        FILE *frv=fopen(rv_fname,"r");
        if( frv ) {
            char line[160], rdate[16], rtime[16];
            double rdial, rfo, rdt, rdrift;
            int yy,mo,dd,hh,mn;
            long now=0;
            if( sscanf(date,"%2d%2d%2d",&yy,&mo,&dd)==3 &&
                sscanf(uttime,"%2d%2d",&hh,&mn)==2 )
                now=rv_minutes(yy,mo,dd,hh,mn);
            while( fgets(line,sizeof line,frv) && nrv<RVMAX ) {
                if( sscanf(line,"%15s %15s %lf %lf %lf %lf",
                           rdate,rtime,&rdial,&rfo,&rdt,&rdrift) != 6 ) continue;
                if( fabs(rdial-dialfreq) > 1.0e-4 ) continue;   /* other band */
                if( sscanf(rdate,"%2d%2d%2d",&yy,&mo,&dd)!=3 ) continue;
                if( sscanf(rtime,"%2d%2d",&hh,&mn)!=2 ) continue;
                if( strlen(rdate)!=6 || strlen(rtime)!=4 ) continue;
                long tmin=rv_minutes(yy,mo,dd,hh,mn);
                if( now - tmin > RVAGE || now - tmin < 0 ) continue;
                int dupe=0;
                for (i=0; i<nrv; i++)
                    if( fabsf(rvp[i].freq-(float)rfo) < 0.2f ) { dupe=1; break; }
                if( dupe ) continue;
                strcpy(rvp[nrv].rdate,rdate);
                strcpy(rvp[nrv].rtime,rtime);
                rvp[nrv].freq=(float)rfo;
                rvp[nrv].dtv=(float)rdt;
                rvp[nrv].drift=(float)rdrift;
                rvp[nrv].stale=0;
                nrv++;
            }
            fclose(frv);
        }
    }

    // Compute corrected fmin, fmax, accounting for dial frequency error.  This
    // belongs outside the pass loop: applied once per pass, three passes moved
    // the search window by three times the error.
    fmin += dialfreq_error;        // dialfreq_error is in units of Hz
    fmax += dialfreq_error;

    /* Noise-blanker sweep.  A fixed setting is a bad idea: on clean gaussian
       noise, blanking 10% of the samples costs about 2.6 dB.  Trying several and
       keeping every decode any of them finds costs nothing on a quiet band and
       is worth several dB on an impulsive one.  Same shape as the NB=Auto that
       FST4 and FST4W already have. */
    static const int nbtab[6]={0,3,6,10,15,20};
    static const float nbktab[6]={0.0f,8.0f,6.0f,5.0f,4.0f,3.0f};
    int nnb = (g_nbauto && !isc2) ? 6 : 1;
    if( nnb>1 && wspr_type==2 ) g_oobmode=1;
    long lastnh=0;   /* hit count of the previous executed round; 0 = unblanked */
    for (int inb=0; inb<nnb; inb++) {
      if( g_nbauto && !isc2 ) {
        if( g_oobmode ) { g_nbk=nbktab[inb]; g_nbpct=0; }
        else g_nbpct = nbtab[inb];
        if( g_oobmode && inb>0 && g_oobenv!=NULL ) {
            /* Threshold hit sets are nested, so an equal count means an
               identical set, an identical buffer, and identical decodes:
               skip the round rather than repeat it. */
            float th=g_nbk*g_oobsig;
            long nh=0, ii;
            for(ii=0;ii<114L*12000L;ii++) if(g_oobenv[ii]>th) nh++;
            if( nh==lastnh ) continue;
            lastnh=nh;
        }
        npoints = readwavfile(ptr_to_infile, wspr_type, idat, qdat);
        if( npoints == 1 ) break;
      }
    //*************** main loop starts here *****************
    for (ipass=0; ipass<npasses; ipass++) {
        if(ipass==1 && ndecodes_pass == 0 && npasses>2) ipass=2;
        /* An extra pass exists to search what new subtraction has uncovered.
           If the previous full-depth pass decoded nothing, nothing was
           subtracted, the buffers are identical, and the next pass would
           provably repeat the same failures: stop. */
        if(ipass>=3 && ndecodes_pass == 0) break;
        if(ipass < 2) {
            nblocksize=1;
            maxdrift=4;
            minsync2=0.12;
        }
        if(ipass >= 2 ) {
            nblocksize=nbtrials;
            maxdrift=0;    // no drift for smaller frequency estimator variance
            minsync2=0.10;
        }
        if(g_minsync1>=0.0f) minsync1=g_minsync1;
        if(g_minsync2>=0.0f) minsync2=g_minsync2;
        ndecodes_pass=0;   // still needed?
        
        /* one independent 512-point transform per half-symbol step */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (i=0; i<nffts; i++) {
            fftwf_complex lin[512], lout[512];
            int jj,kk;
            for(jj=0; jj<512; jj++ ) {
                kk=i*128+jj;
                lin[jj][0]=idat[kk] * w[jj];
                lin[jj][1]=qdat[kk] * w[jj];
            }
            fftwf_execute_dft(PLAN3,lin,lout);
            for (jj=0; jj<512; jj++ ) {
                kk=jj+256;
                if( kk>511 )
                    kk=kk-512;
                ps[jj][i]=lout[kk][0]*lout[kk][0]+lout[kk][1]*lout[kk][1];
            }
        }
        
        // Compute average spectrum
        for (i=0; i<512; i++) psavg[i]=0.0;
        for (i=0; i<nffts; i++) {
            for (j=0; j<512; j++) {
                psavg[j]=psavg[j]+ps[j][i];
            }
        }
        
        /* The average spectrum has been taken, and every remaining reader of
           ps[][] wants its square root, so take it once here rather than a
           million times inside the coarse search below. */
        for (i=0; i<nffts; i++) {
            for (j=0; j<512; j++) {
                ps[j][i]=sqrtf(ps[j][i]);
            }
        }

        // Smooth with 7-point window and limit spectrum to +/-150 Hz
        int window[7]={1,1,1,1,1,1,1};
        float smspec[411];
        for (i=0; i<411; i++) {
            smspec[i]=0.0;
            for(j=-3; j<=3; j++) {
                k=256-205+i+j;
                smspec[i]=smspec[i]+window[j+3]*psavg[k];
            }
        }
        
        // Sort spectrum values, then pick off noise level as a percentile
        float tmpsort[411];
        for (j=0; j<411; j++) {
            tmpsort[j]=smspec[j];
        }
        qsort(tmpsort, 411, sizeof(float), floatcomp);
        
        // Noise level of spectrum is estimated as 123/411= 30'th percentile
        float noise_level = tmpsort[122];
        
        /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
         * We know from experience that threshold snr is near -7dB in wspr bandwidth,
         * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
         * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */
        
        float min_snr, snr_scaling_factor;
        min_snr = pow(10.0,-8.0/10.0); //this is min snr in wspr bw
        if( wspr_type == 2 ) {
            snr_scaling_factor=26.3;
        } else {
            snr_scaling_factor=35.3;
        }
        for (j=0; j<411; j++) {
            smspec[j]=smspec[j]/noise_level - 1.0;
            if( smspec[j] < min_snr) smspec[j]=0.1*min_snr;
            continue;
        }
        
        // Find all local maxima in smoothed spectrum.
        for (i=0; i<200; i++) {
            candidates[i].freq=0.0;
            candidates[i].snr=0.0;
            candidates[i].drift=0.0;
            candidates[i].shift=0;
            candidates[i].sync=0.0;
            candidates[i].prior=0;
        }
        
        int npk=0;
        unsigned char candidate;
        for(j=1; j<410; j++) {
            candidate = (smspec[j]>smspec[j-1]) &&
            (smspec[j]>smspec[j+1]) &&
            (npk<200);
            if ( candidate ) {
                candidates[npk].freq = (j-205)*df;
                candidates[npk].snr  = 10*log10(smspec[j])-snr_scaling_factor;
                npk++;
            }
        }
        if( more_candidates ) {
            for(j=0; j<411; j=j+3) {
                candidate = (smspec[j]>min_snr) && (npk<200);
                if ( candidate ) {
                    candidates[npk].freq = (j-205)*df;
                    candidates[npk].snr  = 10*log10(smspec[j])-snr_scaling_factor;
                    npk++;
                }
            }
        }
        
        // Don't waste time on signals outside of the range [fmin,fmax].
        i=0;
        for( j=0; j<npk; j++) {
            if( candidates[j].freq >= fmin && candidates[j].freq <= fmax ) {
                candidates[i]=candidates[j];
                i++;
            }
        }
        npk=i;
        
        // bubble sort on snr
        int pass;
        struct cand tmp;
        for (pass = 1; pass <= npk - 1; pass++) {
            for (k = 0; k < npk - pass ; k++) {
                if (candidates[k].snr < candidates[k+1].snr) {
                    tmp = candidates[k];
                    candidates[k]=candidates[k+1];
                    candidates[k+1] = tmp;
                }
            }
        }
        
        t0=clock();
        
        /* Make coarse estimates of shift (DT), freq, and drift
         
         * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
         to nominal start time, which is 2 seconds into the file
         
         * Calculates shift relative to the beginning of the file
         
         * Negative shifts mean that signal started before start of file
         
         * The program prints DT = shift-2 s
         
         * Shifts that cause sync vector to fall off of either end of the data
         vector are accommodated by "partial decoding", such that missing
         symbols produce a soft-decision symbol value of 128
         
         * The frequency drift model is linear, deviation of +/- drift/2 over the
         span of 162 symbols, with deviation equal to 0 at the center of the
         signal vector.
         */
        
        /* Each candidate is refined into its own slot and nothing else is
           written, so this whole search is independent per candidate. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for(j=0; j<npk; j++) {                              //For each candidate...
            int idrift,ifr,if0,ifd,k0,kindex,k;
            float smax,ss,pow,p0,p1,p2,p3;
            smax=-1e30;
            if0=candidates[j].freq/df+256;
            for (ifr=if0-2; ifr<=if0+2; ifr++) {                      //Freq search
                for( k0=-10; k0<22; k0++) {                             //Time search
                    for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {  //Drift search
                        ss=0.0;
                        pow=0.0;
                        for (k=0; k<162; k++) {                             //Sum over symbols
                            ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/(2.0*df);
                            kindex=k0+2*k;
                            if( kindex >= 0 && kindex < nffts ) {
                                p0=ps[ifd-3][kindex];   /* already square-rooted */
                                p1=ps[ifd-1][kindex];
                                p2=ps[ifd+1][kindex];
                                p3=ps[ifd+3][kindex];
                                
                                ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                                pow=pow+p0+p1+p2+p3;
                            }
                        }
                        sync1=ss/pow;
                        if( sync1 > smax ) {                  //Save coarse parameters
                            smax=sync1;
                            candidates[j].shift=128*(k0+1);
                            candidates[j].drift=idrift;
                            candidates[j].freq=(ifr-256)*df;
                            candidates[j].sync=sync1;
                        }
                    }
                }
            }
        }
        tcandidates += (float)(clock()-t0)/CLOCKS_PER_SEC;
        
        /*
         Refine the estimates of freq, shift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].
         
         Function sync_and_demodulate has three modes of operation
         mode is the last argument:
         
         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.
         
         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
         */
        wt0 = wsecs();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
        for (j=0; j<npk; j++) {
            /* everything below is per-candidate scratch, so declare it here
               rather than sharing the enclosing scope's variables */
            unsigned char lsym[162*2];
            float cf1, cdrift, csync, cfstep;
            int cshift, cifmin, cifmax, clagmin, clagmax, clagstep;

            cf1=candidates[j].freq;
            cdrift=candidates[j].drift;
            cshift=candidates[j].shift;
            csync=candidates[j].sync;
            
            // coarse-grid lag and freq search, then if sync>minsync1 continue
            cfstep=0.0; cifmin=0; cifmax=0;
            clagmin=cshift-128;
            clagmax=cshift+128;
            clagstep=64;
            sync_and_demodulate(idat, qdat, npoints, lsym, &cf1, cifmin, cifmax, cfstep, &cshift,
                                clagmin, clagmax, clagstep, &cdrift, symfac, &csync, 0);
            
            cfstep=0.25; cifmin=-2; cifmax=2;
            sync_and_demodulate(idat, qdat, npoints, lsym, &cf1, cifmin, cifmax, cfstep, &cshift,
                                clagmin, clagmax, clagstep, &cdrift, symfac, &csync, 1);
            
            if(ipass < 2) {
                // refine drift estimate
                cfstep=0.0; cifmin=0; cifmax=0;
                float driftp,driftm,syncp,syncm;
                driftp=cdrift+0.5;
                sync_and_demodulate(idat, qdat, npoints, lsym, &cf1, cifmin, cifmax, cfstep, &cshift,
                                    clagmin, clagmax, clagstep, &driftp, symfac, &syncp, 1);
                
                driftm=cdrift-0.5;
                sync_and_demodulate(idat, qdat, npoints, lsym, &cf1, cifmin, cifmax, cfstep, &cshift,
                                    clagmin, clagmax, clagstep, &driftm, symfac, &syncm, 1);
                
                if(syncp>csync) {
                    cdrift=driftp;
                    csync=syncp;
                } else if (syncm>csync) {
                    cdrift=driftm;
                    csync=syncm;
                }
            }
            
            // fine-grid lag and freq search
            if( csync > minsync1 ) {
                
                clagmin=cshift-32; clagmax=cshift+32; clagstep=16;
                sync_and_demodulate(idat, qdat, npoints, lsym, &cf1, cifmin, cifmax, cfstep, &cshift,
                                    clagmin, clagmax, clagstep, &cdrift, symfac, &csync, 0);
                
                // fine search over frequency
                cfstep=0.05; cifmin=-2; cifmax=2;
                sync_and_demodulate(idat, qdat, npoints, lsym, &cf1, cifmin, cifmax, cfstep, &cshift,
                                    clagmin, clagmax, clagstep, &cdrift, symfac, &csync, 1);
                
                candidates[j].freq=cf1;
                candidates[j].shift=cshift;
                candidates[j].drift=cdrift;
                candidates[j].sync=csync;
            }
        }
        tsync1 += (float)(wsecs()-wt0);
        
        int nwat=0; 
        int idupe;
        for ( j=0; j<npk; j++) {
            idupe=0;
            for (k=0;k<nwat;k++) {
                if( fabsf(candidates[j].freq - candidates[k].freq) < 0.05  &&
                   abs(candidates[j].shift - candidates[k].shift) < 16 ) {
                    idupe=1;
                    break;
                }
            }
            if( idupe == 1 ) {
                if(candidates[j].sync > candidates[k].sync) candidates[k]=candidates[j];
            } else if ( candidates[j].sync > minsync2 ) {
                candidates[nwat]=candidates[j];
                nwat++;
            }
        }

        /* Return-visitor recall: a station decoded in the last half hour is
           usually still there, so at full depth its last known frequency and
           DT get a decode attempt even when candidate detection offers
           nothing.  Acceptance is unchanged: the attempt still has to satisfy
           Fano or the hash-gated OSD, so this is only extra search, aimed
           where stations are known to live. */
        if( g_recall && wspr_type==2 && ipass==2 ) {
            for (j=0; j<nrv && j<20 && nwat<200; j++) {
                int nearby=0;
                for (k=0; k<nwat && !nearby; k++)
                    if( fabsf(candidates[k].freq - rvp[j].freq) < 0.5f ) nearby=1;
                for (k=0; k<uniques && !nearby; k++)
                    if( fabsf(allfreqs[k] - rvp[j].freq) < 0.5f ) nearby=1;
                if( nearby ) continue;
                candidates[nwat].freq=rvp[j].freq;
                candidates[nwat].snr=-33.0;
                candidates[nwat].shift=(int)(375.0f*(1.0f+rvp[j].dtv));
                candidates[nwat].drift=rvp[j].drift;
                candidates[nwat].sync=1.0f;
                candidates[nwat].prior=1;
                nwat++;
            }
        }
        
        struct dec_ctx ctx;
        ctx.idat=idat; ctx.qdat=qdat; ctx.npoints=npoints;
        ctx.mettab=mettab; ctx.delta=delta; ctx.ndepth=ndepth;
        ctx.symfac=symfac; ctx.stackdecoder=stackdecoder;
        ctx.maxcycles=maxcycles; ctx.nbits=nbits; ctx.stacksize=stacksize;
        ctx.minrms=minrms; ctx.hashtab=hashtab; ctx.loctab=loctab;

        int ii;
        int blocksize, bitmetric;
        int osd_decode;
        for (j=0; j<nwat; j++) {
            memset(symbols,0,sizeof(char)*nbits*2);
            memset(callsign,0,sizeof(char)*13);
            memset(grid,0,sizeof(char)*5);
            memset(call_loc_pow,0,sizeof(char)*23);
            f1=candidates[j].freq;
            shift1=candidates[j].shift;
            drift1=candidates[j].drift;
            not_decoded=1;
            osd_decode=0;
            
            /* Flatten the (demodulator setting, DT jitter) grid into one list so
               it can be run in parallel.  The attempts are independent; we keep
               the lowest-numbered one that succeeded, which is exactly where the
               sequential search would have stopped, so the answer does not
               depend on how many workers ran. */
            {
            /* Stock wsprd walks the demodulator settings in order and, within
               each one, walks the DT offsets in order until an attempt either
               decodes or produces an OSD codeword that will not unpack -- the
               latter abandons the rest of that setting's offsets.  So the
               settings stay sequential, there being only a handful of them,
               and the offsets, of which there are up to 129, run in parallel
               and keep the lowest-numbered attempt that reached either
               outcome.  That is the attempt stock would have stopped at. */
            int nidt = quickmode ? 1 : (candidates[j].prior ? 5 : (128/iifac + 1));
            int found = 0;
            struct trial_result best;
            best.decoded = 0;
            wt0 = wsecs();
            for (int tb=0; tb<nblocksize && !found; tb++) {
                int firstev = nidt;
                struct trial_result evr;
                evr.decoded = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
                for (int td=0; td<nidt; td++) {
                    int seen;
#ifdef _OPENMP
#pragma omp atomic read
#endif
                    seen = firstev;
                    if( td > seen ) continue;   /* an earlier offset settled it */
                    int jj = (td+1)/2;
                    if( td%2 == 1 ) jj = -jj;
                    jj = iifac*jj;
                    int tid = 0;
#ifdef _OPENMP
                    tid = omp_get_thread_num();
#endif
                    struct trial_result r;
                    run_trial(&ctx, candidates[j].freq, candidates[j].shift,
                              candidates[j].drift, bstab[tb], bmtab[tb], jj,
                              bwtab[tb], stacks ? stacks[tid] : NULL, &r);
                    if( r.decoded || (r.abandoned && !g_keepdt) ) {
#ifdef _OPENMP
#pragma omp critical (wsprd_trialwin)
#endif
                        { if( td < firstev ) { firstev = td; evr = r; } }
                    }
                }
                if( firstev < nidt && evr.decoded ) {
                    if( candidates[j].prior && !evr.osd_decode ) {
                        /* an injected prior may only decode through the
                           hash-gated deep search: the station it claims to be
                           is in the table, and a bare Fano accept here is the
                           wrong-codeword channel */
                        ;
                    } else {
                        found = 1; best = evr;
                    }
                }
            }
            ttrials += (float)(wsecs()-wt0);
            if( found ) {
                not_decoded    = 0;
                blocksize      = best.blocksize;
                bitmetric      = best.bitmetric;
                ii             = best.jitter;
                metric         = best.metric;
                cycles         = best.cycles;
                nhardmin       = best.nhardmin;
                osd_decode     = best.osd_decode;
                memcpy(decdata,best.decdata,11);
                memcpy(symbols,best.symbols,162);
            }
            }
            
            if( !not_decoded ) {
                ndecodes_pass++;
                
                for(i=0; i<11; i++) {
                    
                    if( decdata[i]>127 ) {
                        message[i]=decdata[i]-256;
                    } else {
                        message[i]=decdata[i];
                    }
                    
                }
                
                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
                if( g_guard ) {
                    memcpy(hashsnap,hashtab,32768*13);
                    memcpy(locsnap,loctab,32768*5);
                }
                noprint=unpk_(message,hashtab,loctab,call_loc_pow,callsign);
                if( subtraction && !noprint ) {
                    if( get_wspr_channel_symbols(call_loc_pow, hashtab, loctab, channel_symbols) ) {
                        if(!osd_decode) nhardmin=count_hard_errors(symbols,channel_symbols);
                        if( g_guard && !osd_decode && nhardmin>24 &&
                            strchr(call_loc_pow,'/')!=NULL && call_loc_pow[0]!='<' ) {
                            /* A bare type 2 from Fano that disagrees with a
                               third of the demodulated symbols: the signature
                               of an accepted wrong codeword.  Do not subtract
                               a signal that is not there, do not print it, and
                               put the hash table back the way it was, so the
                               fabricated call cannot vouch for itself through
                               the hash-gated deep search on a later pass. */
                            memcpy(hashtab,hashsnap,32768*13);
                            memcpy(loctab,locsnap,32768*5);
                            continue;
                        }
                        subtract_signal2(idat, qdat, npoints, f1, shift1, drift1, channel_symbols);
                    } else {
                        break;
                    }
                }
                
                // Remove dupes (same callsign and freq within 4 Hz)
                int dupe=0;
                for (i=0; i<uniques; i++) {
                    if(!strcmp(callsign,allcalls[i]) &&
                       (fabs(f1-allfreqs[i]) <4.0)) dupe=1;
                }
                if( (verbose || !dupe) && !noprint && uniques < MAXUNIQUE) {
                    strcpy(allcalls[uniques],callsign);
                    allfreqs[uniques]=f1;
                    uniques++;
                    
                    // Add an extra space at the end of each line so that wspr-x doesn't
                    // truncate the power (TNX to DL8FCL!)
                    
                    if( wspr_type == 15 ) {
                        freq_print=dialfreq+(1500+112.5+f1/8.0)/1e6;
                        dt_print=shift1*8*dt-1.0;
                    } else {
                        freq_print=dialfreq+(1500+f1)/1e6;
                        dt_print=shift1*dt-1.0;
                    }
                    
                    strcpy(decodes[uniques-1].date,date);
                    strcpy(decodes[uniques-1].time,uttime);
                    decodes[uniques-1].sync=candidates[j].sync;
                    decodes[uniques-1].snr=candidates[j].snr;
                    decodes[uniques-1].dt=dt_print;
                    decodes[uniques-1].freq=freq_print;
                    strcpy(decodes[uniques-1].message,call_loc_pow);
                    decodes[uniques-1].drift=drift1;
                    decodes[uniques-1].cycles=cycles;
                    decodes[uniques-1].jitter=ii;
                    decodes[uniques-1].blocksize=blocksize<0?blocksize:blocksize+3*bitmetric;
                    decodes[uniques-1].metric=metric;
                    decodes[uniques-1].nhardmin=nhardmin;
                    decodes[uniques-1].ipass=ipass;
                    decodes[uniques-1].decodetype=osd_decode;
                }
            }
        }
        
        if( ipass == 0 && writec2 && inb == 0 ) {
            char c2filename[15];
            double carrierfreq=dialfreq;
            int wsprtype=2;
            strcpy(c2filename,"000000_0001.c2");
            printf("Writing %s\n",c2filename);
            writec2file(c2filename, wsprtype, carrierfreq, idat, qdat);
        }
    }
    }
    
    // sort the result in order of increasing frequency
    struct result temp;
    for (j = 1; j <= uniques - 1; j++) {
        for (k = 0; k < uniques - j ; k++) {
            if (decodes[k].freq > decodes[k+1].freq) {
                temp = decodes[k];
                decodes[k]=decodes[k+1];;
                decodes[k+1] = temp;
            }
        }
    }
    
    for (i=0; i<uniques; i++) {
        printf("%4s %3.0f %4.1f %10.6f %2d  %-s \n",
               decodes[i].time, decodes[i].snr,decodes[i].dt, decodes[i].freq,
               (int)decodes[i].drift, decodes[i].message);
        fprintf(fall_wspr,
                "%6s %4s %3.0f %5.2f %11.7f  %-22s %2d %5.2f %2d %2d %4d %2d %3d %5u %5d\n",
                decodes[i].date, decodes[i].time, decodes[i].snr,
                decodes[i].dt, decodes[i].freq, decodes[i].message,
                (int)decodes[i].drift, decodes[i].sync,
                decodes[i].ipass+1,decodes[i].blocksize,decodes[i].jitter,
                decodes[i].decodetype,decodes[i].nhardmin,decodes[i].cycles/81,
                decodes[i].metric);
        fprintf(fwsprd,
                "%6s %4s %3d %3.0f %4.1f %10.6f  %-22s %2d %5u %4d\n",
                decodes[i].date, decodes[i].time, (int)(10*decodes[i].sync),
                decodes[i].snr, decodes[i].dt, decodes[i].freq,
                decodes[i].message, (int)decodes[i].drift, decodes[i].cycles/81,
                decodes[i].jitter);
        
    }
    printf("<DecodeFinished>\n");
    
    fftwf_free(fftin);
    fftwf_free(fftout);
    
    if ((fp_fftwf_wisdom_file = fopen(wisdom_fname, "w"))) {
        fftwf_export_wisdom_to_file(fp_fftwf_wisdom_file);
        fclose(fp_fftwf_wisdom_file);
    }
    
    ttotal += (float)(clock()-t00)/CLOCKS_PER_SEC;
    twall  += (float)(wsecs()-wtstart);
    
    fprintf(ftimer,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n\n",
            treadwav,tcandidates,tsync0,tsync1,tsync2,tfano,tosd,ttotal);
    
    fprintf(ftimer,"Code segment        Seconds   Frac\n");
    fprintf(ftimer,"-----------------------------------\n");
    fprintf(ftimer,"readwavfile        %7.2f %7.2f\n",treadwav,treadwav/ttotal);
    fprintf(ftimer,"Coarse DT f0 f1    %7.2f %7.2f\n",tcandidates,
            tcandidates/ttotal);
    fprintf(ftimer,"sync_and_demod(0)  %7.2f %7.2f\n",tsync0,tsync0/ttotal);
    fprintf(ftimer,"sync_and_demod(1)  %7.2f %7.2f\n",tsync1,tsync1/ttotal);
    fprintf(ftimer,"sync_and_demod(2)  %7.2f %7.2f\n",tsync2,tsync2/ttotal);
    fprintf(ftimer,"Stack/Fano decoder %7.2f %7.2f\n",tfano,tfano/ttotal);
    fprintf(ftimer,"OSD        decoder %7.2f %7.2f\n",tosd,tosd/ttotal);
    fprintf(ftimer,"trial grid (wall)  %7.2f\n",ttrials);
    fprintf(ftimer,"sync refine (wall) %7.2f\n",tsync1);
    fprintf(ftimer,"WALL CLOCK TOTAL   %7.2f  (%d workers)\n",twall,nworkers);
    fprintf(ftimer,"-----------------------------------\n");
    fprintf(ftimer,"Total              %7.2f %7.2f\n",ttotal,1.0);
    
    fclose(fall_wspr);
    fclose(fwsprd);
    //  fclose(fdiag);
    fclose(ftimer);
    fftwf_destroy_plan(PLAN1);
    fftwf_destroy_plan(PLAN2);
    fftwf_destroy_plan(PLAN3);
    
    if( g_recall && wspr_type==2 ) {
        for (i=0; i<uniques; i++) {
            double fo=(decodes[i].freq - dialfreq)*1.0e6 - 1500.0;
            for (j=0; j<nrv; j++)
                if( fabsf(rvp[j].freq-(float)fo) < 0.2f ) rvp[j].stale=1;
        }
        FILE *frv=fopen(rv_fname,"w");
        if( frv ) {
            int nout=0;
            for (i=0; i<uniques && nout<RVMAX; i++) {
                double fo=(decodes[i].freq - dialfreq)*1.0e6 - 1500.0;
                fprintf(frv,"%6s %4s %11.6f %8.2f %6.2f %5.1f\n",
                        date, uttime, dialfreq, fo, decodes[i].dt, decodes[i].drift);
                nout++;
            }
            for (j=0; j<nrv && nout<RVMAX; j++) {
                if( rvp[j].stale ) continue;
                fprintf(frv,"%6s %4s %11.6f %8.2f %6.2f %5.1f\n",
                        rvp[j].rdate, rvp[j].rtime, dialfreq, rvp[j].freq,
                        rvp[j].dtv, rvp[j].drift);
                nout++;
            }
            fclose(frv);
        }
    }

    if( usehashtable ) {
        fhash=fopen(hash_fname,"w");
        for (i=0; i<32768; i++) {
            if( strncmp(hashtab+i*13,"\0",1) != 0 ) {
                fprintf(fhash,"%5d %s %s\n",i,hashtab+i*13,loctab+i*5);
            }
        }
        fclose(fhash);
    }
    
    free(g_rawbuf);
    free(g_oobenv);
    free(hashtab);
    free(loctab);
    free(symbols);
    free(decdata);
    free(channel_symbols);
    free(callsign);
    free(call_loc_pow);
    free(idat);
    free(qdat);
    if(stacks){ for(i=0;i<nworkers;i++) free(stacks[i]); free(stacks); }
    
    return 0;
}
