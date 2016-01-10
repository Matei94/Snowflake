/*** INCLUDES ************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "helpers.h"

/*************************************************************************************************/



/*** DEFINES *************************************************************************************/

#define NRMAX 1002
#define NCMAX 1002

#define KAPPAMAX 64

/* Parameters */
#define RHO_MIN    0.5
#define RHO_MAX    0.9
#define RINIT      0
#define RHORINIT   1
#define BETA_MIN   1.05
#define BETA_MAX   3.0
#define ALPHA_MIN  0.0
#define ALPHA_MAX  0.3
#define THETA_MIN  0.002
#define THETA_MAX  0.009595
#define KAPPA_MIN  0.01
#define KAPPA_MAX  0.05
#define MU_MIN     0.0
#define MU_MAX     0.01
#define GAM        0.0000515
#define SIGMA_MIN -0.5
#define SIGMA_MAX  0.0
#define NR         250
#define SP         2

#define GRAPHICSFILE "snowflake.ppm"

/*************************************************************************************************/



/*** VARIABLES ***********************************************************************************/

double adif [NRMAX][NCMAX]; /* diffusion field              */
int    apic [NRMAX][NCMAX]; /* indicator of snowflake sites */
double afr  [NRMAX][NCMAX]; /* boundary mass                */
double alm  [NRMAX][NCMAX]; /* crystal mass                 */
int    ash  [NRMAX][NCMAX]; /* rings pallette               */

/* Parameters */
double rho;
int rinit;
double rhorinit;
double beta;
double alpha;
double theta;
double kappa;
double mu;
double gam;
double sigma;
int nr, nc;
int sp;

int centeri, centerj;
int parupdate;
int parash;
int rold, rnew;
int frchange;
int stop;

/* Graphics */
Display *td;
Window tw;
GC tgc;
XEvent te;
KeySym tk;
XSizeHints th;
int ts;
unsigned long tf, tb;
char tbuf[8];
int kc;
int noac, pq;

Colormap cmap;
XColor clr[KAPPAMAX];
XColor clron[128];
XColor clroff[128];
XColor othp[20];

extern int red[125], green[125], blue[125];

/*************************************************************************************************/



/*** FUNCTIONS ***********************************************************************************/

/**
 * Create hexagonal boundary
 */
void createbdry() {
  int i,j;

  for (j=2; j<nc; j++) {
    ash[j-1][j]=ash[j][j-1]; ash[j-2][j]=ash[j][j-2];
    adif[j-1][j]=adif[j][j-1]; adif[j-2][j]=adif[j][j-2];
    afr[j-1][j]=afr[j][j-1]; afr[j-2][j]=afr[j][j-2];
    apic[j-1][j]=apic[j][j-1]; apic[j-2][j]=apic[j][j-2];
    alm[j-1][j]=alm[j][j-1]; alm[j-2][j]=alm[j][j-2];
  }

  for (i=2; i<nr; i++) {
    ash[i][0]=ash[i-1][2];
    adif[i][0]=adif[i-1][2];
    afr[i][0]=afr[i-1][2];
    apic[i][0]=apic[i-1][2];
    alm[i][0]=alm[i-1][2];
  }

  ash[0][2]=ash[2][0]; ash[0][1]=ash[2][0]; ash[1][0]=ash[2][0];
  adif[0][2]=adif[2][0]; adif[0][1]=adif[2][0]; adif[1][0]=adif[2][0];
  afr[0][2]=afr[2][0]; afr[0][1]=afr[2][0]; afr[1][0]=afr[2][0];
  apic[0][2]=apic[2][0]; apic[0][1]=apic[2][0]; apic[1][0]=apic[2][0];
  alm[0][2]=alm[2][0]; alm[0][1]=alm[2][0]; alm[1][0]=alm[2][0];

  for (i=1; i<=nr-2; i++){
    j=nr-i;
    ash[i][j]=ash[i][j-1];
    adif[i][j]=adif[i][j-1];
    afr[i][j]=afr[i][j-1];
    apic[i][j]=apic[i][j-1];
    alm[i][j]=alm[i][j-1];
  }

  ash[nr-1][1]=ash[nr-2][1];
  adif[nr-1][1]=adif[nr-2][1];
  afr[nr-1][1]=afr[nr-2][1];
  apic[nr-1][1]=apic[nr-2][1];
  alm[nr-1][1]=alm[nr-2][1];

  ash[nr-2][0]=ash[nr-3][2];
  adif[nr-2][0]=adif[nr-3][2];
  afr[nr-2][0]=afr[nr-3][2];
  apic[nr-2][0]=apic[nr-3][2];
  alm[nr-2][0]=alm[nr-3][2];

  ash[nr-1][0]=ash[nr-3][2];
  adif[nr-1][0]=adif[nr-3][2];
  afr[nr-1][0]=afr[nr-3][2];
  apic[nr-1][0]=apic[nr-3][2];
  alm[nr-1][0]=alm[nr-3][2];
}


/**
 * Big image
 */
void buildbig() {
  int i, j;

  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) {
      if ((i>=1)&&(j>i)&&(i+j<=nr)) {
        ash[i][j]  = ash[j][i];
        adif[i][j] = adif[j][i];
        apic[i][j] = apic[j][i];
        afr[i][j]  = afr[j][i];
        alm[i][j]  = alm[j][i];
      }
    }
  }
}


/**
 * Check mass
 */
void checkmass() {
  int i, j;
  double totalmass;
  totalmass = 0.0;

  buildbig();

  for (i=2; i<=nr-2; i++) {
    for (j=1; i+j <= nr-1; j++) {
      totalmass += adif[i][j]+afr[i][j]+alm[i][j];
    }
  }

  totalmass+=adif[1][1]+afr[1][1]+alm[i][j];
}


/**
 * Algorithm initialization
 */
void initialize() {
  int i, j, k;
  double x;

  pq = 0;
  stop = false;
  parupdate = 0;

  centeri = 1;
  centerj = 1;

  rold = 0;
  rnew = 0;

  for (i=1; i<nr; i++) {
    for (j=1; ((j<=i)&&(i+j<=nr-1)); j++) {
      x = drand48();
      if ((norminf(i-centeri,j-centerj) <= rinit)   &&
          (seminorm(i-centeri, j-centerj) <= rinit) &&
          (x <= rhorinit)) {
        adif[i][j] = 0.0;
        apic[i][j] = 1;
        afr[i][j]  = 0;
        ash[i][j]  = 0;
        alm[i][j]  = 1.0;
        k = norminf(i-centeri, j-centerj);
        if (k>rnew) {
          rnew=k;
        }
      } else {
        adif[i][j] = rho;
        apic[i][j] = 0;
        afr[i][j]  = 0.0;
        ash[i][j]  = 0;
        alm[i][j]  = 0.0;
      }
    }
  }

  rold=rnew;
  parash=1;
  createbdry();
  buildbig();
}


/**
 * Dynamics
 */
void dynamicsdif() {
  double b[NRMAX][NCMAX];
  int i,j;
  int id,iu,jl,jr;
  int count;
  double masscorrection;
  int nrhalf;

  for (i=1;i<nr;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
      b[i][j]=0.0;
    }
  }

  nrhalf=nr/2;
  if (nr%2==0) {
    masscorrection=(1.0/7.0)*(adif[nr-2][2]+adif[nr-3][3]-2.0*adif[nrhalf][nr-nrhalf]);
  } else {
    masscorrection=(1.0/7.0)*(adif[nr-2][2]+adif[nr-3][3]-adif[nrhalf][nr-nrhalf]-adif[nrhalf+1][nr-nrhalf-1]);
  }

  for (i=1;i<nr;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
      if (apic[i][j]==0) {
        id=(i+1);iu=(i-1);
        jr=(j+1); jl=(j-1);
        count=0;
        if (apic[id][j]==0) count++;
        if (apic[iu][j]==0) count++;
        if (apic[i][jl]==0) count++;
        if (apic[i][jr]==0) count++;
        if (apic[iu][jr]==0) count++;
        if (apic[id][jl]==0) count++;

        if (count==0) {
          b[i][j]=adif[i][j];
        } else {
          b[i][j] = (1.0-(double)count/7.0) * adif[i][j] +
                    (adif[id][j]*(1.0-apic[id][j])+adif[iu][j]*(1.0-apic[iu][j]) +
                     adif[i][jl]*(1.0-apic[i][jl])+adif[i][jr]*(1.0-apic[i][jr]) +
                     adif[iu][jr]*(1.0-apic[iu][jr])+adif[id][jl] * (1.0-apic[id][jl])) / 7.0;
        }
      }
    }
  }

  for (i=1;i<nr;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
      if (apic[i][j]==0) adif[i][j]=b[i][j];
    }
  }

  adif[nr-2][1]-=masscorrection;
  createbdry();
}


/**
 * Dynamics
 */
void dynamicspop() {
  double x;
  int i,j;

  for (i=1;i<nr;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {

      x = drand48();
      if (x<0.5) {
        adif[i][j]=adif[i][j]*(1+sigma);
      } else {
        adif[i][j]=adif[i][j]*(1-sigma);
      }
    }
  }

  createbdry();
}


/**
 * Dynamics
 */
void dynamicsunfre() {
  double y,afrij;
  int i,j;

  int iup;

  iup = centeri+rnew+1;
  frchange = false;

  for (i=1;i<=iup;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
       if (apic[i][j]==0) {

        afrij=afr[i][j];
        y=afrij*mu;
        afr[i][j]=afr[i][j]-y;
        adif[i][j]=adif[i][j]+y;

        afrij=alm[i][j];
        if (afrij>0.0){
           y=afrij*gam;
           alm[i][j]=alm[i][j]-y;
           adif[i][j]=adif[i][j]+y;
        }
      }
    }
  }

  createbdry();
}


/**
 * Dynamics
 */
void dynamicsfre() {
  int bpic[NRMAX][NCMAX];
  int i,j,k;
  int id,iu,jl,jr;
  int count;
  double difmass;
  int iup;

  iup = centeri+rnew+1;
  frchange = false;

  for (i=1;i<=iup;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
      bpic[i][j]=apic[i][j];
    }
  }

  for (i=1;i<=iup;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
      if (apic[i][j]==0) {
        id=i+1;iu=i-1;
        jr=j+1; jl=j-1;
        count=0;
        if (apic[id][j]==1) count++;
        if (apic[iu][j]==1) count++;
        if (apic[i][jl]==1) count++;
        if (apic[i][jr]==1) count++;
        if (apic[iu][jr]==1) count++;
        if (apic[id][jl]==1) count++;

        if (count>=1) {
          difmass=adif[i][j]+adif[id][j]*(1-apic[id][j])+adif[iu][j]*(1-apic[iu][j])
          +adif[i][jl]*(1-apic[i][jl])+adif[i][jr]*(1-apic[i][jr])
          +adif[iu][jr]*(1-apic[iu][jr])+adif[id][jl]*(1-apic[id][jl]);

          if (count<=2) {
            if (afr[i][j]>=beta) {
              bpic[i][j]=1;
            }
          }

          if (count>=3) {
            if ((afr[i][j]>=1.0) || ((difmass<=theta)&&(afr[i][j]>=alpha)) ) {
              bpic[i][j]=1;
            }
          }

          if (count>=4)  {
            bpic[i][j]=1;
          }
        }
      }
    }
  }

  for (i=1;i<=iup;i++) {
    for (j=1;((j<=i)&&(i+j<=nr-1));j++) {
      if (apic[i][j]!=bpic[i][j]) {
        apic[i][j]=bpic[i][j];
        alm[i][j]+=afr[i][j];
        afr[i][j]=0.0;
        k=norminf(i-centeri, j-centerj);
        if (k>rnew) {
          rnew=k;
        }
        if (rnew>2*nr/3) {
          stop=true;
        }
        ash[i][j]=parash;
        frchange=true;
      }
    }
  }

  parupdate = 1-parupdate;
  if (rnew-rold==1) {
    parash=parash+1;
    rold=rnew;
  }

  createbdry();
}

void dynamicsfre1() {
  int i,j;
  int id,iu,jl,jr;
  int count;
  double offset;

  int iup;

  iup = centeri+rnew+1;
  frchange = false;

  for (i=1; i<=iup; i++) {
    for (j=1; ((j<=i)&&(i+j<=nr-1)); j++) {

         if (apic[i][j]==0){

            id=i+1;iu=i-1;
            jr=j+1; jl=j-1;
            count=0;
            if (apic[id][j]==1) count++;
            if (apic[iu][j]==1) count++;
            if (apic[i][jl]==1) count++;
            if (apic[i][jr]==1) count++;
            if (apic[iu][jr]==1) count++;
            if (apic[id][jl]==1) count++;

            if (count>=1){
               offset=(1.0-kappa)*adif[i][j];
               afr[i][j]=afr[i][j]+offset;
               offset=adif[i][j]-offset; adif[i][j]=0;  alm[i][j]+=offset;
            }
         }
      }
   }
   createbdry();
}


/**
 * All dynamics
 */
void dynamics() {
  dynamicsdif();
  dynamicsfre1();
  dynamicsfre();
  dynamicsunfre();

  if (sigma > 0.0) {
    dynamicspop();
  }
}


/**
 * Draw whole
 */
void picturebig() {
  int i, j, k;
  double y;

  buildbig();
  for (i=1; i < nr; i++) {
    for (j=1; j < nc; j++) {
      if (apic[i][j] == 0) {
        k = floor(63.0 * (adif[i][j]/(rho)));
        XSetForeground(td, tgc, clroff[k].pixel);
        XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j*sp+30, i*sp+60, sp, sp);
      } else {
        y = alm[i][j] + adif[i][j];
        k = floor((33.0*y-alpha)/(beta-alpha));
        if (k>32) {
          k=32;
        }

        XSetForeground(td,tgc,clron[k].pixel);
        XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j*sp+30, i*sp+60, sp, sp);
      }
    }
  }
}


/**
 * Draw in rings
 */
void picturerings() {
  int i, j, k;

  buildbig();
  for (i = 1; i < nr; i++) {
    for (j = 1; j < nc; j++) {
      if (apic[i][j] == 0) {
        k = floor(63.0 * (adif[i][j]/(rho)));
        XSetForeground(td, tgc, clroff[k].pixel);
        XFillRectangle(te.xexpose.display, te.xexpose.window, tgc,j*sp+30,i*sp+60,sp,sp);
      } else {
        k = ash[i][j];
        k = k % KAPPAMAX;
        XSetForeground(td, tgc, clr[k].pixel);
        XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j*sp+30, i*sp+60, sp, sp);
        if (alm[i][j]>1+0.5*(beta-1.0)) {
          if (alm[i][j] >= 1+0.2 * (beta-1.0)) k=12;
          if (alm[i][j] >= 1+0.5 * (beta-1.0)) k=13;
          if (alm[i][j] >= 1+0.7 * (beta-1.0)) k=14;
          if (alm[i][j] >= beta)               k=15;

          XSetForeground(td,tgc,othp[k].pixel);
          XFillRectangle(te.xexpose.display, te.xexpose.window, tgc, j*sp+30, i*sp+60, sp, sp);
        }
      }
    }
  }
}


/**
 * Draw buttons
 */
void drawbuttons() {
  char quitstring[]  ="QUIT";
  char pausestring[] ="STOP";
  char playstring[]  ="START";

  XSetForeground(td, tgc, tf);
  XSetBackground(td, tgc, tb);

  XDrawRectangle(te.xexpose.display, te.xexpose.window,tgc,
                 20, 50, nc*sp + 20, nr*sp+20);
  XDrawRectangle(te.xexpose.display, te.xexpose.window,tgc,
                 10, 10, 50, 20);
  XDrawRectangle(te.xexpose.display, te.xexpose.window,tgc,
                 65, 10, 50, 20);
  XDrawRectangle(te.xexpose.display, te.xexpose.window,tgc,
                 120, 10, 50, 20);

  XDrawImageString(te.xexpose.display, te.xexpose.window, tgc,
                   20,25, quitstring, strlen(quitstring));
  XDrawImageString(te.xexpose.display,te.xexpose.window, tgc,
                   75,25, pausestring,strlen(pausestring));
  XDrawImageString(te.xexpose.display, te.xexpose.window, tgc,
                   130,25, playstring,strlen(playstring));
}


/**
 * Sections orientation
 */
void transform(int i, int j, int *i1, int *j1) {
  int x1, y1, z1, n1;
  n1 = nc-2;
  x1 = j-n1;
  y1 = n1-i;
  while ((x1<0)||(y1>0)) {
    if ((y1>0)&&(y1<=x1)) {
      x1 = x1 - y1;
      y1 = -y1;
    } else if ((x1>0)&&(x1<=y1)) {
      z1 = x1;
      x1 = y1;
      y1 = z1;
    } else {
      x1 = -x1;
      y1 = -y1;
    }
  }
  *i1 = -y1 + 1;
  *j1 = x1 + 1;
}


/**
 * Save to graphicsfile
 */
void savesnowflake() {
  int i, j, i1, j1, k;
  double y;

  FILE *picf = fopen(GRAPHICSFILE, "w");
  fprintf(picf, "P3\n");

  fprintf(picf, "#rho:%lf\n", rho);
  fprintf(picf, "#h:%d\n",rinit);
  fprintf(picf, "#p:%lf\n", rhorinit);
  fprintf(picf, "#beta:%lf\n", beta);
  fprintf(picf, "#alpha:%lf\n", alpha );
  fprintf(picf, "#theta:%lf\n",theta);
  fprintf(picf, "#kappa:%lf\n", kappa);
  fprintf(picf, "#mu:%lf\n", mu);
  fprintf(picf, "#gam:%lf\n",gam);
  fprintf(picf, "#sigma:%lf\n",sigma);

  fprintf(picf, "%d %d\n", 2*(nc-2)+1, 2*(nr-2)+1);
  fprintf(picf, "255\n");

  buildbig();

  for (i = 0; i <= 2*(nr-2); i++) {
    for (j = 0; j <= 2*(nc-2); j++) {
      transform(i, j, &i1, &j1);
      if (pq % 2 == 1) {
        if (apic[i1][j1] == 0) {
          k = floor(63.0*(adif[i1][j1]/(rho)));
          fprintf(picf, "%d %d %d ",
                  clroff[k].red   * 255 / 65535,
                  clroff[k].green * 255 / 65535,
                  clroff[k].blue  * 255 / 65535);
        } else {
          y = alm[i1][j1]+adif[i1][j1];
          k = floor((33.0*y-alpha)/(beta-alpha)); if (k>32) k=32;
          fprintf(picf, "%d %d %d ",
                  clron[k].red   * 255 / 65535,
                  clron[k].green * 255 / 65535,
                  clron[k].blue  * 255 / 65535);

        }
      } else {
        if (apic[i1][j1] == 0) {
          k=floor(63.0*(adif[i1][j1]/(rho)));
          fprintf(picf, "%d %d %d ",
                  clroff[k].red   * 255 / 65535,
                  clroff[k].green * 255 / 65535,
                  clroff[k].blue  * 255 / 65535);
        } else {
          if (alm[i1][j1]>1+0.5*(beta-1.0)) {
            if (alm[i1][j1] >= 1+0.2*(beta-1.0)) k=12;
            if (alm[i1][j1] >= 1+0.5*(beta-1.0)) k=13;
            if (alm[i1][j1] >= 1+0.7*(beta-1.0)) k=14;
            if (alm[i1][j1] >= beta)             k=15;
            fprintf(picf, "%d %d %d ",
                    othp[k].red   * 255 / 65535,
                    othp[k].green * 255 / 65535,
                    othp[k].blue  * 255 / 65535);
          } else {
            k = ash[i1][j1];
            k = k % KAPPAMAX;
            fprintf(picf,"%d %d %d ",
                    clr[k].red   * 255 / 65535,
                    clr[k].green * 255 / 65535,
                    clr[k].blue  * 255 / 65535);
          }
        }
      }
    }

    fprintf(picf, "\n");
  }

  fclose(picf);
}


double myrand(double min, double max) {
  double f = (double)rand() / RAND_MAX;
  return min + f * (max - min);
}


/**
 * Parameter initialization
 */
void set_parameters() {
  srand(time(NULL));

  rho      = myrand(RHO_MIN, RHO_MAX);
  rinit    = RINIT;
  rhorinit = RHORINIT;
  beta     = myrand(BETA_MIN, BETA_MAX);
  alpha    = myrand(ALPHA_MIN, ALPHA_MAX);
  theta    = myrand(THETA_MIN, THETA_MAX);
  kappa    = myrand(KAPPA_MIN, KAPPA_MAX);
  mu       = myrand(MU_MIN, MU_MAX);
  gam      = GAM;
  sigma    = myrand(SIGMA_MIN, SIGMA_MAX);
  nr       = NR;
  nc       = nr;
  sp       = SP;
}


/**
 * Graphics initialization
 */
void init_graphics() {
  int i;

  td = XOpenDisplay("");
  ts = DefaultScreen(td);
  tb = XWhitePixel(td, ts);
  tf = XBlackPixel(td, ts);

  th.x = 0;
  th.y = 0;

  th.width = nc * sp + 100;
  th.height = nr * sp + 100;

  th.flags = PPosition | PSize;

  tw = XCreateSimpleWindow(td, DefaultRootWindow(td),th.x,th.y, th.width, th.height,7, tf,tb);
  XSetWindowBorderWidth(td, tw, 100);
  XSetStandardProperties(td, tw, "digital snowflake", "sn", None, NULL, 0, &th);
  cmap = DefaultColormap(td, ts);

  braquecolors64();

  for (i = 0; i < KAPPAMAX; i++) {
    clr[i].red   = red[i]   * 65535 / 255;
    clr[i].green = green[i] * 65535 / 255;
    clr[i].blue  = blue[i]  * 65535 / 255;
    XAllocColor(td, cmap, &clr[i]);
  }

  bluecolors33();

  for (i = 0; i <= 32; i++) {
    clron[i].red   = red[i]   * 65535 / 255;
    clron[i].green = green[i] * 65535 / 255;
    clron[i].blue  = blue[i]  * 65535 / 255;
    XAllocColor(td, cmap, &clron[i]);
  }

  offcolors();

  for (i = 0; i <= 63; i++) {
    clroff[63-i].red   = red[i]   * 65535 / 255;
    clroff[63-i].green = green[i] * 65535 / 255;
    clroff[63-i].blue  = blue[i]  * 65535 / 255;
    XAllocColor(td, cmap, &clroff[63-i]);
  }

  XAllocNamedColor(td, cmap, "orange",         &othp[0],  &othp[0]);
  XAllocNamedColor(td, cmap, "gray90",         &othp[1],  &othp[1]);
  XAllocNamedColor(td, cmap, "gray80",         &othp[2],  &othp[2]);
  XAllocNamedColor(td, cmap, "gray70",         &othp[3],  &othp[3]);
  XAllocNamedColor(td, cmap, "gray60",         &othp[4],  &othp[4]);
  XAllocNamedColor(td, cmap, "gray50",         &othp[5],  &othp[5]);
  XAllocNamedColor(td, cmap, "gray40",         &othp[6],  &othp[6]);
  XAllocNamedColor(td, cmap, "gray30",         &othp[7],  &othp[7]);
  XAllocNamedColor(td, cmap, "gray25",         &othp[8],  &othp[8]);
  XAllocNamedColor(td, cmap, "gray20",         &othp[9],  &othp[9]);
  XAllocNamedColor(td, cmap, "black",          &othp[10], &othp[10]);
  XAllocNamedColor(td, cmap, "azure",          &othp[11], &othp[11]);
  XAllocNamedColor(td, cmap, "lightblue2",     &othp[12], &othp[12]);
  XAllocNamedColor(td, cmap, "lightblue3",     &othp[13], &othp[13]);
  XAllocNamedColor(td, cmap, "lightblue4",     &othp[14], &othp[14]);
  XAllocNamedColor(td, cmap, "cornflowerblue", &othp[15], &othp[15]);
  XAllocNamedColor(td, cmap, "white",          &othp[16], &othp[16]);
  XAllocNamedColor(td, cmap, "palegreen",      &othp[17], &othp[17]);
  XAllocNamedColor(td, cmap, "red",            &othp[18], &othp[18]);

  tgc = XCreateGC(td, tw, 0, 0);

  XSetBackground(td, tgc, tb);

  XSelectInput(td, tw, (ButtonPressMask|ExposureMask));

  XMapRaised(td, tw);

  XNextEvent(td, &te);
  drawbuttons();

  initialize();
  picturebig();
}

/*************************************************************************************************/



/*** MAIN ****************************************************************************************/

int main(int argc, char *argv[]) {
  set_parameters();
  init_graphics();

  Window rw, cw;
  int posx, posy;
  int rootx, rooty;
  unsigned int kgb;

  pq = 0;

  bool keep_running = true;
  while (keep_running) {
    XNextEvent(td, &te);
    switch (te.type) {
      case ButtonPress:
        XQueryPointer(td, tw, &rw, &cw, &rootx, &rooty, &posx, &posy, &kgb);

        /* QUIT */
        if ((posx >= 10) && (posx <= 60) && (posy >= 10) && (posy <= 30)) {
          savesnowflake();
          keep_running = false;
        }

        /* STOP */
        else if ((posx >= 65) && (posx <= 115) && (posy >= 10) && (posy <= 30)) {
          if (pq % 2 == 0) {
            picturerings();
          } else {
            picturebig();
          }
        }

        /* START */
        else if ((posx>=120) && (posx<=170) && (posy>=10) && (posy<=30)) {
          while((XEventsQueued(td, QueuedAfterReading) == 0) && (pq!=-1) && (stop == false)) {
            noac=0;
            pq++;
            dynamics();
            if (pq % 10 == 0) {
              picturebig();
            }
          }
        }

        break;

      case Expose:
        if (te.xexpose.count == 0) {
          picturebig();
          drawbuttons();
        }

        break;
    }
  }

  /* Clean graphics */
  XFreeGC(td, tgc);
  XDestroyWindow(td, tw);
  XCloseDisplay(td);

  return 0;
}

/*************************************************************************************************/
