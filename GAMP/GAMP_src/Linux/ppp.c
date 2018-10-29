/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*-----------------------------------------------------------------------------*/
#include "gamp.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

const static int glo_fcn[] = {-7,-6,-4,-3,-2,-1,0,1,2,3,4,5,6}; /* glonass frequency channel number */

/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
	double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
	int i,j;
	const char *type;
	double sec,ex[3]={0.0,0.0,0.0};
	int week;

	/* unit vector of sun direction (ecef) */
	sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
	normv3(rsun,esun);

	for (i=0;i<n;i++) {
		sec=time2gpst(obs[0].time,&week);
		week=calEclips(obs[i].sat,rs+i*6,rs+i*6+3,rsun,sec,ex,nav);

		if (week && 1) {
			sprintf(PPP_Glo.chMsg,"*** WARNING: %s ecliType=%d\n",PPP_Glo.sFlag[obs[i].sat-1].id,week);
			outDebug(OUTWIN,OUTFIL,OUTTIM);
			PPP_Glo.ecliF[obs[i].sat-1]=4.0;
		}
		else
			PPP_Glo.ecliF[obs[i].sat-1]=1.0;

		type=nav->pcvs[obs[i].sat-1].type;

		if ((r=norm(rs+i*6,3))<=0.0) continue;

		/* only block IIA */
		if (*type&&!strstr(type,"BLOCK IIA")) continue;

		/* sun-earth-satellite angle */
		cosa=dot(rs+i*6,esun,3)/r;
		cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
		ang=acos(cosa);

		/* test eclipse */
		if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

		sprintf(PPP_Glo.chMsg,"*** WARNING: eclipsing sat excluded sat=%s\n",PPP_Glo.sFlag[obs[i].sat-1].id);
		outDebug(OUTWIN,OUTFIL,OUTTIM);

		for (j=0;j<3;j++) rs[j+i*6]=0.0;
	}
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
	if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
	return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
	double *yaw)
{
	*yaw=yaw_nominal(beta,mu);
	return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
	const double *rs, double *exs, double *eys)
{
	double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
	double yaw,cosy,siny,erpv[5]={0};
	int i;

	sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);

	/* beta and orbit angle */
	matcpy(ri,rs,6,1);
	ri[3]-=OMGE*ri[1];
	ri[4]+=OMGE*ri[0];
	cross3(ri,ri+3,n);
	cross3(rsun,n,p);
	if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
		!normv3(p,ep)) return 0;
	beta=PI/2.0-acos(dot(esun,en,3));
	E=acos(dot(es,ep,3));
	mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
	if      (mu<-PI/2.0) mu+=2.0*PI;
	else if (mu>=PI/2.0) mu-=2.0*PI;

	/* yaw-angle of satellite */
	if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;

	/* satellite fixed x,y-vector */
	cross3(en,es,ex);
	cosy=cos(yaw);
	siny=sin(yaw);
	for (i=0;i<3;i++) {
		exs[i]=-siny*en[i]+cosy*ex[i];
		eys[i]=-cosy*en[i]-siny*ex[i];
	}
	return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
	const double *rs, const double *rr, double *phw)
{
	double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
	double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
	int i;

	if (opt<=0) return 1; /* no phase windup */

	if (norm(rr,3)<=0.0) return 0;

	/* satellite yaw attitude model */
	if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;

	/* unit vector satellite to receiver */
	for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
	if (!normv3(r,ek)) return 0;

	/* unit vectors of receiver antenna */
	ecef2pos(rr,pos);
	xyz2enu(pos,E);
	exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
	eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

	/* phase windup effect */
	cross3(ek,eys,eks);
	cross3(ek,eyr,ekr);
	for (i=0;i<3;i++) {
		ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
		dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
	}
	cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
	if      (cosp<-1.0) cosp=-1.0;
	else if (cosp> 1.0) cosp= 1.0;
	//acos£¨-1£© invalid
	if (fabs(fabs(cosp)-1.0)<1.0e-10) return 0;

	ph=acos(cosp)/2.0/PI;
	cross3(ds,dr,drs);
	if (dot(ek,drs,3)<0.0) ph=-ph;

	*phw=ph+floor(*phw-ph+0.5); /* in cycle */
	return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int freq, int type,
	const prcopt_t *opt)
{
	double a=opt->err[1],b=opt->err[2];
	double c=1.0,fact=1.0;
	double sinel=sin(el);

	fact=EFACT_GPS;
	c=type?opt->err[0]:1.0;   /* type=0:phase,1:code */

	if (sys==SYS_GLO) {
		fact=EFACT_GLO;
		if (type) c=PPP_Glo.prcOpt_Ex.errRatioGLO;
	}
	else if (sys==SYS_CMP) {
		if (type) c=PPP_Glo.prcOpt_Ex.errRatioBDS;
	}
	else if (sys==SYS_GAL) {
		if (type) c=PPP_Glo.prcOpt_Ex.errRatioGAL;
	}
	else if (sys==SYS_QZS) {
		if (type) c=PPP_Glo.prcOpt_Ex.errRatioQZS;
	}

	if (opt->ionoopt==IONOOPT_IF12) fact*=3.0;

	return SQR(fact*c)*(SQR(a)+SQR(b/sinel));
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
	int j;
	rtk->x[i]=xi;
	for (j=0;j<rtk->nx;j++) {
		rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
	}
}
/* L1/L2 geometry-free phase measurement -------------------------------------*/
extern double gfmeas(const obsd_t *obs, const nav_t *nav)
{
	const double *lam=nav->lam[obs->sat-1];

	if (lam[0]==0.0||lam[1]==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;

	return lam[0]*obs->L[0]-lam[1]*obs->L[1];
}
/* L1/L2 geometry-free phase measurement considering PCO, PCV and phase windup -------------------------------------*/
static double gfmeas_new(const obsd_t *obs, const nav_t *nav)
{
	const double *lam=nav->lam[obs->sat-1];
	double L[NFREQ]={0},gf;
	int i,sat;

	if (lam[0]==0.0||lam[1]==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;

	for (i=0;i<NFREQ;i++) {
		/* antenna phase center and phase windup correction */
		L[i]=obs->L[i]*lam[i];
	}

	sat=obs->sat;
	for (i=0;i<NFREQ;i++) {
		/* antenna phase center and phase windup correction */
		L[i]=obs->L[i]*lam[i]-PPP_Glo.ssat_Ex[sat-1].dantr[i]-PPP_Glo.ssat_Ex[sat-1].phw*lam[i];
	}
	gf=L[0]-L[1];

	return gf;
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
	const prcopt_t *opt, const double *dantr,
	const double *dants, double phw, double *L, double *P,
	double *Lc, double *Pc)
{
	const double *lam=nav->lam[obs->sat-1];
	double C1,C2;
	int i,sys;

	for (i=0;i<NFREQ;i++) {
		L[i]=P[i]=0.0;
		if (lam[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
		if (testsnr(0,0,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) continue;

		/* antenna phase center and phase windup correction */
		L[i]=obs->L[i]*lam[i]-dants[i]-dantr[i]-phw*lam[i];
		P[i]=obs->P[i]       -dants[i]-dantr[i];

		/* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
		if (obs->code[i]==CODE_L1C) {
			P[i]+=nav->cbias[obs->sat-1][1];
		}
		else if (obs->code[i]==CODE_L2C||obs->code[i]==CODE_L2X||
			obs->code[i]==CODE_L2L||obs->code[i]==CODE_L2S) {
				P[i]+=nav->cbias[obs->sat-1][2];
#if 0
				L[i]-=0.25*lam[i]; /* 1/4 cycle-shift */
#endif
		}
	}
	/* iono-free LC */
	*Lc=*Pc=0.0;
	sys=PPP_Glo.sFlag[obs->sat-1].sys;
	i=1;
	if (lam[0]==0.0||lam[i]==0.0) return;

	C1= SQR(lam[i])/(SQR(lam[i])-SQR(lam[0]));
	C2=-SQR(lam[0])/(SQR(lam[i])-SQR(lam[0]));

	/* P1-P2 dcb correction (P1->Pc,P2->Pc) */
	if (P[0]!=0.0) P[0]-=C2*nav->cbias[obs->sat-1][0];
	if (P[1]!=0.0) P[1]+=C1*nav->cbias[obs->sat-1][0];

	if (L[0]!=0.0&&L[i]!=0.0) *Lc=C1*L[0]+C2*L[i];
	if (P[0]!=0.0&&P[i]!=0.0) *Pc=C1*P[0]+C2*P[i];
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk)
{
	int i;

	/* fixed mode */
	if (rtk->opt.mode==PMODE_PPP_FIXED) {
		for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
		return;
	}
	/* initialize position for first epoch */
	if (norm(rtk->x,3)<=0.0) {
		for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
	}
	/* static ppp mode */
	if (rtk->opt.mode==PMODE_PPP_STATIC) {
		/*for (i=0;i<3;i++) {
		rtk->P[i*(1+rtk->nx)]+=SQR(rtk->opt.prn[5])*fabs(rtk->tt);
		}*/
		return;
	}
	/* kinmatic mode without dynamics */
	for (i=0;i<3;i++) {
		initx(rtk,rtk->sol.rr[i],VAR_POS,i);
	}
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	prcopt_t *opt=&rtk->opt;
	double dtr,tt=fabs(rtk->tt),tdif,isb_est,isb_conv;
	int i,ii,jj,sat_r,frq,ic,isys,sys=rtk->opt.navsys;

	if (sys==SYS_GLO) {
		dtr=rtk->sol.dtr[1];
		if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
		ic=IC(1,opt);
		initx(rtk,CLIGHT*dtr,VAR_CLK,ic);

		if (opt->gloicb==GLOICB_LNF) {   //linear function of frequency number
			ic=IICB(1,opt);
			if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
		}
		else if (opt->gloicb==GLOICB_QUAD) {   //quadratic polynomial function of frequency number
			ic=IICB(1,opt);
			if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
			ic=IICB(2,opt);
			if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
		}
		else if (opt->gloicb==GLOICB_1SAT) {   //isb+icb every sat
			for (ii=0;ii<n;ii++) {
				sat_r=obs[ii].sat-NSATGPS;
				ic=IICB(sat_r,opt);
				if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
			}
		}
		else if (opt->gloicb==GLOICB_1FRQ) {   //isb+icb every frq
			for (ii=0;ii<n;ii++) {
				frq=get_glo_fcn(obs[ii].sat,nav);
				for (jj=1;jj<14;jj++) {
					if (glo_fcn[jj-1]==frq) break;
				}
				ic=IICB(jj,opt);
				if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
			}
		}

		return;
	}
	else if (sys==SYS_CMP) {
		dtr=rtk->sol.dtr[2];
		if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
		ic=IC(2,opt);
		initx(rtk,CLIGHT*dtr,VAR_CLK,ic);

		return;
	}
	else if (sys==SYS_GAL) {
		dtr=rtk->sol.dtr[3];
		if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
		ic=IC(3,opt);
		initx(rtk,CLIGHT*dtr,VAR_CLK,ic);

		return;
	}

	/* initialize every epoch for clock (white noise) */
	dtr = rtk->sol.dtr[0];
	if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
	ic=IC(0,opt);
	initx(rtk,CLIGHT*dtr,VAR_CLK,ic);

	for (i=1;i<NSYS;i++) {
		if (!(sys&SYS_GLO)&&i==1) continue;
		if (!(sys&SYS_CMP)&&i==2) continue;
		if (!(sys&SYS_GAL)&&i==3) continue;
		if (!(sys&SYS_QZS)&&i==4) continue;

		if (sys&SYS_GLO&&i==1) {
			if (opt->gloicb==GLOICB_OFF||opt->gloicb==GLOICB_LNF||opt->gloicb==GLOICB_QUAD) {  //handling of ISBs
				dtr=rtk->sol.dtr[1];
                ic=IC(1,opt);

				if (opt->gnsisb==GNSISB_CT) {
					// constant
					if (rtk->x[ic]==0.0) {
						if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
						initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
					}
				}
				else if (opt->gnsisb==GNSISB_PWC) {
					// piece-wise constant
					if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
					if (PPP_Glo.iEpoch==1) initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
					tdif=timediff(PPP_Glo.t_30min,PPP_Glo.tNow);
					if (tdif==0.0) {
						isb_est=rtk->x[ic];
						isb_conv=rtk->P[ic+ic*rtk->nx];
						initx(rtk,isb_est,VAR_CLK,ic);
						PPP_Glo.isb_30min=isb_est/CLIGHT*1.0e+9;
					}
				}
				else if (opt->gnsisb==GNSISB_RW) {
					// random walk process
					if (rtk->x[ic]==0.0) {
						if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
						initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
					}
					else {
						rtk->P[ic+ic*rtk->nx]+=SQR(0.001)*tt;
					}
				}
				else if (opt->gnsisb==GNSISB_WN) {
					//white noise process
					if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
					initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
				}
			}

			if (opt->gloicb==GLOICB_LNF) {   //linear function of frequency number
				ic=IICB(1,opt);
				if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
			}
			else if (opt->gloicb==GLOICB_QUAD) {   //quadratic polynomial function of frequency number
				ic=IICB(1,opt);
				if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
				ic=IICB(2,opt);
				if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
			}
			else if (opt->gloicb==GLOICB_1SAT) {   //isb+icb every sat
                dtr=rtk->sol.dtr[1];
                ic=IC(1,opt);
				initx(rtk,0.0,0.0,ic);
				for (ii=0;ii<n;ii++) {
					isys=PPP_Glo.sFlag[obs[ii].sat-1].sys;
					if (isys!=SYS_GLO) continue;
					sat_r=obs[ii].sat-NSATGPS;
					ic=IICB(sat_r,opt);

					if (opt->gnsisb==GNSISB_CT) {
						//constant
						if (rtk->x[ic]==0.0) initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
					}
                    else if (opt->gnsisb==GNSISB_RW) {
                        // random walk process
                        if (rtk->x[ic]==0.0) {
                            if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
                            initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
                        }
                        else {
                            rtk->P[ic+ic*rtk->nx]+=SQR(0.001)*tt;
                        }
                    }
					else if (opt->gnsisb==GNSISB_WN) {
						//white noise process
						initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
					}
				}
			}
			else if (opt->gloicb==GLOICB_1FRQ) {   //isb+icb every frq
				for (ii=0;ii<n;ii++) {
					isys=PPP_Glo.sFlag[obs[ii].sat-1].sys;
					if (isys!=SYS_GLO) continue;
					frq=get_glo_fcn(obs[ii].sat,nav);
					for (jj=1;jj<14;jj++) {
						if (glo_fcn[jj-1]==frq) break;
					}
					ic=IICB(jj,opt);
					//if (rtk->x[ic]==0.0) initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
					if (rtk->x[ic]==0.0) initx(rtk,0.1,VAR_CLK,ic);
				}
			}
		}
		else {
			dtr=rtk->sol.dtr[i];
			ic=IC(i,opt);

			if (opt->gnsisb==GNSISB_CT) {
				// constant
				if (rtk->x[ic]==0.0) {
					if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
					initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
				}
			}
			else if (opt->gnsisb==GNSISB_PWC) {
				// piece-wise constant
				if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
				if (PPP_Glo.iEpoch==1) initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
				tdif=timediff(PPP_Glo.t_30min,PPP_Glo.tNow);
				if (tdif==0.0) {
					isb_est=rtk->x[ic];
					isb_conv=rtk->P[ic+ic*rtk->nx];
					initx(rtk,isb_est,VAR_CLK,ic);
					PPP_Glo.isb_30min=isb_est/CLIGHT*1.0e+9;
				}
			}
			else if (opt->gnsisb==GNSISB_RW) {
				// random walk process
				if (rtk->x[ic]==0.0) {
					if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
					initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
				}  
				else {
					rtk->P[ic+ic*rtk->nx]+=SQR(0.001)*tt;
				}
			}
			else if (opt->gnsisb==GNSISB_WN) {
				//white noise process
				if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
				initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
			}
		}
	}
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
	double var,tt=fabs(rtk->tt);
	int i=IT(&rtk->opt),j;

	if (rtk->x[i]==0.0) {
		var=SQR(0.3);
		initx(rtk,0.15,var,i);

		if (rtk->opt.tropopt>=TROPOPT_ESTG) {
			for (j=i+1;j<i+3;j++) initx(rtk,1E-6,VAR_GRA,j);
		}
	}
	else {
		rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[2])*tt;

		if (rtk->opt.tropopt>=TROPOPT_ESTG) {
			for (j=i+1;j<i+3;j++) {
				rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2]*0.1)*tt;
			}
		}
	}
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const double *lam;
	double ion,g0,g1,dg,dion,elev,sinel,dion_tec,vari_tec,pos[3],tt=fabs(rtk->tt);
	char *p;
	int i,j,k,sat,slip[MAXOBS]={0},gap_resion=GAP_RESION;

	if ((p=strstr(rtk->opt.pppopt,"-GAP_RESION="))) {
		sscanf(p,"-GAP_RESION=%d",&gap_resion);
	}
	for (i=0;i<MAXSAT;i++) {
		j=II(i+1,&rtk->opt);
		if (rtk->x[j]!=0.0&&(int)rtk->ssat[i].outc[0]>gap_resion) {
			rtk->x[j]=0.0;
		}
	}
	ecef2pos(rtk->sol.rr,pos);
	for (i=0;i<n;i++) {
		sat=obs[i].sat;
		j=II(sat,&rtk->opt);
		lam=nav->lam[sat-1];
		elev=rtk->ssat[sat-1].azel[1];
		sinel=sin(elev);
		elev*=R2D;
		if (PPP_Glo.prcOpt_Ex.ionopnoise==3||((PPP_Glo.prcOpt_Ex.ionopnoise==1||
			PPP_Glo.prcOpt_Ex.ionopnoise==2)&&rtk->x[j]==0.0)) {
				if (PPP_Glo.prcOpt_Ex.ion_const) {
					//ionospheric delay derived from GIM
					iontec(obs[i].time,nav,pos,rtk->ssat[sat-1].azel,1,&dion_tec,&vari_tec);
					ion=dion_tec;
				}
				else {
					k=1;
					if (obs[i].P[0]==0.0||obs[i].P[k]==0.0||lam[0]==0.0||lam[k]==0.0) {
						continue;
					}
					ion=(obs[i].P[0]-obs[i].P[k])/(1.0-SQR(lam[k]/lam[0]));
				}
				
				initx(rtk,ion,VAR_IONO,j);

				if ((g1=gfmeas_new(obs+i,nav))==0.0) continue;
				PPP_Glo.ssat_Ex[sat-1].gf=g1;
		}
		else if ((PPP_Glo.prcOpt_Ex.ionopnoise==1||PPP_Glo.prcOpt_Ex.ionopnoise==2)&&rtk->x[j]!=0.0) {
			if (PPP_Glo.prcOpt_Ex.ionopnoise==1)
				rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1])*tt;
			else if (PPP_Glo.prcOpt_Ex.ionopnoise==2) {
				slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
				if (!slip[i]) {
					if ((g1=gfmeas_new(obs+i,nav))==0.0) continue;
					g0=PPP_Glo.ssat_Ex[sat-1].gf;
					PPP_Glo.ssat_Ex[sat-1].gf=g1;
					dg=g1-g0;
					dion=dg*SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));
					rtk->x[j]+=dion;
					if (elev>=30.0)
						rtk->P[j+j*rtk->nx]+=SQR(PPP_Glo.prcOpt_Ex.prn_iono)*tt;
					else
						rtk->P[j+j*rtk->nx]+=SQR(PPP_Glo.prcOpt_Ex.prn_iono/(2*sinel))*tt;
				}
				else
					if (elev>=30.0)
						rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1])*tt;
					else
						rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]/(2*sinel))*tt;
			}
		}
	}
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t *rtk)
{
	int i=ID(&rtk->opt);

	if (rtk->x[i]==0.0) {
		initx(rtk,1E-6,VAR_DCB,i);
	}
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const double *lam;
	double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
	double ion,dantr[NFREQ]={0},dants[NFREQ]={0},tt=fabs(rtk->tt);
	int i,j,k,l,f,sat,slip[MAXOBS]={0},clk_jump=0;

	/* handle day-boundary clock jump */
	clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;

	ecef2pos(rtk->sol.rr,pos);

	for (f=0;f<NF(&rtk->opt);f++) {
		/* reset phase-bias if expire obs outage counter */
		for (i=0;i<MAXSAT;i++) {
			if (++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout||
				rtk->opt.modear==ARMODE_INST||clk_jump) {
					initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));
			}
		}
		for (i=k=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			j=IB(sat,f,&rtk->opt);
			corr_meas(obs+i,nav,rtk->ssat[sat-1].azel,&rtk->opt,dantr,dants,
				0.0,L,P,&Lc,&Pc);

			bias[i]=0.0;

			if (rtk->opt.ionoopt==IONOOPT_IF12) {
				bias[i]=Lc-Pc;
				slip[i]=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
			}
			else if (L[f]!=0.0&&P[f]!=0.0) {
				slip[i]=rtk->ssat[sat-1].slip[f];
				l=1;
				lam=nav->lam[sat-1];
				if (obs[i].P[0]==0.0||obs[i].P[l]==0.0||
					lam[0]==0.0||lam[l]==0.0||lam[f]==0.0) continue;
				ion=(obs[i].P[0]-obs[i].P[l])/(1.0-SQR(lam[l]/lam[0]));
				bias[i]=L[f]-P[f]+2.0*ion*SQR(lam[f]/lam[0]);
			}
			if (rtk->x[j]==0.0||slip[i]||bias[i]==0.0) continue;

			offset+=bias[i]-rtk->x[j];
			k++;
		}
		/* correct phase-code jump to ensure phase-code coherency */
		if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
			for (i=0;i<MAXSAT;i++) {
				j=IB(i+1,f,&rtk->opt);
				if (rtk->x[j]!=0.0) rtk->x[j]+=offset/k;
			}
			printf("phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
				time_str(rtk->sol.time,0),k,offset/k/CLIGHT);
		}
		for (i=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			j=IB(sat,f,&rtk->opt);

			rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[0])*tt;

			if (bias[i]==0.0||(rtk->x[j]!=0.0&&!slip[i])) continue;

			/* reinitialize phase-bias if detecting cycle slip */
			initx(rtk,bias[i],VAR_BIAS,IB(sat,f,&rtk->opt));

			/* reset fix flags */
			for (k=0;k<MAXSAT;k++) rtk->ambc[sat-1].flags[k]=0;
		}
	}
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	/* temporal update of position */
	udpos_ppp(rtk);

	/* temporal update of clock */
	udclk_ppp(rtk,obs,n,nav);

	/* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
		udtrop_ppp(rtk);
	}
	/* temporal update of ionospheric parameters */
	if (rtk->opt.ionoopt==IONOOPT_UC1||rtk->opt.ionoopt==IONOOPT_UC12) {
		udiono_ppp(rtk,obs,n,nav);
	}
	/* temporal update of L5-receiver-dcb parameters */
	if (rtk->opt.nf>=3||rtk->opt.ionoopt==IONOOPT_UC12) {
	//if (rtk->opt.nf>=3) {
		uddcb_ppp(rtk);
	}
	/* temporal update of phase-bias */
	udbias_ppp(rtk,obs,n,nav);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(int sat, const double *rs, const double *rr, const pcv_t *pcv,
	double *dant)
{
	double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
	int i;

	for (i=0;i<3;i++) {
		ru[i]=rr[i]-rs[i];
		rz[i]=-rs[i];
	}
	if (!normv3(ru,eu)||!normv3(rz,ez)) return;

	cosa=dot(eu,ez,3);
	cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
	nadir=acos(cosa);

	antmodel_s(sat,pcv,nadir,dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos, const double *azel, 
	const prcopt_t *opt, const double *x, double *dtdx, double *shd, double *var)
{
	const double zazel[]={0.0,PI/2.0};
	double zhd,zwd,m_h,m_w,cotz,grad_n,grad_e;

	/* zenith hydrostatic delay */
	zhd=tropmodel(time,pos,zazel,0.0,&zwd,1);
    PPP_Glo.zhd=zhd;

	/* mapping function */
	m_h=tropmapf(time,pos,azel,&m_w);

	*shd=m_h*zhd;

	if (opt->tropopt>=TROPOPT_ESTG&&azel[1]>0.0) {
		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
		cotz=1.0/tan(azel[1]);
		grad_n=m_w*cotz*cos(azel[0]);
		grad_e=m_w*cotz*sin(azel[0]);
		m_w+=grad_n*x[1]+grad_e*x[2];
		dtdx[1]=grad_n*(x[0]);
		dtdx[2]=grad_e*(x[0]);
	}
	dtdx[0]=m_w;
	*var=SQR(0.01);

	//return m_h*zhd+m_w*(x[0]-zhd);
	return m_h*zhd+m_w*(x[0]);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
	const prcopt_t *opt, const double *x, double *dtdx,
	const nav_t *nav, double *dtrp, double *shd, double *var)
{
	double trp[3]={0};
	double zwd;

	if (opt->tropopt==TROPOPT_SAAS) {
		*dtrp=tropmodel(time,pos,azel,REL_HUMI,&zwd,0);
		*dtrp+=zwd;
		*var=SQR(ERR_SAAS);
		return 1;
	}
	if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
		matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
		*dtrp=trop_model_prec(time,pos,azel,opt,trp,dtdx,shd,var);
		return 1;
	}
	return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
	const prcopt_t *opt, int sat, const double *x,
	const nav_t *nav, double *dion, double *var)
{    
	if (opt->ionoopt==IONOOPT_TEC) {
		return iontec(time,nav,pos,azel,1,dion,var);
	}
	if (opt->ionoopt==IONOOPT_BRDC) {
		*dion=ionmodel(time,nav->ion_gps,pos,azel);
		*var=SQR(*dion*ERR_BRDCI);
		return 1;
	}
	if (opt->ionoopt==IONOOPT_UC1||opt->ionoopt==IONOOPT_UC12) {
		*dion=x[II(sat,opt)];
		*var=0.0;
		return 1;
	}
	if (opt->ionoopt==IONOOPT_IF12) {
		*dion=*var=0.0;
		return 1;
	}
	return 0;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
	const double *dts, const double *var_rs, const int *svh,
	int *exc, const nav_t *nav, const double *x, rtk_t *rtk, 
	double *v, double *H, double *R, double *azel)
{
	const double *lam;
	prcopt_t *opt=&rtk->opt;
	double y,r,cdtr,bias,C,rr[3],pos[3],e[3],dtdx[3]={0},dr[3]={0},L[NFREQ],P[NFREQ],Lc,Pc;
	double var[MAXOBS*4],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,shd=0.0,dcb,gravitationalDelayModel=0.0;
	double dion_tec=0.0,vari_tec=0.0,C1,C2;
	double dantr[NFREQ]={0},dants[NFREQ]={0};
	double ve[MAXOBS*2*NFREQ]={0},vmax=0;
	char str[32];
	int ne=0,obsi[MAXOBS*2*NFREQ]={0},frqi[MAXOBS*2*NFREQ],maxobs,maxfrq,rej;
	int i,j,k,ii,ic=0,id=0,jj,kk,sat,sat_r,prn,sys,frq,nv=0,nx=rtk->nx,stat=1,nx_nv;

	time2str(obs[0].time,str,2);

	for (i=0;i<MAXSAT;i++) {
		for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;
		rtk->ssat[i].azel[0]=rtk->ssat[i].azel[1]=0.0;

		if (!post) {
			for (j=0;j<NFREQ;j++) {
				rtk->ssat[i].resc_pri[j]=0.0;
				rtk->ssat[i].resp_pri[j]=0.0;
			}
		}
		else {
			for (j=0;j<NFREQ;j++) {
				rtk->ssat[i].resc_pos[j]=0.0;
				rtk->ssat[i].resp_pos[j]=0.0;
			}
		}
	}

	for (i=0;i<3;i++) rr[i]=x[i];
	if (norm(rr,3)<=100.0) return 0;
	/* earth tides correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time),rr,opt->tidecorr==1?1:7,&nav->erp,
			opt->odisp[0],dr);
		for (i=0;i<3;i++) rr[i]+=dr[i];
	}
	
	ecef2pos(rr,pos);

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		lam=nav->lam[sat-1];
		if (lam[j/2]==0.0||lam[0]==0.0) continue;

		if ((r=geodist(rs+i*6,rr,e))<=0.0||
			satazel(pos,e,azel+i*2)<opt->elmin) {
				exc[i]=1;
				continue;
		}
		if (!(sys=PPP_Glo.sFlag[sat-1].sys)||!rtk->ssat[sat-1].vs||
			satexclude(obs[i].sat,svh[i],opt)||exc[i]) {
				exc[i]=1;
				continue;
		}

		/* tropospheric and ionospheric model */
		if (!model_trop(obs[i].time,pos,azel+i*2,opt,x,dtdx,nav,&dtrp,&shd,&vart)||
			!model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
				continue;
		}

		//ionospheric delay derived from GIM
		if (PPP_Glo.prcOpt_Ex.ion_const) {
			iontec(obs[i].time,nav,pos,azel+i*2,1,&dion_tec,&vari_tec);
			dion_tec*=SQR(lam[0]/lam_carr[0]);
		}

		/* satellite and receiver antenna model */
		satantpcv(sat,rs+i*6,rr,nav->pcvs+sat-1,dants);
		antmodel(sat,&opt->pcvr,opt->antdel,azel+i*2,1,dantr);

		/* phase windup model */
		if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,2,rs+i*6,rr,&rtk->ssat[sat-1].phw)) {
			continue;
		}

		//Gravitational delay correction */
		gravitationalDelayModel=gravitationalDelayCorrection(sys,rr,rs+i*6);

		/* corrected phase and code measurements */
		corr_meas(obs+i,nav,azel+i*2,&rtk->opt,dantr,dants,
			rtk->ssat[sat-1].phw,L,P,&Lc,&Pc);
        rtk->ssat[sat-1].PC=Pc;
        rtk->ssat[sat-1].LC=Lc;

		/* stack phase and code residuals {L1,P1,L2,P2,...} */
		for (j=0;j<2*NF(opt);j++) {
			dcb=bias=0.0;

			if (opt->ionoopt==IONOOPT_IF12) {
				if ((y=j%2==0?Lc:Pc)==0.0) continue;
			}
			else {
				if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;
			}

			//C=SQR(lam[j/2]/lam[0])*ionmapf(pos,azel+i*2)*(j%2==0?-1.0:1.0);
			C=SQR(lam[j/2]/lam[0])*(j%2==0?-1.0:1.0);
			C1= SQR(lam[1])/(SQR(lam[1])-SQR(lam[0]));
			C2=-SQR(lam[0])/(SQR(lam[1])-SQR(lam[0]));

			nx_nv=nx*nv;

			for (k=0;k<nx;k++) H[k+nx_nv]=k<3?-e[k]:0.0;

			/* receiver clock */
			cdtr=0.0;
			if (sys==SYS_GPS) {
				ic=IC(0,opt);
				cdtr=x[ic];
				H[ic+nx_nv]=1.0;
			}
			if (sys==SYS_GLO) {
				if (opt->navsys==SYS_GLO||opt->gloicb==GLOICB_OFF||opt->gloicb==GLOICB_LNF||opt->gloicb==GLOICB_QUAD) {  //handling of ISBs
					ic=IC(0,opt);
					id=IC(1,opt);
					cdtr=x[ic]+x[id];
					H[ic+nx_nv]=1.0;
					H[id+nx_nv]=1.0;
				}
				else if (opt->navsys!=SYS_GLO) {
					ic=IC(0,opt);
					cdtr=x[ic];
					H[ic+nx_nv]=1.0;
				}
				if (opt->gloicb==GLOICB_LNF) {   //linear function of frequency number
					if ((opt->nf==2&&(j/2==1&&j%2==1))||(opt->nf==1&&(j/2==0&&j%2==1))) {
						frq=get_glo_fcn(sat,nav);
						ic=IICB(1,opt);
						cdtr+=frq*x[ic];
						H[ic+nx_nv]=frq;
					}
				}
				else if (opt->gloicb==GLOICB_QUAD) {   //quadratic polynomial function of frequency number
					if ((opt->nf==2&&(j/2==1&&j%2==1))||(opt->nf==1&&(j/2==0&&j%2==1))) {
						frq=get_glo_fcn(sat,nav);
						ic=IICB(1,opt);
						id=IICB(2,opt);
						cdtr+=frq*x[ic];
						H[ic+nx_nv]=frq;
						cdtr+=frq*frq*x[id];
						H[id+nx_nv]=frq*frq;
					}
				}
				else if (opt->gloicb==GLOICB_1SAT) {   //isb+icb every sat
					if ((opt->nf==2&&(j/2==1&&j%2==1))||(opt->nf==1&&(j/2==0&&j%2==1))) {
						sat_r=sat-NSATGPS;
						ic=IICB(sat_r,opt);
						cdtr+=x[ic];
						H[ic+nx_nv]=1.0;
					}
				}
				else if (opt->gloicb==GLOICB_1FRQ) {   //isb+icb every frq
					if ((opt->nf==2&&(j/2==1&&j%2==1))||(opt->nf==1&&(j/2==0&&j%2==1))) {
						frq=get_glo_fcn(sat,nav);
						for (ii=1;ii<14;ii++) {
							if (glo_fcn[ii-1]==frq) break;
						}
						ic=IICB(ii,opt);
						cdtr+=x[ic];
						H[ic+nx_nv]=1.0;
					}
				}
			}
			if (sys==SYS_CMP) {
				ic=IC(0,opt);
				id=IC(2,opt);
				cdtr=x[ic]+x[id];
				H[ic+nx_nv]=1.0;
				H[id+nx_nv]=1.0;
			}
			if (sys==SYS_GAL) {
				ic=IC(0,opt);
				id=IC(3,opt);
				cdtr=x[ic]+x[id];
				H[ic+nx_nv]=1.0;
				H[id+nx_nv]=1.0;
			}
			if (sys==SYS_QZS) {
				ic=IC(0,opt);
				id=IC(4,opt);
				cdtr=x[ic]+x[id];
				H[ic+nx_nv]=1.0;
				H[id+nx_nv]=1.0;
			}

			if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
				for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
					ic=IT(opt);
					H[ic+k+nx_nv]=dtdx[k];
				}
			}
			if (opt->ionoopt==IONOOPT_UC1||opt->ionoopt==IONOOPT_UC12) {
				ic=II(sat,opt);
				if (x[ic]==0.0) continue;
				H[ic+nx_nv]=C;

				//if (j%2==1) {  /* receiver-dcb */
				//	if (j/2==0) {  // for P1 observation
				//		dcb+=rtk->x[ID(opt)];
				//		H[ID(opt)+nx_nv]=C2;
				//	}
				//	else {  // for P2 observation
				//		dcb+=rtk->x[ID(opt)];
				//		H[ID(opt)+nx_nv]=-C1;
				//	}
				//}
				if (sys!=SYS_GLO&&PPP_Glo.prcOpt_Ex.ion_const&&(j/2==1&&j%2==1)) {
					dcb+=rtk->x[ID(opt)];
					H[ID(opt)+nx_nv]=1.0;
				}
			}
			if (j/2==2&&j%2==1) {  /* L5-receiver-dcb */
				id=ID(opt);
				dcb+=x[id];
				H[id+nx_nv]=1.0;
			}
			if (j%2==0) {  /* phase bias */
				ic=IB(sat,j/2,opt);
				if ((bias=x[ic])==0.0) continue;
				H[ic+nx_nv]=1.0;
			}
			/* residual */
			v[nv]=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias-gravitationalDelayModel);

			if (j%2==0) rtk->ssat[sat-1].resc[j/2]=v[nv];
			else        rtk->ssat[sat-1].resp[j/2]=v[nv];

			/* variance */
			var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],j/2,j%2,opt)+
				vart+SQR(C)*vari+var_rs[i];
			//if (sys==SYS_GLO&&j%2==1) var[nv]+=VAR_GLO_IFB;
			var[nv]*=PPP_Glo.ecliF[sat-1];

			if (sys==SYS_CMP) {  //to reduce weight of the BDS GEO satellites
				prn=PPP_Glo.sFlag[sat-1].prn;
				if (prn>=1&&prn<=5)               {var[nv]*=100.0;}
				else if (prn>=6&&prn<=10)         {var[nv]*=1.0;}
				else if (prn>=11&&prn<=MAXPRNCMP) {var[nv]*=1.0;}
				else                              {var[nv]*=1.0;}
			}
			if (sys==SYS_QZS) {
				var[nv]*=25.0;
			}

			/* reject satellite by pre-fit residuals */
			if (!post&&opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
				sprintf(PPP_Glo.chMsg,"*** WARNING: outlier (%d) rejected %s sat=%s %s%d res=%9.4f el=%4.1f\n",
					post,str,PPP_Glo.sFlag[sat-1].id,j%2?"P":"L",j/2+1,v[nv],azel[1+i*2]*R2D);
				outDebug(OUTWIN,OUTFIL,OUTTIM);
				exc[i]=1; rtk->ssat[sat-1].rejc[j%2]++;
				continue;
			}
			/* record large post-fit residuals */
			if (post&&fabs(v[nv])>sqrt(var[nv])*THRES_REJECT) {
				obsi[ne]=i; frqi[ne]=j; ve[ne]=v[nv]; ne++;
			}

			if (j%2==0) rtk->ssat[sat-1].vsat[j/2]=1;

			if (!post) {
				if (j==0)      rtk->ssat[sat-1].resc_pri[0]=v[nv];
				else if (j==1) rtk->ssat[sat-1].resp_pri[0]=v[nv];
				else if (j==2) rtk->ssat[sat-1].resc_pri[1]=v[nv];
				else if (j==3) rtk->ssat[sat-1].resp_pri[1]=v[nv];
				else if (j==4) rtk->ssat[sat-1].resc_pri[2]=v[nv];
				else if (j==5) rtk->ssat[sat-1].resp_pri[2]=v[nv];
			}
			else {
				if (j==0)      {rtk->ssat[sat-1].resc_pos[0]=v[nv];}
				else if (j==1) {rtk->ssat[sat-1].resp_pos[0]=v[nv];}
				else if (j==2) {rtk->ssat[sat-1].resc_pos[1]=v[nv];}
				else if (j==3) {rtk->ssat[sat-1].resp_pos[1]=v[nv];}
				else if (j==4) {rtk->ssat[sat-1].resc_pos[2]=v[nv];}
				else if (j==5) {rtk->ssat[sat-1].resp_pos[2]=v[nv];}
			}

			nv++;
		}
		rtk->ssat[sat-1].azel[0]=azel[0+i*2];
		rtk->ssat[sat-1].azel[1]=azel[1+i*2];

		//constraint to external ionosphere correction
		if ((opt->ionoopt==IONOOPT_UC1||opt->ionoopt==IONOOPT_UC12)&&PPP_Glo.prcOpt_Ex.ion_const) {
			nx_nv=nx*nv;
			jj=II(sat,opt);
			v[nv]=dion_tec-rtk->x[jj];
			for (kk=0;kk<nx;kk++) H[kk+nx_nv]=kk==jj?1.0:0.0;
			var[nv++]=2.0;
		}
	}

	sys=rtk->opt.navsys;
	//add datum constraint for GLOICB_1SAT
	if (sys==SYS_GLO&&opt->gloicb==GLOICB_1SAT&&PPP_Glo.iEpoch==1) {
		nx_nv=nx*nv;
		for (kk=0;kk<nx;kk++) H[kk+nx_nv]=0.0;
		for (i=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			sat_r=sat-NSATGPS;
			ic=IICB(sat_r,opt);
			H[ic+nx_nv]=1.0;
			v[nv]-=x[ic];
		}
		var[nv]=1e-8;

		nv++;
	}
	//add datum constraint for GLOICB_1FRQ
	if (sys==SYS_GLO&&opt->gloicb==GLOICB_1FRQ&&PPP_Glo.iEpoch==1) {
		nx_nv=nx*nv;
		for (kk=0;kk<nx;kk++) H[kk+nx_nv]=0.0;
		for (i=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			frq=get_glo_fcn(sat,nav);
			for (ii=1;ii<14;ii++) {
				if (glo_fcn[ii-1]==frq) break;
			}
			ic=IICB(ii,opt);
			H[ic+nx_nv]=1.0;
			v[nv]-=x[ic];
		}
		var[nv]=1e-8;

		nv++;
	}

	/* reject satellite with large and max post-fit residual */
	if (post&&ne>0) {
		vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
		for (j=1;j<ne;j++) {
			if (fabs(vmax)>=fabs(ve[j])) continue;
			vmax=ve[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
		}
		sat=obs[maxobs].sat;
		sprintf(PPP_Glo.chMsg,"*** WARNING: outlier (%d) rejected %s sat=%s %s%d res=%9.4f el=%4.1f\n",
			post,str,PPP_Glo.sFlag[sat-1].id,maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
		outDebug(OUTWIN,OUTFIL,OUTTIM);
		exc[maxobs]=1; rtk->ssat[sat-1].rejc[maxfrq%2]++; stat=0;
		ve[rej]=0;
	}

	for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
		R[i+j*nv]=i==j?var[i]:0.0;
	}
	return post?stat:nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
	return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
	const prcopt_t *opt=&rtk->opt;
	int i,j;

	/* test # of valid satellites */
	rtk->sol.ns[0]=0;
	for (i=0;i<n&&i<MAXOBS;i++) {
		for (j=0;j<opt->nf;j++) {
			if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
			rtk->ssat[obs[i].sat-1].lock[j]++;
			rtk->ssat[obs[i].sat-1].outc[j]=0;
			if (j==0) rtk->sol.ns[0]++;
		}
	}
	rtk->sol.stat=rtk->sol.ns[0]<MIN_NSAT_SOL?SOLQ_NONE:stat;

	if (rtk->sol.stat==SOLQ_FIX) {
		for (i=0;i<3;i++) {
			rtk->sol.rr[i]=rtk->xa[i];
			rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
		}
		rtk->sol.qr[3]=(float)rtk->Pa[1];
		rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
		rtk->sol.qr[5]=(float)rtk->Pa[2];
	}
	else {
		for (i=0;i<3;i++) {
			rtk->sol.rr[i]=rtk->x[i];
			rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
		}
		rtk->sol.qr[3]=(float)rtk->P[1];
		rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
		rtk->sol.qr[5]=(float)rtk->P[2];
	}
	rtk->sol.dtr[0]=rtk->x[IC(0,opt)];
	rtk->sol.dtr[1]=rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)];
	rtk->sol.dcb_rcv=rtk->x[ID(opt)];

	for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
		rtk->ssat[obs[i].sat-1].snr[j]=obs[i].SNR[j];
	}
	for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
		if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
		if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
	}
}
/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const prcopt_t *opt=&rtk->opt;
	double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp;
	double rr[3],pos[3],dr[3]={0},enu[3]={0},dtdx[3]={0};
	double dtrp=0.0,shd=0.0,vart=0.0;
	double dantr[NFREQ]={0},dants[NFREQ]={0};
	double cosaz,sinaz,cosel,sinel;
	char str[32];
	int i,j,nv,sat,info,svh[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE;

	time2str(obs[0].time,str,2);

	rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);

	for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;

	/* satellite positions and clocks */
	satposs_rtklib(obs[0].time,obs,n,nav,opt->sateph,rs,dts,var,svh);

	/* exclude measurements of eclipsing satellite (block IIA) */
	testeclipse(obs,n,nav,rs);

	//calculate satellite elevation
	calElev(rtk,obs,n,rs);

	//for modeling ionosphere variations between epochs
	/* initialize position for first epoch */
	if (norm(rtk->x,3)<=0.0)
		for (i=0;i<3;i++) PPP_Glo.rr[i]=rtk->sol.rr[i];
	else
		for (i=0;i<3;i++) PPP_Glo.rr[i]=rtk->x[i];
	for (i=0;i<n;i++) {
		sat=obs[i].sat;
		for (j=0;j<3;j++) PPP_Glo.ssat_Ex[sat-1].rs[j]=rs[j+i*6];

		/* satellite and receiver antenna model */
		satantpcv(sat,rs+i*6,PPP_Glo.rr,nav->pcvs+sat-1,PPP_Glo.ssat_Ex[sat-1].dants);
		antmodel(sat,&opt->pcvr,opt->antdel,rtk->ssat[sat-1].azel,1,PPP_Glo.ssat_Ex[sat-1].dantr);

		/* phase windup model */
		if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,2,rs+i*6,PPP_Glo.rr,&PPP_Glo.ssat_Ex[sat-1].phw)) {
			continue;
		}
	}

	//to compute the elements of initialized PPP files
	if (PPP_Glo.outFp[OFILE_IPPP]) {
		//epoch time
		PPP_Info.t=obs[0].time;

		for (i=0;i<3;i++) rr[i]=PPP_Glo.crdTrue[i];
		if (norm(rr,3)==0.0) rr[i]=rtk->sol.rr[i];
		ecef2pos(rr,pos);

		//the displacements by earth tides
		tidedisp(gpst2utc(obs[0].time),rr,7,&nav->erp,opt->odisp[0],dr);

		for (i=0;i<MAXSAT;i++) PPP_Info.ssat[i].vs=0;
		for (i=0;i<n&&i<MAXOBS;i++) {
			sat=obs[i].sat;
			PPP_Info.ssat[sat-1].vs=1;

			//satellite position and clock offsets
			matcpy(PPP_Info.ssat[sat-1].satpos,rs+i*6,6,1);
			PPP_Info.ssat[sat-1].satclk=dts[i*2]*CLIGHT;

			//the eclipse satellites
			PPP_Info.ssat[sat-1].flag=0;
			if (norm(PPP_Info.ssat[sat-1].satpos,3)==0.0) PPP_Info.ssat[sat-1].flag=1;

			//sagnac effect
			PPP_Info.ssat[sat-1].dsag=sagnac(rs+i*6,rr);

			//satellite azimuth and elevation
			PPP_Info.ssat[sat-1].azel[0]=rtk->ssat[sat-1].azel[0];
			PPP_Info.ssat[sat-1].azel[1]=rtk->ssat[sat-1].azel[1];

			//line-of-sight (LOS) unit vector
			cosel=cos(PPP_Info.ssat[sat-1].azel[1]); sinel=sin(PPP_Info.ssat[sat-1].azel[1]);
			cosaz=cos(PPP_Info.ssat[sat-1].azel[0]); sinaz=sin(PPP_Info.ssat[sat-1].azel[0]);

			//convert 3D tidal displacements to LOS range
			ecef2enu(pos,dr,enu);
			PPP_Info.ssat[sat-1].dtid=enu[1]*cosel*cosaz+enu[0]*cosel*sinaz+enu[2]*sinel;

			//tropospheric zenith total delays (ZTDs)
			if (!model_trop(obs[i].time,pos,PPP_Info.ssat[sat-1].azel,opt,
				rtk->x,dtdx,nav,&dtrp,&shd,&vart)) continue;
			PPP_Info.ssat[sat-1].dtrp=dtrp;
			PPP_Info.ssat[sat-1].shd=shd;
			PPP_Info.ssat[sat-1].wmap=dtdx[0];

			//satellite and receiver antenna model
			satantpcv(sat,rs+i*6,rr,nav->pcvs+sat-1,dants);
			antmodel(sat,&opt->pcvr,opt->antdel,PPP_Info.ssat[sat-1].azel,1,dantr);
			for (j=0;j<NFREQ;j++) PPP_Info.ssat[sat-1].dant[j]=dants[j]+dantr[j];

			/* phase windup model */
			if (!model_phw(rtk->sol.time,sat,nav->pcvs[sat-1].type,2,rs+i*6,rr,
				&PPP_Info.ssat[sat-1].phw)) continue;

			for (j=0;j<NFREQ;j++) {
				PPP_Info.ssat[sat-1].L[j]=PPP_Info.ssat[sat-1].P[j]=0.0;
				PPP_Info.ssat[sat-1].L[j]=obs[i].L[j];
				PPP_Info.ssat[sat-1].P[j]=obs[i].P[j];
			}
		}
	}

	detecs(rtk,obs,n,nav);

	/* temporal update of ekf states */
	udstate_ppp(rtk,obs,n,nav);

	nv=n*rtk->opt.nf*2+MAXSAT+3;
	xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
	v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv);
	for (i=0;i<MAX_ITER;i++) {
		matcpy(xp,rtk->x,rtk->nx,1);
		matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

		/* prefit residuals */
		if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,exc,nav,xp,rtk,v,H,R,azel))) {
			sprintf(PPP_Glo.chMsg,"*** WARNING: %s ppp (%d) no valid obs data\n",str,i+1);
			outDebug(OUTWIN,OUTFIL,OUTTIM);
			break;
		}
		/* measurement update of ekf states */
		if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {
			sprintf(PPP_Glo.chMsg,"*** ERROR: %s ppp (%d) filter error info=%d\n",str,i+1,info);
			outDebug(OUTWIN,OUTFIL,OUTTIM);
			break;
		}
		/* postfit residuals */
		if (ppp_res(i+1,obs,n,rs,dts,var,svh,exc,nav,xp,rtk,v,H,R,azel)) {
			matcpy(rtk->x,xp,rtk->nx,1);
			matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
			stat=SOLQ_PPP;
			break;
		}
	}
	if (i>=MAX_ITER) {
		sprintf(PPP_Glo.chMsg,"*** WARNING: %s ppp (%d) iteration overflows\n",str,i);
		outDebug(OUTWIN,OUTFIL,OUTTIM);
	}
	if (stat==SOLQ_PPP) {
		/* ambiguity resolution in ppp */
		/*if (ppp_ar(rtk,obs,n,exc,nav,azel,xp,Pp)&&
		ppp_res(9,obs,n,rs,dts,var,svh,dr,exc,nav,xp,rtk,v,H,R,azel)) {

		matcpy(rtk->xa,xp,rtk->nx,1);
		matcpy(rtk->Pa,Pp,rtk->nx,rtk->nx);

		for (i=0;i<3;i++) std[i]=sqrt(Pp[i+i*rtk->nx]);
		if (norm(std,3)<MAX_STD_FIX) stat=SOLQ_FIX;
		}
		else {*/
		rtk->nfix=0;
		//}
		/* update solution status */
		update_stat(rtk,obs,n,stat);
	}
	free(rs); free(dts); free(var); free(azel);
	free(xp); free(Pp); free(v); free(H); free(R);
}
