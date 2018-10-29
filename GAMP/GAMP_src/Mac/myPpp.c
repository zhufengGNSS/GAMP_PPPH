#include "gamp.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))

//calculate threshold values for cycle slip detection
extern int calCsThres(prcopt_t *opt, const double sample)
{
	int b=0;

	if (sample>0.0) {
		if (PPP_Glo.prcOpt_Ex.bUsed_gfCs==1&&fabs(PPP_Glo.prcOpt_Ex.csThresGF)<0.01) {
			if (sample<=1.0)        PPP_Glo.prcOpt_Ex.csThresGF=0.05;
			else if (sample<=20.0)  PPP_Glo.prcOpt_Ex.csThresGF=(0.10)/(20.0-0.0)*sample+0.05;
			else if (sample<=60.0)  PPP_Glo.prcOpt_Ex.csThresGF=0.15;
			else if (sample<=100.0) PPP_Glo.prcOpt_Ex.csThresGF=0.25;
			else                    PPP_Glo.prcOpt_Ex.csThresGF=0.35;

			b=1;
		}
		if (PPP_Glo.prcOpt_Ex.bUsed_mwCs==1&&fabs(PPP_Glo.prcOpt_Ex.csThresMW)<0.01) {
			if (sample<=1.0)        PPP_Glo.prcOpt_Ex.csThresMW=2.5;
			else if (sample<=20.0)  PPP_Glo.prcOpt_Ex.csThresMW=(2.5)/(20.0-0.0)*sample+2.5;
			else if (sample<=60.0)  PPP_Glo.prcOpt_Ex.csThresMW=5.0;
			else                    PPP_Glo.prcOpt_Ex.csThresMW=7.5;

			b=1;
		}

		return b;
	}
	else {
		//sample<=0.0
		PPP_Glo.prcOpt_Ex.csThresGF=0.15;
		PPP_Glo.prcOpt_Ex.csThresMW=5.0;
		b=0;
	}

	return b;
}

//pseudorange observation checking
extern void obsScan_SPP(const prcopt_t *popt, obsd_t *obs, const int nobs, int *nValid)
{
	double dt;
	int i,j,n,sat,sys;

	for (i=n=0;i<nobs;i++) {
		sat=obs[i].sat;
		sys=PPP_Glo.sFlag[sat-1].sys;

		/* exclude satellites */
		if (!(sys&popt->navsys)) continue;
		if (popt->exsats[sat-1])	 continue;

		dt=0.0;
		for (j=0;j<NFREQ;j++) {
			dt+=obs[i].P[j]*obs[i].P[j];
		}
		if (dt==0.0)	 continue;

		obs[n++]=obs[i];
	}

	if (nValid) *nValid=n;
}
extern void obsScan_PPP(const prcopt_t *popt, obsd_t *obs, const int nobs, int *nValid)
{
	int i,n,sat,f2;

	for (i=n=0;i<nobs&&i<MAXOBS;i++) {
		sat=obs[i].sat;

		f2=1;
		//if (NFREQ>=3&&(PPP_Glo.sFlag[sat-1].sys&(SYS_GAL|SYS_SBS))) f2=2;
		if (popt->mode>=PMODE_PPP_KINEMA) {
			if (obs[i].L[0]*obs[i].L[f2]==0.0) continue;
		}

		if (fabs(obs[i].P[0]-obs[i].P[f2])>=200.0) continue;

		obs[n]=obs[i];
		n++;
	}

	if (nValid) *nValid=n;
}

/* mutipath correct-------------------------------------------------------------
* BeiDou satellite-induced code pseudorange variations correct
* args   : rtk_t *rtk       IO  rtk control/result struct
           obsd_t *obs      IO  observation data
           int    n         I   number of observation data
		   nav_t  *nav      I   navigation messages
* note   : 
* 
* -----------------------------------------------------------------------------*/
extern void BDmultipathCorr(rtk_t *rtk, obsd_t *obs, int n)
{
	int i,j,sat,prn,b;
	double dmp[3],elev,a;
	const static double IGSOCOEF[3][10]={		/* m */
		{-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},	//B1
		{-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},	//B2
		{-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32},	//B3
	};
	const static double MEOCOEF[3][10]={		/* m */
		{-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},	//B1
		{-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},	//B2
		{-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47},	//B3
	};

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
        if (PPP_Glo.sFlag[sat-1].sys!=SYS_CMP) continue;

		prn=PPP_Glo.sFlag[sat-1].prn;
		if (prn<=5) continue;

		elev=rtk->ssat[sat-1].azel[1]*R2D;

		if (elev<=0.0) continue;

		for (j=0;j<3;j++) dmp[j]=0.0;

		a=elev*0.1;
		b=(int)a;

		if (prn>=6&&prn<11) { // IGSO(C06, C07, C08, C09, C10)
			if (b<0) {
				for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][0];
			}
			else if (b>=9) {
				for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][9];
			}
			else {
				for (j=0;j<3;j++) dmp[j]=IGSOCOEF[j][b]*(1.0-a+b)+IGSOCOEF[j][b+1]*(a-b);
			}
		}
		else if (prn>=11) {   // MEO(C11, C12, C13, C14)
			if (b<0) {
				for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][0];
			}
			else if (b>=9) {
				for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][9];
			}
			else {
				for (j=0;j<3;j++) dmp[j]=MEOCOEF[j][b]*(1.0-a+b)+MEOCOEF[j][b+1]*(a-b);
			}
		}

		for (j=0;j<3;j++) obs[i].P[j]+=dmp[j];
	}
}

/* L1/L2 wide-lane phase measurement -----------------------------------------*/
extern double wlAmbMeas(const obsd_t *obs, const nav_t *nav)
{
	int i=0,j=1;
	const double *lam=nav->lam[obs->sat-1];
	double P1,P2,P1_C1,P2_C2,lam1,lam2,res;

	if (obs->L[i]==0.0) return 0.0;
	if (obs->L[j]==0.0) return 0.0;
	if (obs->P[i]==0.0) return 0.0;
	if (obs->P[j]==0.0) return 0.0;
	if (lam[i]*lam[j]==0.0) return 0.0;

	P1=obs->P[i];
	P2=obs->P[j];
	P1_C1=nav->cbias[obs->sat-1][1];
	P2_C2=nav->cbias[obs->sat-1][2];

	if (obs->code[0]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
	if (obs->code[1]==CODE_L2C) P2+=P2_C2; /* C2->P2 */

	lam1=lam[i];
	lam2=lam[j];
	res=(obs->L[i]-obs->L[j])-(lam2-lam1)/(lam1+lam2)*(P1/lam1+P2/lam2);

	return res;
}

/* detect cycle slip by widelane jump -----------------------------*/
static void detslp_mw(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	int i,j,nd=0,sat;
	int bSlip[MAXSAT],bLowElev=0;
	double wl0,wl1,elev,el,thres,thres0,delta[50],dmw[MAXSAT],dtmp,fact;
	double std_ex,ave_ex;
	const double rad_20=20.0*D2R;

	for (i=0;i<MAXSAT;i++) {
		dmw[i]=0.0;
		bSlip[i]=0;
	}

	//the thresh values is suitable for 30s interval
	fact=1.0;
	if (PPP_Glo.sample>=29.5) {
		if (PPP_Glo.delEp<=2) fact=1.0;
		else if (PPP_Glo.delEp<=4) fact=1.25;
		else if (PPP_Glo.delEp<=6) fact=1.5;
		else	 fact=2.0;
	}

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		if (timediff(PPP_Glo.tNow,PPP_Glo.ssat_Ex[sat-1].tLast)>PPP_Glo.sample) {
			PPP_Glo.ssat_Ex[sat-1].mw[1]=0.0;
			PPP_Glo.ssat_Ex[sat-1].mwIndex=0;
		}

		wl1=wlAmbMeas(obs+i,nav);
		wl0=PPP_Glo.ssat_Ex[sat-1].mw[1];

		/*sprintf(PPP_Glo.chMsg,"sat=%s, mw0=%13.3f, mw1=%13.3f\n",PPP_Glo.sFlag[sat-1].id,wl1,
			PPP_Glo.ssat_Ex[sat-1].mw[1]);
		outDebug(1,1,1);*/

		if (wl1==0.0||wl0==0.0) continue;

		elev=rtk->ssat[sat-1].azel[1];
		el=elev;

		bLowElev=0;

		if (elev<rtk->opt.elmin) {
			el=rtk->opt.elmin;
			bLowElev=1;
		}

		dtmp=el*R2D;

		if (el>=rad_20) thres=PPP_Glo.prcOpt_Ex.csThresMW;
		else	 thres=-PPP_Glo.prcOpt_Ex.csThresMW*0.1*dtmp+3*PPP_Glo.prcOpt_Ex.csThresMW;

		dmw[sat-1]=wl1-wl0;

		thres0=MIN(thres*fact,6.0);
		if (fabs(wl1-wl0)>MIN(thres*fact,6.0)) {
			bSlip[sat-1]=1;
			delta[nd++]=wl1-wl0;

			if (bLowElev) continue;

			sprintf(PPP_Glo.chMsg, "%s MW CS mw_new=%13.3f, mw_old=%13.3f, diff_abs=%13.3f, thres=%13.3f, elev=%4.1f\n", 
				PPP_Glo.sFlag[sat-1].id,wl1,wl0,wl1-wl0,thres0,elev);
			outDebug(OUTWIN,OUTFIL,OUTTIM);
		}
	}

	if (2*nd+1<=n&&nd<n-3&&nd<=3) {

	}
	else if (nd>2) {
		i=findGross(0,0,delta,nd,0,&std_ex,&ave_ex,NULL,4.0,0.3,0.2);

		if (i<=1&&std_ex<=1.0&&fabs(ave_ex)<=10.0 ) {
			sprintf(PPP_Glo.chMsg,"*** WARNING: MW CS abnormally at this epoch %4.1f  %4.2f %d %d\n",ave_ex,std_ex,nd,n);
			outDebug(OUTWIN,OUTFIL,OUTTIM);

			for (i=0;i<n&&i<MAXOBS;i++) {
				sat=obs[i].sat;
				dmw[sat-1]=0.0;
				bSlip[sat-1]=0;
			}
		}
	}

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;

		if (!bSlip[sat-1]) {
			if (fabs(dmw[sat-1])>=1.0) 
				PPP_Glo.ssat_Ex[sat-1].obsStat.iCycle=1;
		}
		else {
			for (j=0;j<rtk->opt.nf;j++) 
				rtk->ssat[sat-1].slip[j]|=1;

			PPP_Glo.ssat_Ex[sat-1].obsStat.iCycle=2;
		}
	}
}

/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	double g0,g1,elev,thres=0.0,deltaEpoch=1.0;
	double elev0,thres0;
	int i,j,sat;

	deltaEpoch=fabs(rtk->tt/PPP_Glo.sample);

	for (i=0;i<n&&i<MAXOBS;i++) 
	{
		sat=obs[i].sat;
		elev=rtk->ssat[sat-1].azel[1]*R2D;
		elev0=elev;

		g1=gfmeas(obs+i,nav);

		if (g1==0.0) continue;

		if (elev<rtk->opt.elmin*R2D) elev=rtk->opt.elmin*R2D;


		if (elev>=15.0) thres=PPP_Glo.prcOpt_Ex.csThresGF;
		else thres=-PPP_Glo.prcOpt_Ex.csThresGF/15.0*elev+2*PPP_Glo.prcOpt_Ex.csThresGF;

		g0=rtk->ssat[sat-1].gf;

		thres0=MIN(thres*deltaEpoch,1.5);
		if (g0!=0.0&&fabs(g1-g0)>MIN(thres*deltaEpoch,1.5)) 
		{
			for (j=0;j<rtk->opt.nf;j++) 
				rtk->ssat[sat-1].slip[j]|=1;

			sprintf(PPP_Glo.chMsg,"%s GF CS gf_new=%13.3f, gf_old=%13.3f, diff_abs=%13.3f, thres=%13.3f, elev=%4.1f\n", 
				PPP_Glo.sFlag[sat-1].id,g1,g0,g1-g0,thres0,elev0);
			outDebug(OUTWIN,OUTFIL,OUTTIM);
		}
	}
}

extern void detecs(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	int i,j,b1,b2;
	double dt;

	for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
		rtk->ssat[i].slip[j]=0;
	}

	dt=1.0;
	if (PPP_Glo.sample>0.0) dt=fabs(rtk->tt/PPP_Glo.sample);
	PPP_Glo.delEp=myRound(dt);

	if (PPP_Glo.sample<=1.5) {
		if (fabs(rtk->tt)<=10.0)	 PPP_Glo.delEp=1;
		else if (fabs(rtk->tt)<=15.0) PPP_Glo.delEp=2;
		else if (fabs(rtk->tt)<=22.0) PPP_Glo.delEp=3;
	}

	if (PPP_Glo.delEp<=0) {
		b1=(PPP_Glo.iEpoch==1)&&(PPP_Glo.revs==0) ;
		b2=(PPP_Glo.nEpoch-1==PPP_Glo.iEpoch)&&(PPP_Glo.revs==1) ;

		if (!b1&&!b2) {
			sprintf(PPP_Glo.chMsg,"deltaEpoch=%d\n",PPP_Glo.delEp);
			outDebug(OUTWIN,OUTFIL,OUTTIM);
		}
		PPP_Glo.delEp=1;
	}

	/* detect cycle slip by LLI */
	//detslp_ll(rtk,obs,n);

	/* detect slip by Melbourne-Wubbena linear combination jump */
	if (PPP_Glo.prcOpt_Ex.bUsed_mwCs) detslp_mw(rtk,obs,n,nav);

	/* detect cycle slip by geometry-free phase jump */
	if (PPP_Glo.prcOpt_Ex.bUsed_gfCs) detslp_gf(rtk,obs,n,nav);
}

//save the information of current epoch
extern void keepEpInfo(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	int i,j,sat;
	prcopt_t *opt=&rtk->opt;
	double wl0,wl1,var0,var1,gf;

	for (i=0;i<MAXSAT;i++) {
		rtk->ssat[i].gf=0.0;
		PPP_Glo.ssat_Ex[i].mw[0]=0.0;
	}

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;

		PPP_Glo.ssat_Ex[sat-1].tLast=PPP_Glo.tNow;

		//if ( rtk->ssat[sat-1].azel[1]<rtk->opt.elmin ) 
		//	continue;

		if ((gf=gfmeas(obs+i,nav))!=0.0)
			rtk->ssat[sat-1].gf=gf;

		if ((wl1=wlAmbMeas(obs+i,nav))==0.0)
			continue;

		wl0=PPP_Glo.ssat_Ex[sat-1].mw[1];
		PPP_Glo.ssat_Ex[sat-1].mw[0]=wl1;

		if (PPP_Glo.ssat_Ex[sat-1].mwIndex>0) {
			j=PPP_Glo.ssat_Ex[sat-1].mwIndex;
			var0=PPP_Glo.ssat_Ex[sat-1].mwVar_c;
			var1=(wl1-wl0)*(wl1-wl0)-var0;
			var1=var0 + var1/j;

			PPP_Glo.ssat_Ex[sat-1].mw[1]=(wl0*j+wl1)/(j+1);
			PPP_Glo.ssat_Ex[sat-1].mwIndex++;
			PPP_Glo.ssat_Ex[sat-1].mwVar_c=var1;
		}
		else {
			PPP_Glo.ssat_Ex[sat-1].mw[1]=wl1;
			PPP_Glo.ssat_Ex[sat-1].mwIndex++;
			PPP_Glo.ssat_Ex[sat-1].mwVar_c=0.25;
		}

		j=IB(sat,0,opt);

		PPP_Glo.ssat_Ex[sat-1].arc.ifArc_m=rtk->x[j];
		PPP_Glo.ssat_Ex[sat-1].arc.ifVarArc_m=rtk->P[j*rtk->nx+j];

		PPP_Glo.ssat_Ex[sat-1].arc.mwArc_c=PPP_Glo.ssat_Ex[sat-1].mw[1];
		PPP_Glo.ssat_Ex[sat-1].arc.mwArcVar_c=PPP_Glo.ssat_Ex[sat-1].mwVar_c;
	}
}

//calculate satellite elevation
extern void calElev(rtk_t *rtk, const obsd_t *obs, int n, double *rs)
{
	int i,sat;
	double rr[3]={0},pos[3],r,e[3];

	for (i=0;i<MAXSAT;i++)
		rtk->ssat[i].azel[1]=0.0;

	if ( 0.0!=PPP_Glo.crdTrue[0] ) {
		for (i=0;i<3;i++)
			rr[i]=PPP_Glo.crdTrue[i];
	}
	else {
		for (i=0;i<3;i++) 
			rr[i]=rtk->x[i];

		if (rr[0]==0.0) {
			for (i=0;i<3;i++) 
				rr[i]=rtk->sol.rr[i];
		}
	}

	if (norm(rr,3)<=100.0) return ;

	ecef2pos(rr,pos);

	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;

		if (!PPP_Glo.sFlag[sat-1].sys) continue;

		/* geometric distance/azimuth/elevation angle */
		if ((r=geodist(rs+i*6,rr,e))<0) continue;
		satazel(pos,e,rtk->ssat[sat-1].azel);
	}
}

#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
//From GLAB
/*****************************************************************************
* Name        : gravitationalDelayCorrection
* Description : Obtains the gravitational delay correction for the effect of 
*               general relativity (red shift) to the GPS signal
* Parameters  :
* Name                           |Da|Unit|Description
* double  *receiverPosition       I  m    Position of the receiver
* double  *satellitePosition      I  m    Position of the satellite
* Returned value (double)         O  m    Gravitational delay correction
*****************************************************************************/
extern double gravitationalDelayCorrection (const int sys, const double *receiverPosition, 
	                                        const double *satellitePosition) 
{
	double	receiverModule;
	double	satelliteModule;
	double	distance;
	double  MU=MU_GPS;

	receiverModule=sqrt(receiverPosition[0]*receiverPosition[0]+receiverPosition[1]*receiverPosition[1]+
		receiverPosition[2]*receiverPosition[2]);
	satelliteModule=sqrt(satellitePosition[0]*satellitePosition[0]+satellitePosition[1]*satellitePosition[1]+
		satellitePosition[2]*satellitePosition[2]);
	distance=sqrt((satellitePosition[0]-receiverPosition[0])*(satellitePosition[0]-receiverPosition[0])+
		(satellitePosition[1]-receiverPosition[1])*(satellitePosition[1]-receiverPosition[1])+
		(satellitePosition[2]-receiverPosition[2])*(satellitePosition[2]-receiverPosition[2]));

	switch (sys) {
	case SYS_GPS:
		MU=MU_GPS;
		break;
	case SYS_GLO:
		MU=MU_GLO;
		break;
	case SYS_GAL:
		MU=MU_GAL;
		break;
	case SYS_CMP:
		MU=MU_CMP;
		break;
	default:
		MU=MU_GPS;
		break;
	}

	return 2.0*MU/(CLIGHT*CLIGHT)*log((satelliteModule+receiverModule+distance)/(satelliteModule+receiverModule-distance));
}

