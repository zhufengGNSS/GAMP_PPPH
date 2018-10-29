/*------------------------------------------------------------------------------
* gampPos.c : post-processing positioning
*-----------------------------------------------------------------------------*/
#include "gamp.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))

#define MAXPRCDAYS  100          /* max days of continuous processing */

/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss={0};        /* receiver antenna parameters */
static pcvs_t pcvsr={0};        /* satellite antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sta_t stas[MAXRCV];      /* station infomation */


/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* input obs data, navigation messages and sbas correction -------------------*/
int inputobs(obsd_t *obs, obs_t obss, int revs, int *iobsu, int *iepoch)
{
	int i,nu,n=0;

	if (!revs) {    /* input forward data */
		if ((nu=nextobsf(&obss,iobsu,1))<=0) 
			return -1;
		for (i=0;i<nu&&n<MAXOBS;i++) 
			obs[n++]=obss.data[*iobsu+i];
		*iobsu+=nu;
		*iepoch+=1;
	}
	else {        /* input backward data */
		if ((nu=nextobsb(&obss,iobsu,1))<=0) 
			return -1;
		for (i=0;i<nu&&n<MAXOBS;i++)
			obs[n++]=obss.data[*iobsu-nu+1+i];
		*iobsu-=nu;
		*iepoch-=1;
	}
	return n;
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
static void readpreceph(char **infile, int n, const prcopt_t *prcopt, nav_t *nav)
{
    int i;
    
    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    
    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */
    for (i=0;i<n;i++) {
		if (strlen(infile[i])<=3) continue;
        readrnxc(infile[i],nav);
    }
}
/* free prec ephemeris -----------------------------------------*/
static void freepreceph(nav_t *nav)
{
    int i;
    
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec); nav->tec =NULL; nav->nt=nav->ntmax=0;
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char *infile[MAXINFILE],
	                  const int *index, int n, const prcopt_t *prcopt,
	                  obs_t *obs, nav_t *nav, sta_t *sta)
{
	int i,j,ind=0,nobs=0,rcv=1,nep;

	obs->data=NULL; obs->n =obs->nmax =0;
	nav->eph =NULL; nav->n =nav->nmax =0;
	nav->geph=NULL; nav->ng=nav->ngmax=0;
	PPP_Glo.nEpoch=0;

	for (i=0;i<n;i++) {
		if (index[i]!=ind) {
			if (obs->n>nobs) rcv++;
			ind=index[i]; nobs=obs->n; 
		}
		/* read rinex obs and nav file */
		nep=readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt,obs,nav,rcv<=2?sta+rcv-1:NULL);
	}
	if (obs->n<=0) {
		printf("*** ERROR: no obs data!\n");
		return 0;
	}
	if (nav->n<=0&&nav->ng<=0) {
		printf("*** ERROR: no nav data!\n");
		return 0;
	}

	/* sort observation data */
	PPP_Glo.nEpoch=sortobs(obs);

	/* delete duplicated ephemeris */
	uniqnav(nav);

	/* set time span for progress display */
	if (ts.time==0||te.time==0) {
		for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
		for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
		if (i<j) {
			if (ts.time==0) ts=obs->data[i].time;
			if (te.time==0) te=obs->data[j].time;
			settspan(ts,te);
		}
	}

	if (prcopt->navsys&SYS_GLO) {
		if (nav->ng<=0) {
			printf("*** ERROR: nav->ng<=0!\n");
		}
	}

	return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{   
	if (obs->data) {free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;}
	if (nav->eph ) {free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;}
	if (nav->geph) {free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;}
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
	const pcvs_t *pcvr, const sta_t *sta)
{
	pcv_t *pcv;
	double pos[3],del[3],dt;
	int i,j,k=0,sys;
	char id[64];

	/* set satellite antenna parameters */
	for (i=0;i<MAXSAT;i++) {
		sys=PPP_Glo.sFlag[i].sys;
		if (!(sys&popt->navsys)) continue;
		if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
			satno2id(i+1,id);
			continue;
		}
		nav->pcvs[i]=*pcv;

		if (pcv->dzen==0.0) j=10;
		else j=myRound((pcv->zen2-pcv->zen1)/pcv->dzen);

		if (sys==SYS_GPS) k=0;
		else if (sys==SYS_GLO) k=0+1*NFREQ;
		else if (sys==SYS_CMP) k=0+2*NFREQ;
		else if (sys==SYS_GAL) k=0+3*NFREQ;
		//double dt=norm(pcv->var[0],j);
		dt=norm(pcv->var[k],j);
		if (dt<=0.0001) {
			sprintf(PPP_Glo.chMsg,"%s ATTENTION! PRELIMINARY PHASE CENTER CORRECTIONS!\n",
				PPP_Glo.sFlag[pcv->sat-1].id);
			outDebug(OUTWIN,OUTFIL,0);
		}
	}
	for (i=0;i<1;i++) {
		if (!strcmp(popt->anttype,"*")) { /* set by station parameters */
			strcpy(popt->anttype,sta[i].antdes);
			if (sta[i].deltype==1) { /* xyz */
				if (norm(sta[i].pos,3)>0.0) {
					ecef2pos(sta[i].pos,pos);
					ecef2enu(pos,sta[i].del,del);
					for (j=0;j<3;j++) popt->antdel[j]=del[j];
				}
			}
			else { /* enu */
				for (j=0;j<3;j++) popt->antdel[j]=stas[i].del[j];
			}
		}
		if (!(pcv=searchpcv(0,popt->anttype,time,pcvr))) {
			*popt->anttype='\0';
			continue;
		}
		strcpy(popt->anttype,pcv->type);
		popt->pcvr=*pcv;
	}
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    
    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}
/* write header to output file -----------------------------------------------*/
static int outhead(char **outfile, const prcopt_t *popt, const solopt_t *sopt, FILE* fp[], int sum)
{
	int i;

	for (i=0; i<sum; i++) {
		if (!outfile[i]||strlen(outfile[i])<=4) continue;

		createdir((const char *)outfile[i]);

		if (!(fp[i]=fopen(outfile[i],"w"))) {
			printf("error : open output file %s",outfile[i]);
			return 0;
		}

		switch (i) {
		case 0:
			//outsolhead(fp[0],sopt);
			break;
		case 1:
			//outCsInfoHead(fp[1],sopt);
			break;
		case 2:
			//outResiHead(fp[2],popt,sopt);
			break;
		case 3:
			//outResiHead(fp[3],popt,sopt);
			break;
		case 4:
			//outResiHead(fp[4],popt,sopt);
			break;
		case 5:
			//outAllAmbHead(fp[5],sopt,12);
			break;
		case 6:
			//outAllAmbHead(fp[6],sopt,7);
			break;
		case 7:
			//outAllAmbHead(fp[7],sopt,7);
			break;
		case 8:
			//outAllAmbHead(fp[8],sopt,7);
			break;
		case 9:
			//outAllAmbHead(fp[9],sopt,7);
			break;
		//case 10:
		//    outallambhead(fp[10],sopt,7);
		//    break;
		//case 11:
		//    outallambhead(fp[11],sopt,12);
		//    break;
		}
		fclose(fp[i]);
	}

	return 1;
}
/* write header to initialized files for PPP processing -----------------------------------------------*/
static void outiPppHead(FILE *fp, rtk_t rtk)
{
	unsigned char buff[MAXSOLMSG+1];
	char *p=(char *)buff;
	int i,j,n;
	int glo_frq[24]={0};
	int ig=0,ir=0,ic=0,ie=0,lg=0,lc=0,le=0;
	gtime_t tts,tte;
	double t1,t2;
	int w1,w2;
	char s8[32],s9[32];

	for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
	for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;
	tts=obss.data[i].time;
	tte=obss.data[j].time;
	t1=time2gpst(tts,&w1);
	t2=time2gpst(tte,&w2);
	time2str(tts,s8,1);
	time2str(tte,s9,1);

	/* get GLONASS frequency number */
	for (i=33;i<57;i++) {
		for (j=0;j<navs.ng;j++) {
			if (navs.geph[j].sat!=i) continue;
			glo_frq[i-33]=navs.geph[j].frq;
		}
	}

	for (i=0;i<obss.n;i++) {
		if (satsys(obss.data[i].sat,NULL)==SYS_GPS) {
			if (obss.data[i].P[0]!=0.0&&obss.data[i].P[1]!=0.0&&obss.data[i].P[2]!=0.0) {
				lg=1;
				ig=i;
				break;
			}
			else if (obss.data[i].P[0]!=0.0 && obss.data[i].P[1]!=0.0) {
				lg=0;
				ig=i;
			}
		}
	}
	for (i=0;i<obss.n;i++) {
		if (satsys(obss.data[i].sat,NULL)==SYS_GLO) {
			if (obss.data[i].P[1]!=0.0) {
				ir=i;
				break;
			}
		}
	}
	for (i=0;i<obss.n;i++) {
		if (satsys(obss.data[i].sat,NULL)==SYS_CMP) {
			if (obss.data[i].P[0]!=0.0&&obss.data[i].P[1]!=0.0&&obss.data[i].P[2]!=0.0) {
				lc=1;
				ic=i;
				break;
			}
			else if (obss.data[i].P[0]!=0.0&&obss.data[i].P[1]!=0.0) {
				lc=0;
				ic=i;
			}
		}
	}
	for (i=0;i<obss.n;i++) {
		if (satsys(obss.data[i].sat,NULL)==SYS_GAL) {
			if (obss.data[i].P[0]!=0.0&&obss.data[i].P[1]!=0.0&&obss.data[i].P[2]!=0.0) {
				le=1;
				ie=i;
				break;
			}
			else if (obss.data[i].P[0]!=0.0&&obss.data[i].P[1]!=0.0) {
				le=0;
				ie=i;
			}
		}
	}

	p+=sprintf(p,"+PPP_HEADER\n");
	p+=sprintf(p,"       STA_NAME: %s\n",stas[0].name);
	p+=sprintf(p,"       RCV_TYPE: %s\n",stas[0].rectype);
	p+=sprintf(p,"       ANT_TYPE: %s\n",stas[0].antdes);
	p+=sprintf(p,"        STA_POS: %14.4f%14.4f%14.4f (m)\n",
		PPP_Glo.crdTrue[0],PPP_Glo.crdTrue[1],PPP_Glo.crdTrue[2]);
	p+=sprintf(p,"       INTERVAL: %6.2f (s)\n",PPP_Glo.sample);
	if (lg) p+=sprintf(p,"       GPS_TYPE: C%2s  L%2s  C%2s  L%2s  C%2s  L%2s\n",
		obss.data[ig].type[0],obss.data[ig].type[0],obss.data[ig].type[1],
		obss.data[ig].type[1],obss.data[ig].type[2],obss.data[ig].type[2]);
	else p+=sprintf(p,"       GPS_TYPE: C%2s  L%2s  C%2s  L%2s\n",obss.data[ig].type[0],
		obss.data[ig].type[0],obss.data[ig].type[1],obss.data[ig].type[1]);
	p+=sprintf(p,"   GLONASS_TYPE: C%2s  L%2s  C%2s  L%2s\n",obss.data[ir].type[0],
		obss.data[ir].type[0],obss.data[ir].type[1],obss.data[ir].type[1]);
	if (lc) p+=sprintf(p,"    BEIDOU_TYPE: C%2s  L%2s  C%2s  L%2s  C%2s  L%2s\n",
		obss.data[ic].type[0],obss.data[ic].type[0],obss.data[ic].type[1],
		obss.data[ic].type[1],obss.data[ic].type[2],obss.data[ic].type[2]);
	else p+=sprintf(p,"    BEIDOU_TYPE: C%2s  L%2s  C%2s  L%2s\n",obss.data[ic].type[0],
		obss.data[ic].type[0],obss.data[ic].type[1],obss.data[ic].type[1]);
	if (le) p+=sprintf(p,"   GALILEO_TYPE: C%2s  L%2s  C%2s  L%2s  C%2s  L%2s\n",
		obss.data[ie].type[0],obss.data[ie].type[0],obss.data[ie].type[1],
		obss.data[ie].type[1],obss.data[ie].type[2],obss.data[ie].type[2]);
	else p+=sprintf(p,"   GALILEO_TYPE: C%2s  L%2s  C%2s  L%2s\n",obss.data[ie].type[0],
		obss.data[ie].type[0],obss.data[ie].type[1],obss.data[ie].type[1]);
	p+=sprintf(p,"     TROP MODEL: %s\n","GPT");
	p+=sprintf(p,"       TROP MAP: %s\n","GMF");
	p+=sprintf(p,"      ZTD MODEL: %s\n","SAAS");
	if (rtk.opt.sateph==EPHOPT_BRDC) {
		p+=sprintf(p,"     ORBIT TYPE: %s\n","Broadcast");
		p+=sprintf(p,"     CLOCK TYPE: %s\n","Broadcast");
	}
	else if (rtk.opt.sateph==EPHOPT_PREC) {
		p+=sprintf(p,"     ORBIT TYPE: %s\n","Precise");
		p+=sprintf(p,"     CLOCK TYPE: %s\n","Precise");
	}
	p+=sprintf(p,"     ELEV. MASK: %4.1f (deg)\n",0.0);
	p+=sprintf(p,"      GLONASS K: ");
	for (i=0;i<24;i++) p+=sprintf(p,"%4d",glo_frq[i]);
	p+=sprintf(p,"\n");
	p+=sprintf(p,"     START TIME: %s GPST (week %04d %8.1fs)\n",s8,w1,t1);
	p+=sprintf(p,"       END TIME: %s GPST (week %04d %8.1fs)\n",s9,w2,t2);
	p+=sprintf(p,"-PPP_HEADER\n");
	p+=sprintf(p,"\n");
	p+=sprintf(p,"+PPP_BLOCK\n");
	p+=sprintf(p,"*PRN,s:     SatPosX(m),     SatPosY(m),     SatPosZ(m),    SatClock(m),"
		" Elevation(deg),   Azimuth(deg),          P1(m),          P2(m),          P3(m),"
		"     L1(cycles),     L2(cycles),     L3(cycles),  Trop Delay(m),    Trop Map(m),"
		"      Sagnac(m), Tide Effect(m),      PCV_L1(m),      PCV_L2(m), Windup(cycles)");
	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    return !*outfile?stdout:fopen(outfile,"a");
}
//calculate sampling interval
static double sampledetermine(const prcopt_t *popt)
{
	obsd_t obs[MAXOBS];
	int i,j,nobs,n,m,it[MINNUM];
	gtime_t gt[MINNUM+1];
	double dt[MINNUM],dtmp;

	j=0;
	PPP_Glo.iObsu=PPP_Glo.revs=PPP_Glo.iEpoch=0;

	//MINNUM=30
	while ((nobs=inputobs(obs,obss,PPP_Glo.revs,&PPP_Glo.iObsu,&PPP_Glo.iEpoch))>=0) {
		gt[j++]=obs[0].time;

		if (j>MINNUM) break;
	}

	PPP_Glo.iObsu=PPP_Glo.revs=PPP_Glo.iEpoch=0;

	if (j<=MINNUM) {
		sprintf(PPP_Glo.chMsg,"*** WARNING: the number of epochs is less than %d.\n",MINNUM);
		outDebug(OUTWIN,OUTFIL,0);

		return 30.0;
	}

	for (i=0;i<MINNUM;i++) {
		dt[i]=0.0;
		it[i]=0;
	}

	dt[0]=timediff(gt[1],gt[0]);
	it[0]=1;
	n=1;
	for (i=0;i<MINNUM;i++) {
		dtmp=timediff(gt[i+1],gt[i]);

		for (j=0;j<n;j++) {
			if (fabs(dtmp-dt[j])<1.0e-8) {
				it[j]++;
				break;
			}
		}
		if (j==n) {
			dt[n]=dtmp;
			it[n++]=1;
		}
	}

	for (i=j=m=0;i<MINNUM;i++) {
		if (j<it[i]) {
			j=it[i];
			m=i;
		}
	}

	if (3*j>=MINNUM) {
		return dt[m];
	}
	else {
		for (i=0;i<n;i++) {
			sprintf(PPP_Glo.chMsg,"Sample is %6.2f\n",dt[i]);
			outDebug(OUTWIN,0,0);
		}
		strcpy(PPP_Glo.chMsg,"*** ERROR: sampling may be inaccurate!\n");
		outDebug(OUTWIN,0,0);
	}

	return dt[m];
}
//clock jump repair
static int clkRepair(obsd_t *obs, int n)
{
	int i,sat,validGps,cjGps;
	int bObserved[MAXPRNGPS];
	double delta0=0.0,delta1=0.0,d1,d2,d3,d4,ddd1,ddd2;
	double *lam;
	double CJ_F1,CJ_F2;

	for (i=0;i<MAXPRNGPS;i++) bObserved[i]=0;

	validGps=cjGps=0;

	for (i=0;i<n;i++) {
		sat=obs[i].sat;
		lam=PPP_Glo.lam[sat-1];

		if (sat>MAXPRNGPS) continue;

		if (obs[i].P[0]*obs[i].P[1]*obs[i].L[0]*obs[i].L[1]==0.0) continue;

		if (PPP_Glo.obs0[sat-1][0]*PPP_Glo.obs0[sat-1][1]*PPP_Glo.obs0[sat-1][2]*PPP_Glo.obs0[sat-1][3]==0.0)
			continue;

		validGps++;

		d1=obs[i].P[0]-PPP_Glo.obs0[sat-1][0];
		d2=obs[i].P[1]-PPP_Glo.obs0[sat-1][1];
		d3=(obs[i].L[0]-PPP_Glo.obs0[sat-1][2])*lam[0];
		d4=(obs[i].L[1]-PPP_Glo.obs0[sat-1][3])*lam[1];

		if(fabs(d1-d3)>290000)   //ms clock jump
		{
			delta0+=d1-d3;
			delta1+=d2-d4;
			cjGps++;
		}
	}

	if(cjGps!=0&&cjGps==validGps)
	{
		d1=delta0/cjGps;
		d2=delta1/cjGps;

		CJ_F1=0.0;   //flag for clock jump
		CJ_F2=0.0;
		CJ_F1=d1/CLIGHT*1000.0;
		CJ_F2=myRound(CJ_F1);

		if (fabs(CJ_F1-CJ_F2)<2.5E-2)
		{
			PPP_Glo.clkJump+=(int)CJ_F2;
			sprintf(PPP_Glo.chMsg, "*** WARNING: clock jump=%d(ms)\n",PPP_Glo.clkJump);
			outDebug(OUTWIN,OUTFIL,0);
		}
		else 
		{
			//
		}
	}

	//
	for (i=0;i<n;i++) 
	{
		sat=obs[i].sat;

		if (sat>MAXPRNGPS) continue;

		bObserved[sat-1]=1;

		PPP_Glo.obs0[sat-1][0]=obs[i].P[0];
		PPP_Glo.obs0[sat-1][1]=obs[i].P[1];
		PPP_Glo.obs0[sat-1][2]=obs[i].L[0];
		PPP_Glo.obs0[sat-1][3]=obs[i].L[1];

		ddd1=PPP_Glo.clkJump*CLIGHT/1000.0;
		ddd2=PPP_Glo.clkJump*CLIGHT/1000.0;

		//repair for phase observations
		if(obs[i].L[0]!=0.0) obs[i].L[0]+=ddd1/lam[0];

		if(obs[i].L[1]!=0.0) obs[i].L[1]+=ddd2/lam[1];

		//repair for code observations
		//if( obs[i].P[0]!=0.0 ) {
		//    obs[i].P[0]-=PPP_Glo.clkjump*CLIGHT/1000;
		//}
		//if( obs[i].P[1]!=0.0 ) {
		//    obs[i].P[1]-=PPP_Glo.clkjump*CLIGHT/1000;
		//}
	}

	for (i=0;i<MAXPRNGPS;i++) {
		if (bObserved[i]==0)
			PPP_Glo.obs0[i][0]=PPP_Glo.obs0[i][1]=PPP_Glo.obs0[i][2]=PPP_Glo.obs0[i][3]=0.0;
	}

	return 1;
}
//calculate DOPs {GDOP,PDOP,HDOP,VDOP}
static double calDop(rtk_t *rtk, const obsd_t *obs, const int n)
{
	double azel[MAXSAT*2],dop[4];
	int i,num,sat;

    for (i=0;i<NSYS_USED;i++) rtk->sol.ns[i]=0;
	for (i=num=0;i<n;i++) {
		sat=obs[i].sat;
		if (rtk->ssat[sat-1].vsat[0]==0) continue;
        if (PPP_Glo.sFlag[sat-1].sys==SYS_GPS) rtk->sol.ns[0]++;
        else if (PPP_Glo.sFlag[sat-1].sys==SYS_GLO) rtk->sol.ns[1]++;
        else if (PPP_Glo.sFlag[sat-1].sys==SYS_CMP) rtk->sol.ns[2]++;
        else if (PPP_Glo.sFlag[sat-1].sys==SYS_GAL) rtk->sol.ns[3]++;
        else if (PPP_Glo.sFlag[sat-1].sys==SYS_QZS) rtk->sol.ns[4]++;
		azel[2*num+0]=rtk->ssat[sat-1].azel[0];
		azel[2*num+1]=rtk->ssat[sat-1].azel[1];
		num++;
	}
	dops(num,azel,0.0,dop);

	rtk->sol.dop[1]=dop[1];

	return dop[1];
}
/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
static void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
    sol_t sol0={{0}};
    ssat_t ssat0={0};
    int i;

    rtk->sol=sol0;

    rtk->nx=pppnx(opt);

    rtk->na=3;
    rtk->tt=0.0;
    rtk->x=zeros(rtk->nx,1);
    rtk->P=zeros(rtk->nx,rtk->nx);
    rtk->xa=zeros(rtk->na,1);
    rtk->Pa=zeros(rtk->na,rtk->na);

    rtk->nfix=0;
    for (i=0;i<MAXSAT;i++) rtk->ssat[i]=ssat0;
    rtk->opt=*opt;
}
static void rtkfree(rtk_t *rtk)
{
	rtk->nx=rtk->na=0;
	if (rtk->x)  {free(rtk->x);  rtk->x =NULL;}
	if (rtk->P)  {free(rtk->P);  rtk->P =NULL;}
	if (rtk->xa) {free(rtk->xa); rtk->xa=NULL;}
	if (rtk->Pa) {free(rtk->Pa); rtk->Pa=NULL;}
}
/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by 
* precise positioning
* args   : rtk_t *rtk       IO  rtk control/result struct
*            rtk->sol       IO  solution
*                .time      O   solution time
*                .rr[]      IO  rover position/velocity
*                               (I:fixed mode,O:single mode)
*                .Qr[]      O   rover position covarinace
*                .stat      O   solution status (SOLQ_???)
*                .ns        O   number of valid satellites
*                .age       O   age of differential (s)
*                .ratio     O   ratio factor for ambiguity validation
*            rtk->rb[]      IO  base station position/velocity
*                               (I:relative mode,O:moving-base mode)
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellite
*            rtk->tt        O   time difference between current and previous (s)
*            rtk->x[]       IO  float states pre-filter and post-filter
*            rtk->P[]       IO  float covariance pre-filter and post-filter
*            rtk->xa[]      O   fixed states after AR
*            rtk->Pa[]      O   fixed covariance after AR
*            rtk->ssat[s]   IO  sat(s+1) status
*                .sys       O   system (SYS_???)
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)
*                .resp [f]  O   freq(f+1) pseudorange residual (m)
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
static int rtkpos(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
	gtime_t time;
	int nu;
	char msg[128]="";
	prcopt_t *opt=&rtk->opt;

	rtk->sol.stat=SOLQ_NONE;
	time=rtk->sol.time;  /* previous epoch */
	PPP_Glo.bOKSPP=1;

	/* rover position by single point positioning */
	if (!spp(obs,n,nav,opt,&rtk->sol,NULL,rtk->ssat,msg)) {
		sprintf(PPP_Glo.chMsg,"*** ERROR: point pos error (%s)\n",msg);
		outDebug(OUTWIN,OUTFIL,0);

		PPP_Glo.bOKSPP=0;
		PPP_Glo.nBadEpSPP++;

		//fewer than 4 satellites available, skip to next epoch
		if (n<=4) {
			return -1;
		}
	}
	

	if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);

	nu=n;
	obsScan_PPP(opt,obs,n,&nu);
	if (nu<=4) {
		sprintf(PPP_Glo.chMsg, "*** WARNING: There are only %d satellites observed, skip PPP!\n",nu);
		outDebug(OUTWIN,OUTFIL,0);
		return 0;
	}

	//clock jump repair
	clkRepair(obs,nu);

	/* precise point positioning */
	if (opt->mode>=PMODE_PPP_KINEMA) {
		pppos(rtk,obs,nu,nav);
	}
	else return 1;

	//calculate DOPs
	calDop(rtk,obs,nu);

	//save the information for current epoch
	keepEpInfo(rtk,obs,nu,nav);

	return 1;
}
/* process positioning -------------------------------------------------------*/
static void procpos(rtk_t *rtk,const prcopt_t *popt,const solopt_t *sopt,int mode)
{
	sol_t sol={{0}};
	gtime_t time={0};
	obsd_t obs[MAXOBS];
	int i,j,k=0,nep=0,nobs,n,solstatic,pri[]={0,1,2,3,4,5,1,6};

	solstatic=sopt->solstatic&&popt->mode==PMODE_PPP_STATIC;

	/* processing epoch-wise */
	while ((nobs=inputobs(obs,obss,PPP_Glo.revs,&PPP_Glo.iObsu,&PPP_Glo.iEpoch))>=0) {
		PPP_Glo.tNow=obs[0].time;
		time2epoch(PPP_Glo.tNow,PPP_Glo.ctNow);
		sprintf(PPP_Glo.chTime,"%02.0f:%02.0f:%04.1f%c",PPP_Glo.ctNow[3],
			PPP_Glo.ctNow[4],PPP_Glo.ctNow[5],'\0');
		PPP_Glo.sowNow=time2gpst(PPP_Glo.tNow,NULL);

		k++;
		if (k==1) {
			for (j=0;j<MAXSAT;j++) {
				PPP_Glo.ssat_Ex[j].tLast=PPP_Glo.tNow;
			}
		}
        nep=(int)(30*(60/PPP_Glo.sample));
        if ((k-1)%nep==0) PPP_Glo.t_30min=PPP_Glo.tNow;

		//pseudorange observation checking
		obsScan_SPP(popt,obs,nobs,&n);
		if (n<=3) {
			sprintf(PPP_Glo.chMsg,"*** WARNING: There are only %d satellites observed, skip SPP!\n",n);
			outDebug(OUTWIN,OUTFIL,0);
			continue;
		}

		if (PPP_Glo.prcOpt_Ex.navSys&SYS_CMP) {
			BDmultipathCorr(rtk,obs,n);
		}

		i=rtkpos(rtk,obs,n,&navs);
		if (i==-1) rtk->sol.stat=SOLQ_NONE;
		else if (i==0) continue;

		if (mode==0) {  /* forward/backward */
			outResult(rtk,sopt);

			if (!solstatic&&PPP_Glo.outFp[0]) 
				outsol(PPP_Glo.outFp[0],&rtk->sol,sopt,PPP_Glo.iEpoch);
			else if (time.time==0||pri[rtk->sol.stat]<=pri[sol.stat]) {
				sol=rtk->sol;
				if (time.time==0||timediff(rtk->sol.time,time)<0.0) 
					time=rtk->sol.time;
			}
		}
		
	}
}
/* execute processing session ------------------------------------------------*/
static int execses(prcopt_t *popt, const solopt_t *sopt, filopt_t *fopt)
{
	rtk_t rtk;

	//to determine the thresh values for cycle slip detection
	PPP_Glo.sample=sampledetermine(popt);
	if (fabs(PPP_Glo.prcOpt_Ex.csThresGF)<0.01||fabs(PPP_Glo.prcOpt_Ex.csThresMW)<0.01)
		calCsThres(popt,PPP_Glo.sample);

	rtkinit(&rtk,popt);

	if (PPP_Glo.outFp[OFILE_IPPP]) outiPppHead(PPP_Glo.outFp[OFILE_IPPP],rtk);

	if (popt->soltype==0) {       /* forward */
		PPP_Glo.revs=0;
		PPP_Glo.iObsu=0;
		PPP_Glo.iEpoch=0;
		procpos(&rtk,popt,sopt,0); 
	}
	else if (popt->soltype==1) {  /* backward */
		PPP_Glo.revs=1; 
		PPP_Glo.iObsu=obss.n-1;
		PPP_Glo.iEpoch=PPP_Glo.nEpoch;
		procpos(&rtk,popt,sopt,0);
	}
	else {   /* combined */
		PPP_Glo.solf=(sol_t *)malloc(sizeof(sol_t)*PPP_Glo.nEpoch);
		PPP_Glo.solb=(sol_t *)malloc(sizeof(sol_t)*PPP_Glo.nEpoch);
		if (PPP_Glo.solf&&PPP_Glo.solb) {

		}
		else
			printf("error : memory allocation");

		free(PPP_Glo.solf); PPP_Glo.solf=NULL;
		free(PPP_Glo.solb); PPP_Glo.solb=NULL;
	}

	rtkfree(&rtk);
    
    return 0;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise 
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first 
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extention as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .lex,.LEX            : qzss lex message log files
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              .fcb,.FCB            : satellite fcb
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int gampPos(gtime_t ts, gtime_t te, double ti, double tu, prcopt_t *popt, 
	               const solopt_t *sopt, filopt_t *fopt)
{
    int i,j,stat=0,index[MAXINFILE]={0};

	/* write header to output file */
	if (!outhead(fopt->outf,popt,sopt,PPP_Glo.outFp,MAXOUTFILE)) 
		return 0;
    
	for (i=0;i<MAXOUTFILE;i++) {
		if (fopt->outf[i]&&strlen(fopt->outf[i])>2) 
			PPP_Glo.outFp[i]=openfile(fopt->outf[i]);
		else
			PPP_Glo.outFp[i]=NULL;
	}

	/* set rinex code priority for precise clock */
	if (PMODE_PPP_KINEMA<=popt->mode)
		setcodepri(SYS_GPS,1,popt->sateph==EPHOPT_PREC?"PYWC":"CPYW");

	/* read satellite antenna parameters */
	if (*fopt->antf&&!(readpcv(fopt->antf,&pcvss))) {
		printf("*** ERROR: no sat ant pcv in %s\n",fopt->antf);
		return -1;
	}

	/* read dcb parameters */
	for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
		navs.cbias[i][j]=0.0;
	}
	if (*fopt->p1p2dcbf) 
		readdcb(fopt->p1p2dcbf,&navs);
	if (*fopt->p1c1dcbf) 
		readdcb(fopt->p1c1dcbf,&navs);
	if (*fopt->p2c2dcbf) 
		readdcb(fopt->p2c2dcbf,&navs);
	if (*fopt->mgexdcbf&&(popt->navsys&SYS_CMP||popt->navsys&SYS_GAL))
		readdcb_mgex(fopt->mgexdcbf,&navs,PPP_Glo.prcOpt_Ex.ts);

	/* read erp data */
	if (*fopt->eopf) {
		if (!readerp(fopt->eopf,&navs.erp)) {
			printf("ERROR: no erp data %s\n",fopt->eopf);
		}
	}

	/* read ionosphere data file */
	if (*fopt->ionf&&(popt->ionoopt==IONOOPT_TEC||((popt->ionoopt==IONOOPT_UC1||popt->ionoopt==IONOOPT_UC12)&&
		PPP_Glo.prcOpt_Ex.ion_const)))
		readtec(fopt->ionf,&navs,1);

	for (i=0;i<MAXINFILE;i++) index[i]=i;

	/* read prec ephemeris */
	readpreceph(fopt->inf,MAXINFILE,popt,&navs);

	/* read obs and nav data */
	if (!readobsnav(ts,te,ti,fopt->inf,index,MAXINFILE,popt,&obss,&navs,stas)) {
		freeobsnav(&obss,&navs);
		return 0;
	}

	if (PPP_Glo.nEpoch<=1) {
		strcpy(PPP_Glo.chMsg,"PPP_Glo.nEpoch<=1!\n\0");
		printf("%s",PPP_Glo.chMsg);
		freeobsnav(&obss,&navs);

		return 0;
	}

	//read igs antex only once, and save the elements in 'pcvss'
	/* set antenna paramters */
	setpcv(obss.data[0].time,popt,&navs,&pcvss,&pcvss,stas);

	/* read ocean tide loading parameters */
	if (popt->mode>PMODE_SINGLE&&fopt->blqf) {
		readotl(popt,fopt->blqf,stas);
	}

	//next processing
	stat=execses(popt,sopt,fopt);

	/* free obs and nav data */
	freeobsnav(&obss,&navs);
	/* free prec ephemeris and sbas data */
	freepreceph(&navs);
	/* free antenna parameters */
	if (pcvss.pcv) {free(pcvss.pcv); pcvss.pcv=NULL; pcvss.n=pcvss.nmax=0;}
	if (pcvsr.pcv) {free(pcvsr.pcv); pcvsr.pcv=NULL; pcvsr.n=pcvsr.nmax=0;}
	/* free erp data */
	if (navs.erp.data) {free(navs.erp.data); navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;}

	if (PPP_Glo.outFp[OFILE_IPPP]) fprintf(PPP_Glo.outFp[OFILE_IPPP],"-PPP_BLOCK\n");
	for (i=0;i<MAXOUTFILE;i++) {
		if (i==OFILE_DEBUG) continue;
		if (PPP_Glo.outFp[i]) {
			fclose(PPP_Glo.outFp[i]);
			PPP_Glo.outFp[i]=NULL;
		}
	}

	return stat;
}
