/*------------------------------------------------------------------------------
* spp.c : standard point positioning
*-----------------------------------------------------------------------------*/
#include "gamp.h"

/* constants -----------------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#define NX_SPP      (4+4)       /* # of estimated parameters */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */


/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;

    fact=sys==SYS_GLO?PPP_Glo.prcOpt_Ex.errRatioGLO:(sys==SYS_CMP?PPP_Glo.prcOpt_Ex.errRatioBDS:
		(sys==SYS_GAL?PPP_Glo.prcOpt_Ex.errRatioGAL:(sys==SYS_QZS?PPP_Glo.prcOpt_Ex.errRatioQZS:opt->err[0])));
    varr=(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));

    if (opt->ionoopt==IONOOPT_IF12) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}
/* get tgd parameter (m) -----------------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, double *tgd1, double *tgd2)
{
    int i;
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;

        if (PPP_Glo.sFlag[sat-1].sys==SYS_CMP) {
            if (tgd1) *tgd1=CLIGHT*nav->eph[i].tgd[0];
            if (tgd2) *tgd2=CLIGHT*nav->eph[i].tgd[1];
            return CLIGHT*(nav->eph[i].tgd[0]);
        }

        return CLIGHT*nav->eph[i].tgd[0];
    }
    return 0.0;
}

/* psendorange with code bias correction -------------------------------------*/
//static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt)
static double prange(const obsd_t *obs, const nav_t *nav, const double *azel,
                     const prcopt_t *opt, double *var)
{
    const double *lam=nav->lam[obs->sat-1];
    double PC,P1,P2,P1_P2,P1_C1,P2_C2,gamma,tgd1=0.0,tgd2=0.0;
    int i=0,j=1,sys;
    
    *var=0.0;
    
    if (!(sys=satsys(obs->sat,NULL))) return 0.0;
    
    /* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
    //if (NFREQ>=3&&(sys&(SYS_GAL|SYS_SBS))) j=2;
    
    if (NFREQ<2||lam[i]==0.0||lam[j]==0.0) return 0.0;
    
    gamma=SQR(lam[j])/SQR(lam[i]); /* f1^2/f2^2 */
    P1=obs->P[i];
    P2=obs->P[j];
    P1_P2=nav->cbias[obs->sat-1][0];
    P1_C1=nav->cbias[obs->sat-1][1];
    P2_C2=nav->cbias[obs->sat-1][2];
    
    /* if no P1-P2 DCB, use TGD instead */
	P1_P2=0.0;
    if (P1_P2==0.0&&(sys&(SYS_GPS|SYS_GAL|SYS_QZS|SYS_CMP))) {
        if (sys==SYS_CMP) 
            P1_P2=(1.0-gamma)*gettgd(obs->sat,nav,&tgd1,&tgd2);
        else
            P1_P2=(1.0-gamma)*gettgd(obs->sat,nav,NULL,NULL);
    }
    if (opt->ionoopt==IONOOPT_IF12) { /* dual-frequency */
        if (P1==0.0||P2==0.0) return 0.0;
        if (obs->code[i]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
        if (obs->code[j]==CODE_L2C) P2+=P2_C2; /* C2->P2 */
        
        /* iono-free combination */
        PC=(gamma*P1-P2)/(gamma-1.0);

        if (sys==SYS_CMP) PC=PC+(tgd2-gamma*tgd1)/(gamma-1.0);
    }
    else { /* single-frequency */
        if (P1==0.0) return 0.0;
        if (obs->code[i]==CODE_L1C) P1+=P1_C1;   /* C1->P1 */
        PC=P1-P1_P2/(1.0-gamma);
    }
    
    *var=SQR(ERR_CBIAS);
    
    return PC;
}

/* ionospheric correction ----------------------------------------------------*/
static int ionocorr(gtime_t time, const int sys, const nav_t *nav, const double *pos,
                    const double *azel, const prcopt_t *opt, double *ion,
                    double *var)
{
    /* broadcast model */
    if (opt->ionoopt==IONOOPT_BRDC) {
        //*ion=ionmodel(time,nav->ion_gps,pos,azel);

        if (sys==SYS_GPS) *ion=ionmodel(time,nav->ion_gps,pos,azel);
        //else if ( SYS_CMP==sys )  *ion=ionmodel(time,nav->ion_cmp,pos,azel);
        else *ion=ionmodel(time,nav->ion_gps,pos,azel);

        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
	/* ionex tec model */
	else if (opt->ionoopt==IONOOPT_TEC) {
		return iontec(time,nav,pos,azel,3,ion,var);
	}
    else if (opt->ionoopt==IONOOPT_IF12) {
        *ion=0.0;
        *var=SQR(0.02);
        return 1;
    }

    *ion=0.0;
    *var=opt->ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction ---------------------------------------------------*/
static int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, const prcopt_t *opt, double *trp,
                    double *var)
{
    double trpw=0.0;
    *trp=0.0;

    /* saastamoinen model */
    if (opt->tropopt==TROPOPT_SAAS||opt->tropopt==TROPOPT_EST) {
        *trp=tropmodel(time,pos,azel,REL_HUMI,&trpw,0);
        *trp+=trpw;
        *var=SQR(ERR_SAAS);

        if (*trp>100.0) *trp=100.0;
        else if (*trp<0.05) *trp=0.05;

        return 1;
    }

    /* no correction */
    *trp=0.0;
    *var=opt->tropopt==TROPOPT_OFF?0.0:SQR(ERR_TROP);
    
	return 1;
}

static int getHVR_s(int bMulGnss, const int iter, int *sat, int *ix, double *v, double *var, int nv, 
	                int *ibadsn, double std_ex, double ave_ex, int bElevCVG)
{
    double dValue,factor=1.0,dt;
    int i,j,ind,nbad=2,nBadRes,ibadsum=0,sn[MAXOBS];
	int b,b1,b2,b3,id[100];
	double dtmp,dFactor[100];

    if (!bElevCVG) {
        factor=3.0;
        factor=2.25;
    }

    if (nv>=14)      nbad=5;
    else if (nv>=11) nbad=4;
    else if (nv>=8)  nbad=3;
    else             nbad=2;

    nBadRes=findGross(1,bMulGnss,v,nv,nbad,&std_ex,&ave_ex,sn,5.0,1.0,2.5);

    dValue=factor*std_ex;
    if (std_ex<5000.0) {
        for (i=ibadsum=0;i<nBadRes;i++) {
            b=0;
            j=sn[i];
            ind=j;
            dt=fabs(v[ind]-ave_ex);
            if (fabs(dt)<1.0) continue;

            if (!bElevCVG) {
                if (fabs(dt)<10.0) continue;
            }

            if (dt>800.0)     {if (dt> 10*dValue) b=1;}
            else if (dt>20.0) {if (dt>5.0*dValue) b=1;}
            else if (dt>10.0) {if (dt>6.0*dValue) b=1;}
            else if (dt>3.0)  {if (dt>7.0*dValue) b=1;}

            if (b==0) continue;

            ibadsn[ibadsum]=ix[ind];
            ibadsum++;
        }
    }
    else {
        if (iter>=0) {
            for (i=0;i<100;i++) {
                id[i]=-1;
                dFactor[i]=0.0;
            }

            dtmp=1.0/std_ex;
            for (i=0;i<nBadRes;i++) {
                j=sn[i];
                ind=j;
                id[i]=ix[ind];
                dFactor[i]=(v[ind]-ave_ex)*dtmp;
            }

            for (i=0;i<nBadRes;i++) {
                for (j=i+1;j<nBadRes;j++) {
                    if (fabs(dFactor[i])>=fabs(dFactor[j])) continue;

                    dtmp=dFactor[i];    dFactor[i]=dFactor[j];  dFactor[j]=dtmp;
                    ind=id[i];          id[i]=id[j];            id[j]=ind;
                }
            }

            for (j=0;j<nBadRes;j++) {
                if (fabs(dFactor[j])<=3.0) break;
            }

            if (nBadRes>1) {
                if (j==0) {
                    return 0;
                }

                b1=fabs(dFactor[j-1])>15.0;
                b2=fabs(dFactor[j-1])>fabs(2.0*dFactor[j]);
                b3=fabs(dFactor[j-1])-fabs(dFactor[j])>2.0;

                if (b1&&b2&&b3) {
                    for (i=ibadsum=0;i<=j-1;i++) {
                        ibadsn[ibadsum]=id[i];
                        ibadsum++;
                    }
                }
            }
            else if (nBadRes==1) {
                b1=fabs(dFactor[0])>25.0;
                if (b1) {
                    ibadsn[0]=id[0];
                    ibadsum=1;
                }
            }
        }
    }

    return ibadsum;
}

static int getHVR_spp(int bMulGnss, const int iter, const int sys, int bElevCVG, int bDeleted[MAXSAT], 
	                  int *sat, double *H, double *v, double *var, double *elev, int nv, int nx)
{
	int bbad;
	double *dv,*dvar,std_ex[NSYS_USED],ave_ex[NSYS_USED];
	int ibadsum[NSYS_USED],ibadsn[NSYS_USED][MAXOBS],navsys[8]={SYS_GPS,SYS_GLO,SYS_CMP,SYS_GAL,SYS_QZS,0};
	int i,j,k,newnv=0,ibadsum_,*ix;
	int nv_nx;

    if (nv<=4) return nv;

    bbad=0;
    dv=mat(nv,1);
    dvar=mat(nv,1);
    ix=imat(nv,1);

    //////////////////////////////////////////////////////////////////////////
    for (i=0;i<NSYS_USED;i++) {
        ibadsum[i]=0;
        std_ex[i]=ave_ex[i]=0.0;

        for (j=0;j<MAXOBS;j++) ibadsn[i][j]=-1;

        if (!(sys&navsys[i])) continue;

        for (k=j=0;k<nv;k++) {
            if (navsys[i]==PPP_Glo.sFlag[sat[k]-1].sys) {
                ix[j]=k;
                dv[j]=v[k];
                dvar[j]=var[k];
                j++;
            }
        }
		ibadsum[i]=getHVR_s(bMulGnss,iter,sat,ix,dv,dvar,j,ibadsn[i],std_ex[i],ave_ex[i],bElevCVG);
    }
    //////////////////////////////////////////////////////////////////////////

    free(dv); free(dvar); free(ix);

    //////////////////////////////////////////////////////////////////////////
    ibadsum_=0;
    k=ibadsum[0];
    for (i=0;i<NSYS_USED;i++) {
        ibadsum_+=ibadsum[i];

        if (i==0) continue;

        for (j=0;j<ibadsum[i];j++)
            ibadsn[0][k++]=ibadsn[i][j];
    }
    //////////////////////////////////////////////////////////////////////////
    if (ibadsum_>0) {
        if (ibadsum_>=3&&nv-ibadsum_>=6) bbad=1;
        if (ibadsum_<=2) {
            if (nv<=6) {
                if (nv-ibadsum_>4) {
                    bbad=1;
                }
            }
            else if (nv-ibadsum_>=5) {
                bbad=1;
            }
        }
    }

    if (bbad==0) return nv;

    for (i=j=0;i<nv;i++) {
        v[newnv]=v[i];
        var[newnv]=var[i];

        nv_nx=newnv*nx;
        for (k=0;k<nx;k++) H[nv_nx+k]=H[i*nx+k];

        if (j<ibadsum_) {
            if (ibadsn[0][j]==i) {
                j++;

                continue; 
            }
        }
        newnv++;
    }

    nv=newnv;

    return nv;
}

/* pseudorange residuals -----------------------------------------------------*/
static int rescode(const int iter, int bElevCVG, const obsd_t *obs, int n, const double *rs, const double *dts, 
	               const double *vare, const int *svh, const nav_t *nav, const double *x, const prcopt_t *opt, 
				   double *v, double *H, double *var, double *azel, int *vsat, double *resp, int *nx, int *bDeleted)
{
    int bObserved[5];
    int i,j,nv=0,ns[5]={0},sys,satsn[MAXOBS],sat;
    double r,dion,dtrp,vion,vtrp,rr[3],pos[3],dtr,e[3],P,elev_t[MAXSAT],vmeas,lam_L1;
	double elmin;
	int bMulGNSS;
    
    *nx=NX_SPP;

    for (i=0;i<5;i++)            {bObserved[i]=0; ns[i]=0;}
    for (i=0;i<n;i++)            v[i]=var[i]=0.0;
    for (i=0;i<NX_SPP*(n+5);i++) H[i]=0.0;
    for (i=0;i<MAXSAT;i++)       elev_t[i]=0.0;
    for (i=0;i<3;i++)            rr[i]=x[i]; dtr=x[3];
    
    ecef2pos(rr,pos);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;

        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;

        sys=PPP_Glo.sFlag[sat-1].sys;

        if (bDeleted[sat-1]==0) continue;
        if (!(sys&opt->navsys)) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
            sprintf(PPP_Glo.chMsg,"*** WARNING: duplicated observation data %s sat=%2d\n",
				time_str(obs[i].time,3),obs[i].sat);
			outDebug(OUTWIN,OUTFIL,0);
            i++;
            continue;
        }
        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        satazel(pos,e,azel+i*2);

        if (bElevCVG) {
            if (PPP_Glo.prcOpt_Ex.bElevCheckEx) {
                elmin=0.0;
                elmin=MIN(opt->elmin,2.0*D2R);
                elmin=MAX(opt->elmin/3.0,elmin);
                if (azel[1+i*2]<elmin) continue;
            }
            else {
                if (azel[1+i*2]<opt->elmin) continue;
            }
        }
        else {
            if (azel[1+i*2]<=0.0)
                azel[1+i*2]=MIN(3.0*D2R,opt->elmin*0.75);

            if (iter>=1) {
                if (azel[1+i*2]<opt->elmin) continue;
            }
        }
        
        /* psudorange with code bias correction */
        if ((P=prange(obs+i,nav,azel+i*2,opt,&vmeas))==0.0) continue;

        /* excluded satellite */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;
        
        /* ionospheric corrections */
        if (!ionocorr(obs[i].time,sys,nav,pos,azel+i*2,opt,&dion,&vion)) continue;
		/* GPS-L1 -> L1/B1 */
		if ((lam_L1=nav->lam[obs[i].sat-1][0])>0.0) {
			dion*=SQR(lam_L1/lam_carr[0]);
		}

        /* tropospheric corrections */
        if (!tropcorr(obs[i].time,nav,pos,azel+i*2,opt,&dtrp,&vtrp)) continue;

        /* pseudorange residual */
        v[nv]=resp[i]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        /* design matrix */
        for (j=0;j<4;j++) H[j+nv*NX_SPP]=j<3?-e[j]:1.0;
        
        /* time system and receiver bias offset */
        if (sys==SYS_GLO)      {v[nv]-=x[4];   H[4+nv*NX_SPP]=1.0; ns[1]++; bObserved[1]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[5];   H[5+nv*NX_SPP]=1.0; ns[2]++; bObserved[2]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[6];   H[6+nv*NX_SPP]=1.0; ns[3]++; bObserved[3]=1;}
		else if (sys==SYS_QZS) {v[nv]-=x[7];   H[7+nv*NX_SPP]=1.0; ns[4]++; bObserved[4]=1;}
        else                   {H[4+nv*NX_SPP]=H[5+nv*NX_SPP]=0.0; ns[0]++; bObserved[0]=1;}
        
        vsat[i]=1; resp[i]=v[nv];

        satsn[nv]=obs[i].sat;
        elev_t[sat-1]=azel[1+i*2]*R2D;
        
        /* error variance */
        var[nv]=varerr(opt,azel[1+i*2],sys)+vare[i]+vion+vtrp;

        if (azel[1+i*2]<opt->elmin) var[nv]*=100.0;

        nv++;
    }


    for (i=j=0;i<5;i++) {
        if (bObserved[i]) j++;
    }

    bMulGNSS=0;
    if (j>=2) bMulGNSS=1;

    i=nv;
    nv=getHVR_spp(bMulGNSS,iter,opt->navsys,bElevCVG,bDeleted,satsn,H,v,var,elev_t,nv,*nx);

    for (i=0;i<5;i++) {
        if (ns[i]>0) continue;
        for (j=0;j<nv;j++) 
            H[(3+i)+j*(*nx)]=0.0;
    }

    for (i=0;i<5;i++) {
        if (ns[i]==0) *nx=*nx-1;
    }

    //*nx=*nx-1;

    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg, double dop[4])
{
    double azels[MAXOBS*2],vv;
    int i,ns,bGdopOK,bPdopOK;
	double elmin;
    
    /* chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx) {
        if (vv>chisqr[nv-nx-1]) {
            sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
            return 0;
        }
    }

    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }

    elmin=opt->elmin;

    if (PPP_Glo.prcOpt_Ex.bElevCheckEx) {
        elmin=0.0;
    }

    bGdopOK=1;
    bPdopOK=1;
    dops(ns,azels,elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(PPP_Glo.chMsg,"*** ERROR: gdop error nv=%d gdop=%.1f\n",nv,dop[0]);
        outDebug(OUTWIN,OUTFIL,OUTTIM);

        if (ns>=6) bGdopOK=0;
        else {
            if (dop[0]<=50.0) bGdopOK=1;
            else bGdopOK=0;
        }
    }
    if (dop[1]<=0.0||dop[1]>50.0) {
        sprintf(PPP_Glo.chMsg,"*** ERROR: pdop error nv=%d pdop=%.1f\n",nv,dop[1]);
        outDebug(OUTWIN,OUTFIL,OUTTIM);

        if (ns>=6) bPdopOK=0;
        else {
            if (dop[1]<=60.0) bPdopOK=1;
            else bPdopOK=0;
        }
    }
    if (!bPdopOK||!bGdopOK) {
        //return 0;
    }

    return 1;
}

static int estpos_(int *bDeleted, double *x,const obsd_t *obs, int n, const double *rs, const double *dts, const double *vare, 
                   const int *svh, const nav_t *nav, const prcopt_t *opt, sol_t *sol, double *azel, int *vsat, double *resp, char *msg)
{
    int i,j,k,info,stat=0,nx,nv,bElevCVG;
    double dx[NX_SPP],Q[NX_SPP*NX_SPP],dop[4],*v,*H,*var,sig,d0;

    bElevCVG=0;

    v=mat(n+5,1); H=mat(NX_SPP,n+5); var=mat(n+5,1);

    for (i=0;i<NX_SPP;i++) dx[i]=0.0;

    for (i=0;i<MAXITR;i++) {
        /* pseudorange residuals */
        nv=rescode(i,bElevCVG,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,&nx,bDeleted);

        if (nv<nx) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }

        /* weight by variance */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX_SPP;k++) H[k+j*NX_SPP]/=sig;
        }

        /* least square estimation */
        if ((info=lsqPlus(H,v,NX_SPP,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            sprintf(PPP_Glo.chMsg, "%s\n",msg);
            outDebug(OUTWIN,OUTFIL,OUTTIM);

            break;
        }

        for (j=0;j<NX_SPP;j++) x[j]+=dx[j];

        d0=norm(dx,NX_SPP);

        if (d0<1E4) bElevCVG=1;
        else        bElevCVG=0;

        if (d0<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* glo-gps time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* bds-gps time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* gal-gps time offset (s) */
			sol->dtr[4]=x[7]/CLIGHT; /* qzs-gps time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX_SPP];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX_SPP]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns[0]=(unsigned char)nv;
            sol->rms=sol->dop[0]=sol->dop[1]=sol->dop[2]=sol->dop[3]=0.0;

            /* validate solution */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,nx,msg,dop)))
                sol->stat=SOLQ_SINGLE;

            sol->dop[0]=dop[0];
            sol->dop[1]=dop[1];
            sol->dop[2]=dop[2];
            sol->dop[3]=dop[3];

            break;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);

    free(v); free(H); free(var);
    return stat;
}

/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int nobs, const double *rs, const double *dts, const double *vare, const int *svh, 
                  const nav_t *nav, const prcopt_t *opt, sol_t *sol, double *azel, int *vsat, double *resp, char *msg)
{
    int bDeleted[MAXSAT];
    int i,j,stat=0,n,nb=0,*it,nMin=4;
    double x[NX_SPP]={0},x_[NX_SPP],dt=0.0;

    for (i=0;i<NX_SPP;i++) x[i]=x_[i]=1.0;
    //for (i=0;i<MAXSAT;i++) bDeleted[i]=true;
    
    if (PPP_Glo.crdTrue[0]==0.0) {
        dt=fabs(timediff(PPP_Glo.tNow,sol->time));
        PPP_Glo.delEp=myRound(dt/PPP_Glo.sample);

        if (dt>1800&&PPP_Glo.delEp>100) {
            for (i=0;i<3;i++) x[i]=100.0;
        }
        else {
			for (i=0;i<3;i++) x[i]=sol->rr[i];
			if (norm(x,3)<=10) 
				for (i=0;i<3;i++) x[i]=100.0;
        }
    }
    else {
        for (i=0;i<3;i++) x[i]=PPP_Glo.crdTrue[i];
    }

    sol->time=obs[0].time;

    for (i=0;i<NX_SPP;i++) x_[i]=x[i];

    for (nb=0;nb<=3;nb++) {
        if (nobs-nb<nMin) {
            sprintf(msg,"lack of valid sats ns=%d/%d",nobs,nb);
            break;
        }

        for (i=0,n=1;i<nb;i++)
            n=n*(nobs-i)/(i+1);

		if (nb<=0) it=imat(n,1);
		else       it=imat(n*nb,1);

        comb_j=0;

        select_combination(0,0,nobs,nb,it);
        //////////////////////////////////////////////////////////////////////////

        for (i=stat=0;i<n;i++) {
            for (j=0;j<nobs;j++) bDeleted[obs[j].sat-1]=1;
            for (j=i*nb;j<i*nb+nb;j++) bDeleted[obs[it[j]-1].sat-1]=0;
            for (j=0;j<NX_SPP;j++) x_[j]=x[j];

            stat=estpos_(bDeleted,x_,obs,nobs,rs,dts,vare,svh,nav,opt,sol,azel,vsat,resp,msg);

            if (stat==1) break;
        }

        free(it);

        if (stat==1) break;
    }


    if (stat==1) {
        for (j=0;j<NX_SPP;j++) x[j]=x_[j];
        return 1;
    }
    else
        return 0;

}
static int estpos(const obsd_t *obs, int nobs, const double *rs, const double *dts, const double *vare, const int *svh, const nav_t *nav, const prcopt_t *opt, sol_t *sol, double *azel, int *vsat, double *resp, char *msg);
//static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts, const double *vare, const nav_t *nav, const prcopt_t *opt, sol_t *sol, double *azel, int *vsat, double *resp, char *msg);

/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int spp(const obsd_t *obs, int n, const nav_t *nav,const prcopt_t *opt, 
	sol_t *sol, double *azel, ssat_t *ssat, char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,sat,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {strcpy(msg,"no observation data"); return 0;}
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);

    opt_.sateph =EPHOPT_BRDC;
    opt_.ionoopt=IONOOPT_BRDC;
    opt_.tropopt=TROPOPT_SAAS;

    /* satellite positons, velocities and clocks */
    satposs_rtklib(obs[0].time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    opt_.sateph =EPHOPT_BRDC;
    opt_.ionoopt=IONOOPT_BRDC;
    opt_.tropopt=TROPOPT_SAAS;

    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp_pos[0]=ssat[i].resc_pos[0]=0.0;
            ssat[i].snr[0]=0;
        }

        for (i=0;i<NSYS;i++) sol->ns[i]=0;

        for (i=0;i<n;i++) {
            if (!vsat[i]) continue;
            sat=obs[i].sat;
            ssat[sat-1].vs=1;
            ssat[sat-1].azel[0]=azel_[  i*2];
            ssat[sat-1].azel[1]=azel_[1+i*2];
            ssat[sat-1].resp_pos[0]=resp[i];
            ssat[sat-1].snr[0]=obs[i].SNR[0];

            if (PPP_Glo.sFlag[sat-1].sys==SYS_GPS) sol->ns[0]++;
            else if (PPP_Glo.sFlag[sat-1].sys==SYS_GLO) sol->ns[1]++;
            else if (PPP_Glo.sFlag[sat-1].sys==SYS_CMP) sol->ns[2]++;
            else if (PPP_Glo.sFlag[sat-1].sys==SYS_GAL) sol->ns[3]++;
			else if (PPP_Glo.sFlag[sat-1].sys==SYS_QZS) sol->ns[4]++;
            else { 
                sprintf(PPP_Glo.chMsg,"*** WARNING: unsupported satellite system %d %d!\n", PPP_Glo.sFlag[sat-1].sys,sat); 
                outDebug(OUTWIN,OUTFIL,OUTTIM);
            }
        }
    }

    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}

