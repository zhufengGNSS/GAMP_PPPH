/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*-----------------------------------------------------------------------------*/
#include "gamp.h"

#define SQR(x)      ((x)*(x))

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* extension to sp3-c */
    if (code=='J') return SYS_QZS; /* extension to sp3-c */
    if (code=='C') return SYS_CMP; /* extension to sp3-c */
    if (code=='L') return SYS_LEO; /* extension to sp3-c */
    return SYS_NONE;
}
/* read sp3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
                    double *bfact, char *tsys)
{
    int i,j,k=0,ns=0,sys,prn;
    char buff[1024];
    
    for (i=0;i<22;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;
        
        if (i==0) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
        }
        else if (2<=i&&i<=6) {
            if (i==2) {
                ns=(int)str2num(buff,4,2);
            }
            for (j=0;j<17&&k<ns;j++) {
                sys=code2sys(buff[9+3*j]);
                prn=(int)str2num(buff,10+3*j,2);
                if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
        }
        else if (i==12) {
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (i==14) {
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
    }
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;
    
    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            printf("readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read sp3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
	char *tsys, int index, int opt, nav_t *nav)
{
	peph_t peph;
	gtime_t time;
	double val,std,base;
	int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
	char buff[1024];
	int bvalid=1;

	while (!feof(fp)) {
		if (bvalid) {
			fgets(buff,sizeof(buff),fp);
		}
		bvalid=1;

		if (!strncmp(buff,"EOF",3)) break;

		if (buff[0]!='*'||str2time(buff,3,28,&time)) {
			continue;
		}
		if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
		peph.time =time;
		peph.index=index;

		for (i=0;i<MAXSAT;i++) for (j=0;j<4;j++) {
			peph.pos[i][j]=0.0;
			peph.std[i][j]=0.0f;
			peph.vel[i][j]=0.0;
			peph.vst[i][j]=0.0f;
		}
		for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {
			if (buff[0]=='*') {
				bvalid=0;
				break;
			}

			if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;

			sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
			prn=(int)str2num(buff,2,2);
			if      (sys==SYS_SBS) prn+=100;
			else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */

			if (!(sat=satno(sys,prn))) continue;

			if (buff[0]=='P') {
				pred_c=strlen(buff)>=76&&buff[75]=='P';
				pred_o=strlen(buff)>=80&&buff[79]=='P';
			}
			for (j=0;j<4;j++) {
				/* read option for predicted value */
				if (j< 3&&(opt&1)&& pred_o) continue;
				if (j< 3&&(opt&2)&&!pred_o) continue;
				if (j==3&&(opt&1)&& pred_c) continue;
				if (j==3&&(opt&2)&&!pred_c) continue;

				val=str2num(buff, 4+j*14,14);
				std=str2num(buff,61+j* 3,j<3?2:3);

				if ((val==0.0)&&(j==3)) {
					//
				}

				//if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
				//    peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
				//}
				//if ((base=bfact[j<3?0:1])>0.0) {
				//    peph.std[sat-1][j]=/*(float)*/(pow(base,std)*(j<3?1E-3:1E-12));
				//}

				if (buff[0]=='P') { /* position */
					if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
						peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
						v=1; /* valid epoch */
					}
					if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
						peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
					}
				}
				else if (v) { /* velocity */
					if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
						peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
					}
					if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
						peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
					}
				}
			}
		}
		if (v) {
			if (!addpeph(nav,&peph)) return;
		}
	}
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void combpeph(nav_t *nav, int opt)
{
	int i,j,k,m;

	qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);

	if (opt&4) return;

	for (i=0,j=1;j<nav->ne;j++) {
		if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {
			for (k=0;k<MAXSAT;k++) {
				if (norm(nav->peph[j].pos[k],4)<=0.0) continue;

				if (norm(nav->peph[i].pos[k],4)>0.0) {
					if (nav->peph[j].pos[k][0]*nav->peph[j].pos[k][1]*nav->peph[j].pos[k][2]*nav->peph[j].pos[k][3]==0.0)
						continue;
				}
				for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
				for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
				for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
				for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
			}
		}
		else if (++i<j) nav->peph[i]=nav->peph[j];
	}

	nav->ne=i+1;
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
	FILE *fp;
	gtime_t time={0};
	double bfact[2]={0};
	int ns,index=0,sats[MAXSAT]={0};
	char type=' ',tsys[4]="";
	const char *ext;

	if (!(ext=strrchr(file,'.'))) return;

	if (!strstr(ext+1,"sp3")&&!strstr(ext+1,".SP3")&&
		!strstr(ext+1,"eph")&&!strstr(ext+1,".EPH")) return;

	if (!(fp=fopen(file,"r"))) {
		return;
	}

	if (strstr(file,"cod")||strstr(file,"COD")) index=10;
	else if (strstr(file,"igs")||strstr(file,"IGS")) index=9;
	else if (strstr(file,"igr")||strstr(file,"IGR")) index=8;
	else if (strstr(file,"gfz")||strstr(file,"GFZ")) index=7;
	else if (strstr(file,"esa")||strstr(file,"ESA")) index=6;
	else if (strstr(file,"iac")||strstr(file,"IAC")) index=-1;
	else index=0;

	/* read sp3 header */
	ns=readsp3h(fp,&time,&type,sats,bfact,tsys);

	/* read sp3 body */
	readsp3b(fp,type,sats,ns,bfact,tsys,index,opt,nav);

	fclose(fp);

	/* combine precise ephemeris */
	if (nav->ne>0) combpeph(nav,opt);
}
/* read satellite antenna parameters -------------------------------------------
* read satellite antenna parameters
* args   : char   *file       I   antenna parameter file
*          gtime_t time       I   time
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : only support antex format for the antenna parameter file
*-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav)
{
    pcvs_t pcvs={0};
    pcv_t pcv0={0},*pcv;
    int i;
    
    if (!readpcv(file,&pcvs)) return 0;
    
    for (i=0;i<MAXSAT;i++) {
        pcv=searchpcv(i+1,"",time,&pcvs);
        nav->pcvs[i]=pcv?*pcv:pcv0;
    }
    free(pcvs.pcv);
    return 1;
}
/* read MGEX dcb parameters file --------------------------------------------------*/
extern int readdcb_mgex(const char *file, nav_t *nav, const gtime_t time)
{
	FILE *fp;
	double cbias;
	char buff[256],id[4];
	int sat,sys,type=0,iy,doy1,doy2;
	gtime_t time1,time2;

	if (!(fp=fopen(file,"r"))) {
		printf("*** ERROR: open MGEX DCB files failed!\n");
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {
		if (strstr(buff,"+BIAS/SOLUTION")) type = 1;
		if (strstr(buff,"-BIAS/SOLUTION")) break;

		if (!type) continue;

		if (strncmp(buff+1,"DSB",3) || strncmp(buff+15,"    ",4)) continue;

		iy = (int)str2num(buff,35,4); doy1 = (int)str2num(buff,40,3); doy2 = (int)str2num(buff,55,3);
		if (iy <= 50) iy += 2000;
		time1 = yrdoy2time(iy,doy1); time2 = yrdoy2time(iy,doy2);
		if (!(timediff(time,time1) >= 0.0 && timediff(time,time2) < 0.0)) continue;
		strncpy(id,buff+11,3); id[3] = '\0';
		sat = satid2no(id);
		sys = satsys(sat,NULL);
		if (sys == SYS_GPS) {
			if (!strncmp(buff+25,"C1C  C5Q",8)) {
				cbias = str2num(buff,71,21);
				nav->cbias[sat-1][3] = cbias*1E-9*CLIGHT; /* ns -> m */
			}
		}
		else if (sys == SYS_CMP) {
			if (!strncmp(buff+25,"C2I  C7I",8)) {
				cbias = str2num(buff,71,21);
				nav->cbias[sat-1][0] = cbias*1E-9*CLIGHT; /* ns -> m */
			}
			else if (!strncmp(buff+25,"C2I  C6I",8)) {
				cbias = str2num(buff,71,21);
				nav->cbias[sat-1][3] = cbias*1E-9*CLIGHT; /* ns -> m */
			}
			else if (!strncmp(buff+25,"C7I  C6I",8)) {
				cbias = str2num(buff,71,21);
				nav->cbias[sat-1][4] = cbias*1E-9*CLIGHT; /* ns -> m */
			}
		}
		else if (sys == SYS_GAL) {
			if (!strncmp(buff+25,"C1X  C5X",8)) {
				cbias = str2num(buff,71,21);
				nav->cbias[sat-1][0] = cbias*1E-9*CLIGHT; /* ns -> m */
			}
			else if (!strncmp(buff+25,"C1X  C7X",8)) {
				cbias = str2num(buff,71,21);
				nav->cbias[sat-1][3] = cbias*1E-9*CLIGHT; /* ns -> m */
			}
		}
	}
	fclose(fp);

	return 1;
}
/* read dcb parameters file --------------------------------------------------*/
static int readdcbf(const char *file, nav_t *nav)
{
	FILE *fp;
	double cbias;
	int sat,type=0;
	char buff[256];

	if (!(fp=fopen(file,"r"))) {
		printf("*** ERROR: open CODE DCB files failed!\n");
		return 0;
	}
	while (fgets(buff,sizeof(buff),fp)) {
		if      (strstr(buff,"DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
		else if (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
		else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;

		if (!type) continue;

		if (!(sat=satid2no(buff))||(cbias=str2num(buff,26,9))==0.0) continue;

		nav->cbias[sat-1][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
	}
	fclose(fp);

	return 1;
}
/* read dcb parameters ---------------------------------------------------------
* read differential code bias (dcb) parameters
* args   : char   *file       I   dcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : currently only p1-c1 bias of code *.dcb file
*-----------------------------------------------------------------------------*/
extern int readdcb(const char *file, nav_t *nav)
{
	int i,n;
	char *efiles[MAXEXFILE]={0};

	//for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
	//    nav->cbias[i][j]=0.0;
	//}
	for (i=0;i<MAXEXFILE;i++) {
		if (!(efiles[i]=(char *)malloc(1024))) {
			for (i--;i>=0;i--) free(efiles[i]);
			return 0;
		}
	}
	n=expath(file,efiles,MAXEXFILE);

	for (i=0;i<n;i++) {
		if (strlen(efiles[i])<2) 
			continue;

		readdcbf(efiles[i],nav);
	}
	for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

	return 1;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
    int i,j;
    
    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
static void posWithEarhRotation(const int k, double pos[3], double p[3][NMAX+1], double dt)
{
	double sinl,cosl;
#if 0
	p[0][k]=pos[0];
	p[1][k]=pos[1];
#else
	/* correciton for earh rotation ver.2.4.0 */
	sinl=sin(OMGE*dt);
	cosl=cos(OMGE*dt);
	p[0][k]=cosl*pos[0]-sinl*pos[1];
	p[1][k]=sinl*pos[0]+cosl*pos[1];
#endif
	p[2][k]=pos[2];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
                   double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl;
    int i,j,k,index,sys;
	int id[NMAX+1],kInd,bBadClk;
    
    rs[0]=rs[1]=rs[2]=dts[0]=0.0;
    
    if (nav->ne<NMAX+1||timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE)
        return 0;

    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; 
    else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;

    for (j=k=0;j<NMAX*50;j++) {
        if (index+j>=0&&index+j<nav->ne&&k<=NMAX) {
            id[k]=index+j;
            t[k]=timediff(nav->peph[id[k]].time,time);
            pos=nav->peph[id[k]].pos[sat-1];
            if (norm(pos,3)>0.0) {
                posWithEarhRotation(k,pos,p,t[k]);
                k++;
            }
        }
        if (k==NMAX+1) break;

        if (index-j>=0&&index-j<nav->ne&&k<=NMAX&&j!=0) {
            id[k]=index-j;
            t[k]=timediff(nav->peph[id[k]].time,time);
            pos=nav->peph[id[k]].pos[sat-1];
            if (norm(pos,3)>0.0) {
                posWithEarhRotation(k,pos,p,t[k]);
                k++;
            }
        }
        if (k==NMAX+1) break;
    }
    if (k<=NMAX) return 0;

    for (i=0;i<=NMAX;i++) {
        for (j=i+1;j<=NMAX;j++) {
            if (t[i]<=t[j]) continue;
            sinl=t[j]; t[j]=t[i];   t[i]=sinl;
            k=id[j];   id[j]=id[i]; id[i]=k;
            for (k=0;k<3;k++) {
                sinl=p[k][j];
                p[k][j]=p[k][i];
                p[k][i]=sinl;
            }
        }
    }

    kInd=0;
    for (i=0;i<=NMAX;i++) {
        if (fabs(t[kInd])<=fabs(t[i])) kInd=i;
    }
    index=id[kInd];

    if (t[0]>900.0||t[NMAX]<-900.0) {
        sprintf(PPP_Glo.chMsg,"%s t[0]=%-5.1f t[%d]=%-5.1f\n",PPP_Glo.sFlag[sat-1].id,t[0],NMAX,t[NMAX]);
        outDebug(0,0,0);
        return 0;
    }

    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }

    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);
        
        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];
    
    bBadClk=0;
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) bBadClk=1;
        std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) bBadClk=1;
        std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];    
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
    }
    else {
        bBadClk=1;
    }
    if (varc) *varc=SQR(std);

    sys=PPP_Glo.sFlag[sat-1].sys;
    if (sys==SYS_CMP&&bBadClk) return 1;

    if (bBadClk) {
        //return 0;
    }
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
                   double *varc)
{
    double t[2],c[2],std;
    int i,j,k,index;
    
    if (nav->nc<2||timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE)
        return 1;

    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];

    for (i=index;i>=0;i--) {
        if (nav->pclk[i].clk[sat-1][0]!=0.0) {
            t[0]=timediff(time,nav->pclk[i].time);
            c[0]=nav->pclk[i].clk[sat-1][0];
            break;
        }
    }

    for (i=index+1;i<nav->nc;i++) {
        if (nav->pclk[i].clk[sat-1][0]!=0.0) {
            t[1]=timediff(time,nav->pclk[i].time);
            c[1]=nav->pclk[i].clk[sat-1][0];
            index=i-1;
            break;
        }
    }

    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->pclk[index+i].std[sat-1][0];

        if (std*CLIGHT>0.05) std=std+EXTERR_CLK*fabs(t[i]);
        else                 std=std*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    else  {
        return 0;
    }

    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant)
{
    const double *lam=nav->lam[sat-1];
    const pcv_t *pcv=nav->pcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0};
    double gamma,C1,C2,dant1,dant2;
    int i,j=0,k=1,sys;
    
    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);
    
    /* unit vectors of satellite fixed coordinates */
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);
    
	sys=satsys(sat,NULL);
    //if (NFREQ>=3&&(sys&(SYS_GAL|SYS_SBS))) k=2;
    if (NFREQ<2||lam[j]==0.0||lam[k]==0.0) return;
    
    gamma=SQR(lam[k])/SQR(lam[j]);
    C1=gamma/(gamma-1.0);
    C2=-1.0 /(gamma-1.0);
    
	if (sys==SYS_GPS) {
		j=0;
		k=1;
	}
	else if (sys==SYS_GLO) {
		j=0+NFREQ;
		k=1+NFREQ;
	}
	else if (sys==SYS_CMP) {
		j=0+2*NFREQ;
		k=1+2*NFREQ;
	}
	else if (sys==SYS_GAL) {
		j=0+3*NFREQ;
		k=1+3*NFREQ;
	}
	else if (sys==SYS_QZS) {
		j=0+4*NFREQ;
		k=1+4*NFREQ;
	}
    /* iono-free LC */
    for (i=0;i<3;i++) {
        dant1=pcv->off[j][0]*ex[i]+pcv->off[j][1]*ey[i]+pcv->off[j][2]*ez[i];
        dant2=pcv->off[k][0]*ex[i]+pcv->off[k][1]*ey[i]+pcv->off[k][2]*ez[i];
        dant[i]=C1*dant1+C2*dant2;
    }
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat postion option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                    double *rs, double *dts, double *var)
{
    double rss[3],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;
    
    if (sat<=0||MAXSAT<sat) return 0;
    
    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !pephclk(time,sat,nav,dtss,&varc)) return 0;
    
    time=timeadd(time,tt);
    if (!pephpos(time,sat,nav,rst,dtst,NULL,NULL)||
        !pephclk(time,sat,nav,dtst,NULL)) return 0;

    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rss,sat,nav,dant);
    }

    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else    /* no precise clock */
        dts[0]=dts[1]=0.0;
    
    *var=vare+varc;

    return 1;
}
