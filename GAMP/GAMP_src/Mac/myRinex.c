#include "gamp.h"

//read time tag of the first epoch
static void getFirstTime(double rnxver, FILE *fp, double cts[6])
{
	char oneline[MAXCHARS];
	int i,ilen;

	for (i=0;i<6;i++) cts[i]=0.0;

	while (!feof(fp)) {
		if (rnxver<=2.99) {
			//read time tag of the first epoch
			fgets(oneline,sizeof(oneline),fp);

			i=myRound(str2num(oneline,28,1));

			if (i==2) continue;
			if (i==4) {
				ilen=myRound(str2num(oneline,30,2));
				for (i=0;i<ilen;i++)
					fgets(oneline,sizeof(oneline),fp);
				continue;
			}
			
			if (strstr(oneline,"COMMENT")||strstr(oneline,"comment")) continue;

			if (oneline[3]!=' '||oneline[6]!=' '||oneline[9]!=' '||
				oneline[12]!=' '||oneline[15]!=' ')
				continue;

			if ((oneline[1]==' '&&oneline[2]==' ')||
				(oneline[4]==' '&&oneline[5]==' ')) 
				continue;

			cts[0]=str2num(oneline,1,2);
			if(cts[0]<=79.0&&cts[0]>=0.0) cts[0]+=2000.0;
			else	 cts[0]+=1900.0;
			cts[1]=str2num(oneline,4,2);
			cts[2]=str2num(oneline,7,2);
			cts[3]=str2num(oneline,10,2);
			cts[4]=str2num(oneline,13,2);
			cts[5]=str2num(oneline,15,11);
		}
		else {
			//read time tag of the first epoch
			fgets(oneline,sizeof(oneline),fp);

			if (oneline[0]!='>') continue;

			if (strstr(oneline,"COMMENT")||strstr(oneline,"comment"))
				continue;

			if (oneline[1]!=' '||oneline[6]!=' '||oneline[9]!=' '||
				oneline[12]!=' '||oneline[15]!=' ')
				continue;

			cts[0]=str2num(oneline,2,4);
			cts[1]=str2num(oneline,7,2);
			cts[2]=str2num(oneline,10,2);
			cts[3]=str2num(oneline,13,2);
			cts[4]=str2num(oneline,16,2);
			cts[5]=str2num(oneline,19,10);

			if (norm(cts,6)<=1.0) continue;
		}
		break;
	}
}

//read time tag of the last epoch
static void getLastTime(double rnxver, FILE *fp, double cte[6])
{
	int i,n,len;
	char oneline[MAXCHARS];

	for (i=0;i<6;i++) cte[i]=0.0;

	i=fseek(fp,-20000L,SEEK_END);
	if (i!=0) {
		fseek(fp,0L,SEEK_END);
		len=ftell(fp);
		if (len>8000) 
			i=fseek(fp,-8000L,SEEK_END);
		else
			i=fseek(fp,0L,SEEK_SET);
	}

	if (i!=0) {
		cte[0]=2100.0;
		return;
	}

	while (!feof(fp)) {
		oneline[0]='\0';
		fgets(oneline,sizeof(oneline),fp);

		if (rnxver<=2.99) {
			if (strlen(oneline)<26)	continue;

			if (strstr(oneline,"COMMENT")||strstr(oneline,"comment"))
				continue;

			if (oneline[3]!=' '||oneline[6]!=' '||oneline[9]!=' '||
				oneline[12]!=' '||oneline[15]!=' ')
				continue;

			if ((oneline[1]==' '&&oneline[2]==' ')||(oneline[4]==' '&&oneline[5]==' ')) 
				continue;

			cte[0]=str2num(oneline,1,2);
			if(cte[0]<=79.0&&cte[0]>=0.0) cte[0]+=2000.0;
			else cte[0]+=1900.0;
			cte[1]=str2num(oneline,4,2);
			cte[2]=str2num(oneline,7,2);
			cte[3]=str2num(oneline,10,2);
			cte[4]=str2num(oneline,13,2);
			cte[5]=str2num(oneline,15,11);
		}
		else {
			if (strlen(oneline)<29)	continue;

			if (oneline[0]!='>') continue;

			if (strstr(oneline,"COMMENT")||strstr(oneline,"comment"))
				continue;

			if (oneline[1]!=' '||oneline[6]!=' '||oneline[9]!=' '||
				oneline[12]!=' '||oneline[15]!=' ')
				continue;

			if (myRound(str2num(oneline,30,2))==4) {
				n=myRound(str2num(oneline,32,3));
				for (i=0;i<n;i++) {
					oneline[0]='\0';
					fgets(oneline,sizeof(oneline),fp);
				}
				continue;
			}
			cte[0]=str2num(oneline,2,4);
			cte[1]=str2num(oneline,7,2);
			cte[2]=str2num(oneline,10,2);
			cte[3]=str2num(oneline,13,2);
			cte[4]=str2num(oneline,16,2);
			cte[5]=str2num(oneline,19,10);
		}
	}
}

//read the  information of time tag of the first and last epoch, et al.
extern int getObsInfo(char ofilepath[], char anttype[], char rcvtype[], 
	double delta[3], gtime_t* ts, gtime_t* te, char *sitename, 
	char *filename, char *filename_ful, char *ext)
{
	FILE *fp=NULL;
	char *p,*q;
	char oneline[300]={'\0'};
	int ilen=strlen(ofilepath);
	double cts[6],cte[6],rnxver=2.11;

	if ((p=strrchr(ofilepath,FILEPATHSEP))) 
		xStrMid(PPP_Glo.obsDir,0,p-ofilepath,ofilepath);

	if ((p=strrchr(ofilepath,FILEPATHSEP))&&(q=strrchr(ofilepath,'.'))) {
		if (ext)	          {xStrMid(ext,q-ofilepath+1,ilen-(q-ofilepath+1),ofilepath);}
		if (sitename)     {strncpy(sitename,p+1,4); sitename[4]='\0';}
		if (filename)     {strncpy(filename,p+1,q-p-1);	filename[q-p-1]='\0';}
		if (filename_ful) {strncpy(filename_ful,p+1,ilen-(p-ofilepath));filename_ful[ilen-(p-ofilepath)]='\0';}
	}

	fp=fopen(ofilepath,"r");
	if (fp==NULL) {
		printf("*** ERROR: read file %s error!\n",ofilepath);
		return 0;
	}

	//get antenna type et al. from header of observation file
	while (fgets(oneline,MAXCHARS,fp)) {
		if(strstr(oneline,"ANT # / TYPE")&&anttype!=NULL)	
			xStrMid(anttype,20,20,oneline);
		if(strstr(oneline,"REC # / TYPE / VERS")&&rcvtype!=NULL)	
			xStrMid(rcvtype,20,20,oneline);
		if(strstr(oneline,"ANTENNA: DELTA H/E/N")&&delta!=NULL) {
			delta[0]=str2num(oneline,0,14);
			delta[1]=str2num(oneline,14,14);
			delta[2]=str2num(oneline,28,14);
		}
		if(strstr(oneline,"RINEX VERSION / TYPE") ) {
			rnxver=str2num(oneline,0,9);
		}
		if(strstr(oneline,"END OF HEADER")) 
			break;
	}

	getFirstTime(rnxver,fp,cts);
	getLastTime(rnxver,fp,cte);
	if (norm(cts,6)<=1.0||norm(cte,6)<=1.0) {
		printf("*** ERROR: norm(cts,6)<=1.0||norm(cte,6)<=1.0");
		return 0;
	}
	if (ts) *ts=epoch2time(cts);
	if (te) *te=epoch2time(cte);

	fclose(fp);

	return 1;
}

//find sinex file
extern int findSnxFile(const gtime_t ts, const gtime_t te, char dir[], char snxpaths[])
{
	int weeks,weeke;
	char tmp[MAXSTRPATH];
	char sep=(char)FILEPATHSEP;

	time2gpst(ts,&weeks);
	time2gpst(te,&weeke);

	if (weeks!=weeke) {
		printf("*** ERROR: obs data is across weeks!\n");
		return 0;
	}

	sprintf(tmp,"%s%cigs%4d.snx",dir,sep,weeks);
	
	if ((access(tmp,0))==-1) {  //for windows
	//if((access(tmp,F_OK))==-1) {  //for linux
		printf("*** ERROR: sinex file %s NOT found!\n",tmp);
		return 0;
	}
	else
		strcpy(snxpaths,tmp);

	return 1;
}

//get navigation file
extern int findNavFile(const gtime_t ts, const gtime_t te, char obspath[], 
	char *navpaths[], const int index, const int sys)
{
	gtime_t t;
	double secs,sece,ct[6];
	int weeks,weeke,days,daye,n,i,w,d,doy,index_=0,year;
	char tmp[MAXSTRPATH]={'\0'},dir[MAXSTRPATH]={'\0'},*p,file[MAXSTRPATH],c='n';

	strcpy(tmp,obspath);

	if (sys==SYS_GPS)      c=tmp[strlen(obspath)-1]='n';
	else if (sys==SYS_GLO) c=tmp[strlen(obspath)-1]='g';
	else if (sys==SYS_CMP) c=tmp[strlen(obspath)-1]='c';
	else if (sys==SYS_GAL) c=tmp[strlen(obspath)-1]='e';

	if (access(tmp,0)!=-1) {
		strcpy(navpaths[index],tmp);
		return 1;
	}

	if ((p=strrchr(obspath,FILEPATHSEP))) xStrMid(dir,0,p-obspath+1,obspath);

	secs=time2gpst(ts,&weeks); days=(int)(secs/86400.0);		
	sece=time2gpst(te,&weeke); daye=(int)(sece/86400.0);

	n=(weeke-weeks)*7+daye-days+1;

	if (n>30) {
		printf("The data is too long to be processed! (NAV>30)\n");
		return 0;
	}

	for (i=0;i<n;i++) {
		d=days+i;
		w=weeks;
		if (d>6) {
			w=w+(int)(d/7);
			d=d%7;
		}

		//calculate day of year
		t=gpst2time(w,d*86400.0);
		doy=(int)time2doy(t);
		time2epoch(t,ct);

		year=ct[0]>=2000.0?myRound((ct[0]-2000.0)):myRound((ct[0]-1900.0));

		num2str(doy,tmp,3);

		if (sys==SYS_GPS||sys==SYS_GAL) {
			sprintf(file,"%sbrdc%s0.%02d%c",dir,tmp,year,c);
			if (access(file,0)!=-1)
				strcpy(navpaths[index+index_++],file);
			else {
				sprintf(file,"%sauto%s0.%02d%c",dir,tmp,year,c);
				if (access(file,0)!=-1)
					strcpy(navpaths[index+index_++],file);
			}
		}
		else if (sys==SYS_GLO) {
			sprintf(file,"%sMCCK%s0.%02d%c",dir,tmp,year,c);
			if (access(file,0)!=-1)
				strcpy(navpaths[index+index_++],file);
			else {
				sprintf(file,"%sbrdc%s0.%02d%c",dir,tmp,year,c);
				if (access(file,0)!=-1)
					strcpy(navpaths[index+index_++],file);
			}
		}
		else if (sys==SYS_CMP) {
			sprintf(file,"%sauto%s0.%02d%c",dir,tmp,year,c);
			if (access(file,0)!=-1)
				strcpy(navpaths[index+index_++],file);
			else {
				sprintf(file,"%sbrdc%s0.%02d%c",dir,tmp,year,c);
				if (access(file,0)!=-1)
					strcpy(navpaths[index+index_++],file);
			}
		}
	}

	return index_;
}

//get navigation file (p file) for multi-GNSS
extern int findNavFile_p(const gtime_t ts, const gtime_t te, char obspath[], 
	char *navpaths[], const int index)
{
	int weeks,weeke,days,daye,n,i,w,d,doy,index_=0,year;
	double secs,sece,ct[6];
	char tmp[MAXSTRPATH]={'\0'},dir[MAXSTRPATH]={'\0'},*p,file[MAXSTRPATH],c;
	gtime_t t;

	strcpy(tmp,obspath);

	c=tmp[strlen(obspath)-1]='p';

	if (access(tmp,0)!=-1) {
		strcpy(navpaths[index],tmp);
		return 1;
	}

	if ((p=strrchr(obspath,FILEPATHSEP))) xStrMid(dir,0,p-obspath+1,obspath);

	secs=time2gpst(ts,&weeks); days=(int)(secs/86400.0);		
	sece=time2gpst(te,&weeke); daye=(int)(sece/86400.0);

	n=(weeke-weeks)*7+daye-days+1;

	if (n>30) {
		printf("The data is too long to be processed! (NAV>30)\n");
		return 0;
	}

	for (i=0;i<n;i++) {
		d=days+i; w=weeks;
		if (d>6) {
			w=w+(int)(d/7);
			d=d%7;
		}

		//calculate day of year
		t=gpst2time(w,d*86400.0);
		doy=(int)time2doy(t);
		time2epoch(t,ct);

		year=ct[0]>=2000.0?myRound((ct[0]-2000.0)):myRound((ct[0]-1900.0));

		num2str(doy,tmp,3);

		sprintf(file,"%sbrdc%s0.*%c",dir,tmp,c);
		if (access(file,0)!=-1)
			strcpy(navpaths[index+index_++],file);
		else {
			sprintf(file,"%sbrdm%s0.%02d%c",dir,tmp,year,c);
			if (access(file,0)!=-1)
				strcpy(navpaths[index+index_++],file);
		}
	}

	return index_;
}

//find precise clock files
extern int findClkFile(const gtime_t ts, const gtime_t te, char dir[], 
	char *clkpaths[], const int index)
{
	int weeks,weeke,days,daye,n,i,j,w,d,index_=0;
	double secs,sece;
	char tmp[200];
	char sep=(char)FILEPATHSEP;
	char ac[][4]={"igs","esa","cod","gfz","grg","jpl","gbm","com","grm","wum"};

	//secs=time2gpst(timeadd(ts,-301),&weeks); days=(int)(secs/86400.0);
	//sece=time2gpst(timeadd(te,+301),&weeke); daye=(int)(sece/86400.0);
	secs=time2gpst(ts,&weeks); days=(int)(secs/86400.0);
	sece=time2gpst(te,&weeke); daye=(int)(sece/86400.0);

	n=(weeke-weeks)*7+daye-days+1;

	if (n>30) {
		printf("The data is too long to be processed! (CLK>30)\n");
		return 0;
	}

	for (i=0;i<n;i++) {
		d=days+i; w=weeks;
		if (d>6) {
			w=w+(int)(d/7);
			d=d%7;
		}

		for (j=0;j<10;j++) {
			sprintf(tmp,"%s%c%s%04d%d.clk",dir,sep,ac[j],w,d);
			if (access(tmp,0)!=-1) {
				strcpy(clkpaths[index+index_++],tmp);
				break;
			}
		}
		if (j==10) {
			printf("*** ERROR: precise clock files NOT found!\n");
			return 0;
		}
	}

	return index_;
}

//find SP3 files
extern int findSp3File(const gtime_t ts, const gtime_t te, char dir[], 
	char *sp3paths[], const int index)
{
	int weeks,weeke,days,daye,n,i,j,k,w,d,index_=0;
	double secs,sece;
	char tmp[200],suffix[2][10]={"eph","sp3"};
	char sep=(char)FILEPATHSEP;
	char ac[][4]={"igs","esa","cod","gfz","grg","jpl","gbm","com","grm","wum"};
	gtime_t t;

	t=timeadd(ts,-900*10);
	secs=time2gpst(t,&weeks);
	t=timeadd(te,+900*10);
	sece=time2gpst(t,&weeke);
	days=(int)(secs/86400.0);
	daye=(int)(sece/86400.0);

	n=(weeke-weeks)*7+daye-days+1;

	if (n>30+2) {
		printf("The data is too long to be processed! (SP3>32)\n");
		return 0;
	}

	for (i=0;i<n;i++) {
		d=days+i; w=weeks;
		if (d>6) {
			w=w+(int)(d/7);
			d=d%7;
		}

		for (j=0;j<2;j++) {
			for (k=0;k<10;k++) {
				sprintf(tmp,"%s%c%s%04d%d.%s",dir,sep,ac[k],w,d,suffix[j]);
				if (access(tmp,0)!=-1) {
					strcpy(sp3paths[index+index_++],tmp);
					break;
				}
			}
			if (k<10) break;
		}
		if (j==2) {
			printf("*** ERROR: SP3 files NOT found!\n");
			return 0;
		}
	}

	return index_;
}

//find ION file
extern int findIonFile(const gtime_t ts, const gtime_t te, char dir[], char filepath[])
{
	double ct[6];
	int doy,year;
	char tmp[MAXSTRPATH]={'\0'},doy_str[MAXSTRPATH]={'\0'};
	char sep=(char)FILEPATHSEP;

	//calculate day of year
	doy=(int)time2doy(ts);
	time2epoch(ts,ct);
	year=ct[0]>=2000.0?myRound((ct[0]-2000.0)):myRound((ct[0]-1900.0));
	num2str(doy,doy_str,3);

	sprintf(tmp,"%s%cCODG%s0.%02dI",dir,sep,doy_str,year);
	if (access(tmp,0)!=-1) {
		strcpy(filepath,tmp);
	}
	else {
		printf("*** ERROR: ION file NOT found!\n");
		return 0;
	}

	return 1;
}

//find EOP file
extern int findErpFile(const gtime_t ts, const gtime_t te, char dir[], char filepath[])
{
	int j,weeks,weeke;
	char tmp[200];
	char sep=(char)FILEPATHSEP;
	char ac[][4]={"igs","esa","cod","gfz","grg","jpl","gbm","com","wum"};

	time2gpst(ts,&weeks);
	time2gpst(te,&weeke);

	if (weeks!=weeke) {
		return 0;
	}

	for (j=0;j<9;j++) {
		sprintf(tmp,"%s%c%s%04d7.erp",dir,sep,ac[j],weeke);
		if (access(tmp,0)!=-1) {
			strcpy(filepath,tmp);
			break;
		}
	}
	if (j==9) {
		printf("*** ERROR: EOP file NOT found!\n");
		return 0;
	}

	return 1;
}

//find P1-P2 DCB file
extern int findP1P2DcbFile(const gtime_t ts, const gtime_t te, char dir[], char dcbpath[])
{
	int bRes=0;
	char tmp[MAXSTRPATH];
	char sep=(char)FILEPATHSEP;
	double dcts[6],dcte[6];
	int i;
	gtime_t t0,t1,gt;
	double ct0[6],ct1[6];

	time2epoch(timeadd(ts,100),dcts);
	time2epoch(timeadd(te,-100),dcte);

	if (dcts[0]==dcte[0]&&dcts[1]==dcte[1]) {
		sprintf(tmp,"%s%cP1P2%02d%02d.DCB",dir,sep,(int)dcts[0]%100,(int)dcts[1]);
		if (access(tmp,0)!=-1) {
			strcpy(dcbpath,tmp);
			bRes=1;
		}
		else {
			if (dcts[1]==1.0) {dcts[1]=12.0; dcts[0]-=1.0;}
			else	 {dcts[1]-=1.0;}
			if (dcte[1]==12.0) {dcte[1]=1.0; dcte[0]+=1.0;}
			else	 {dcte[1]+=1.0;}
		}
	}

	if (bRes) return 1;

	if (dcts[0]!=dcte[0]||dcts[1]!=dcte[1]) {
		for (i=0;i<6;i++) { 
			ct0[i]=dcts[i];
			ct1[i]=dcte[i]; 
		}
		ct0[3]=ct0[4]=ct0[5]=ct1[3]=ct1[4]=ct1[5]=0.0;

		t0=epoch2time(ct0);
		t1=epoch2time(ct1);

		ct0[0]=ct0[1]=0;
		for (gt=t0;timediff(t1,gt)>=0;) {
			time2epoch(gt, ct1);
			if ((ct1[0]!=ct0[0])||(ct1[1]!=ct0[1])) {
				sprintf(tmp,"%s%cP1P2%02d%02d.DCB",dir,sep,(int)ct1[0]%100,(int)ct1[1]);
				if (access(tmp,0)!=-1) {
					strcpy(dcbpath,tmp);
					bRes=1;
					break;
				}
				ct0[0]=ct1[0];
				ct0[1]=ct1[1];
			}
			gt=timeadd(gt,25*86400);
		}
	}
	if (!bRes) printf("*** WARNING: P1-P2 DCB file NOT found!\n");

	return bRes;
}

//find P1-C1 DCB file
extern int findP1C1DcbFile(const gtime_t ts, const gtime_t te, char dir[], char dcbpath[])
{
	int bRes=0;
	char tmp[MAXSTRPATH];
	char sep=(char)FILEPATHSEP;
	double dcts[6],dcte[6];
	int i;
	gtime_t t0,t1,gt;
	double ct0[6],ct1[6];

	time2epoch(timeadd(ts,100),dcts);
	time2epoch(timeadd(te,-100),dcte);


	if (dcts[0]==dcte[0]&&dcts[1]==dcte[1]) {
		sprintf(tmp,"%s%cP1C1%02d%02d.DCB",dir,sep,(int)dcts[0]%100,(int)dcts[1]);
		if (access(tmp,0)!=-1) {
			strcpy(dcbpath,tmp);
			bRes=1;
		}
		else {
			if (dcts[1]==1.0) {dcts[1]=12.0; dcts[0]-=1.0;}
			else	 {dcts[1]-=1.0;}
			if (dcte[1]==12.0) {dcte[1]=1.0; dcte[0]+=1.0;}
			else	 {dcte[1]+=1.0;}
		}
	}

	if (bRes) return 1;

	if (dcts[0]!=dcte[0]||dcts[1]!=dcte[1]) {
		for (i=0;i<6;i++) { 
			ct0[i]=dcts[i];
			ct1[i]=dcte[i]; 
		}
		ct0[3]=ct0[4]=ct0[5]=ct1[3]=ct1[4]=ct1[5]=0.0;

		t0=epoch2time(ct0);
		t1=epoch2time(ct1);

		ct0[0]=ct0[1]=0;
		for (gt=t0;timediff(t1,gt)>=0;) {
			time2epoch(gt, ct1);
			if ((ct1[0]!=ct0[0])||(ct1[1]!=ct0[1])) {
				sprintf(tmp,"%s%cP1C1%02d%02d.DCB",dir,sep,(int)ct1[0]%100,(int)ct1[1]);
				if (access(tmp,0)!=-1) {
					strcpy(dcbpath,tmp);
					bRes=1;
					break;
				}
				ct0[0]=ct1[0];
				ct0[1]=ct1[1];
			}
			gt=timeadd(gt,25*86400);
		}
	}
	if (!bRes) printf("*** WARNING: P1-C1 DCB file NOT found!\n");

	return bRes;
}

//find P2-C2 DCB file
extern int findP2C2DcbFile(const gtime_t ts, const gtime_t te, char dir[], char dcbpath[])
{
	int bRes=0;
	char tmp[MAXSTRPATH];
	char sep=(char)FILEPATHSEP;
	double dcts[6],dcte[6];
	int i;
	gtime_t t0,t1,gt;
	double ct0[6],ct1[6];

	time2epoch(timeadd(ts,100),dcts);
	time2epoch(timeadd(te,-100),dcte);


	if (dcts[0]==dcte[0]&&dcts[1]==dcte[1]) {
		sprintf(tmp,"%s%cP2C2%02d%02d.DCB",dir,sep,(int)dcts[0]%100,(int)dcts[1]);
		if (access(tmp,0)!=-1) {
			strcpy(dcbpath,tmp);
			bRes=1;
		}
		else {
			if (dcts[1]==1.0) {dcts[1]=12.0; dcts[0]-=1.0;}
			else	 {dcts[1]-=1.0;}
			if (dcte[1]==12.0) {dcte[1]=1.0; dcte[0]+=1.0;}
			else	 {dcte[1]+=1.0;}
		}
	}

	if (bRes) return 1;

	if (dcts[0]!=dcte[0]||dcts[1]!=dcte[1]) {
		for (i=0;i<6;i++) { 
			ct0[i]=dcts[i];
			ct1[i]=dcte[i]; 
		}
		ct0[3]=ct0[4]=ct0[5]=ct1[3]=ct1[4]=ct1[5]=0.0;

		t0=epoch2time(ct0);
		t1=epoch2time(ct1);

		ct0[0]=ct0[1]=0;
		for (gt=t0;timediff(t1,gt)>=0;) {
			time2epoch(gt, ct1);
			if ((ct1[0]!=ct0[0])||(ct1[1]!=ct0[1])) {
				sprintf(tmp,"%s%cP2C2%02d%02d.DCB",dir,sep,(int)ct1[0]%100,(int)ct1[1]);
				if (access(tmp,0)!=-1) {
					strcpy(dcbpath,tmp);
					bRes=1;
					break;
				}
				ct0[0]=ct1[0];
				ct0[1]=ct1[1];
			}
			gt=timeadd(gt,25*86400);
		}
	}
	if (!bRes) printf("*** WARNING: P2-C2 DCB file NOT found!\n");

	return bRes;
}

//find MGEX DCB file
extern int findMGEXDcbFile(const gtime_t ts, const gtime_t te, char dir[], char dcbpath[])
{
	double ct[6];
	int doy,year;
	char tmp[MAXSTRPATH]={'\0'},doy_str[MAXSTRPATH]={'\0'};
	char sep=(char)FILEPATHSEP;

	//calculate DOY
	doy=(int)time2doy(ts);
	time2epoch(ts,ct);
	year=(int)ct[0];
	num2str(doy,doy_str,3);

	sprintf(tmp,"%s%cCAS0MGXRAP_%04d%s0000_01D_01D_DCB.BSX",dir,sep,year,doy_str);
	if (access(tmp,0)!=-1) {
		strcpy(dcbpath,tmp);
	}
	else {
		printf("*** ERROR: MGEX file %s NOT found!\n",tmp);
		return 0;
	}

	return 1;
}

//find BLQ file
extern int findBlqFile(const gtime_t ts, const gtime_t te, char dir[], char filepath[])
{
	char tmp[MAXSTRPATH];
	char sep=(char)FILEPATHSEP;

	sprintf(tmp,"%s%c%s",dir,sep,"ocnload.blq\0");
	if (access(tmp,0)!=-1) {
		strcpy(filepath,tmp);
		return 1;
	}
	else {
		printf("*** WARNING: BLQ file NOT found!\n");
		return 0;
	}
}

//find IGS antex file
extern int findAtxFile(const gtime_t ts, const gtime_t te, char dir[], char filepath[])
{
	char tmp0[500]={'\0'},tmp[500]={'\0'};
	char sep=(char)FILEPATHSEP;
	double ct_atx_1[6]={2006,11, 5,0,0,0.00000000};
	double ct_atx_2[6]={2011, 4,17,0,0,0.00000000};
	double ct_atx_3[6]={2017, 1,29,0,0,0.00000000};
	gtime_t gt_atx_1=epoch2time(ct_atx_1);
	gtime_t gt_atx_2=epoch2time(ct_atx_2);
	gtime_t gt_atx_3=epoch2time(ct_atx_3);
	double dt1=timediff(ts,gt_atx_1);
	double dt2=timediff(ts,gt_atx_2);
	double dt3=timediff(ts,gt_atx_3);

	if (dt1<0)      strcpy (tmp0,"igs01.atx\0");
	else if (dt2<0)	sprintf(tmp0,"igs05.atx%c",'\0');
	else	 if (dt3<0)	sprintf(tmp0,"igs08.atx%c",'\0');
	else			    sprintf(tmp0,"igs14.atx%c",'\0');

	sprintf(tmp,"%s%c%s",dir,sep,tmp0);
	if (access(tmp,0)!=-1) {
		strcpy(filepath,tmp);
		return 1;
	}
	else {
		printf("*** WARNING: IGS antex file NOT found!\n");
		return 0;
	}
}
