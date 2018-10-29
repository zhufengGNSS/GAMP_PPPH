#include "gamp.h"

#define MAX_N 100

int rcd[MAX_N];   //for record
int comb_j=0;

extern int myRound(const double dNum)
{
	int iNum;

	if (dNum>=0)
		iNum=(int)(dNum+0.5);
	else 
		iNum=(int)(dNum-0.5);

	return iNum;
}

extern void select_combination(const int l, const int p, const int n, const int m, int *sn)
{
	int i;

	if (l==m) {
		for (i=0;i<m;i++)
			sn[comb_j*m+i]=rcd[i];
		comb_j++;
		return;
	}

	for (i=p;i<=n-(m-l);i++) {
		rcd[l]=i+1;
		select_combination(l+1,i+1,n,m,sn);
	}
}

static double calStd_ex(const double *v, const int n, int *bused, double *ave)
{
	int i,j;
	double dave=0.0,std=0.0;

	for (i=j=0;i<n;i++) {
		if (bused) {
			if (bused[i]==0) continue;
		}

		dave+=v[i];
		j++;
	}

	if (j<=0) {
		return -1.0;
	}

	dave/=j;

	for (i=0;i<n;i++) {
		if (bused) {
			if (bused[i]==0) continue;
		}

		std+=(v[i]-dave)*(v[i]-dave);
	}

	std=sqrt(std/j);

	if (ave) *ave=dave;

	return std;
}

static int findGross_(const int nb, double *dv, const int nv, double *std_ex, double *ave_ex, int *ibadsn,
	                  const double ratio, const double minv, const double minstd)
{
	int bbad=0,*bused;
	int i,j,n,*it,*ibadsn_t;
	int j0,j9;
	double dstd_min=1.0e9,dt0=0.0,dt1=0.0,dave_min=0.0;

	if ((nv-nb<=nb)||nv<=0||ratio<=1.0) {
		if (std_ex)	*std_ex=-1.0;
		if (ave_ex)	*ave_ex=0.0;
		return 0;
	}

	bused=imat(nv,1);
	ibadsn_t=imat(nv,1);

	for (i=0,n=1;i<nb;i++)
		n=n*(nv-i)/(i+1);

	if (nb<=0) it=imat(n,1);
	else       it=imat(n*nb,1);

	comb_j=0;
	select_combination(0,0,nv,nb,it);

	for (i=0;i<n;i++) {
		j0=i*nb;
		j9=j0+nb;
		for (j=0;j<nv;j++)	bused[j]=1;
		for (j=j0;j<j9;j++)	bused[it[j]-1]=0;

		dt0=calStd_ex(dv,nv,bused,&dt1);

		if (dt0<dstd_min&&dt0>0.0) {
			dstd_min=dt0;
			dave_min=dt1;
		}

		if (dt0<minstd && dt0>0.0) break;
	}
	free(bused); free(it);

	for (i=j=0;i<nv;i++) {
		bbad=0;
		dt0=fabs(dv[i]-dave_min);
		if ( dt0>ratio*dstd_min&&dstd_min>1.0e-8) {
			if (minv>0.0) {
				if (fabs(dt0)>minv) bbad=1;
			}
			else
				bbad=1;
		}
		if (bbad) {
			if (j<nv)
				ibadsn_t[j]=i;
			else {
				//
			}
			j++;
		}
	}
	if (ibadsn) {
		for (i=0;i<j;i++) ibadsn[i]=ibadsn_t[i];
	}

	if (std_ex) *std_ex=dstd_min;
	if (ave_ex) *ave_ex=dave_min;

	free(ibadsn_t);

	return j;
}

extern int findGross(int ppp, int bMulGnss, double *v, const int nv, const int nbad, 
	                 double *std_ex, double *ave_ex, int *ibadsn, const double ratio, 
					 const double minv, const double stdmin)
{
	int i,j,badn=0,*ibadsn_t,badn_min=0;
	double dstd_min=1.0e9,dstd=0.0,dave=0.0,dave_min=0.0;
	int kk=4;

	if (nv<=1) return 0;

	ibadsn_t=imat(nv,1);

	if (bMulGnss) {
		kk--;
	}

	if (kk<=1) kk=1;

	for (i=0;i<=nbad;i++) {
		if (ppp&&(nv<i+kk||nv<2*i+1)&&i) continue;

		badn=findGross_(i,v,nv,&dstd,&dave,ibadsn_t,ratio,minv,stdmin);

		if (dstd>0.0&&dstd<dstd_min) {
			dstd_min=dstd;
			dave_min=dave;
			badn_min=badn;

			if (ibadsn)
				for (j=0;j<badn;j++) ibadsn[j]=ibadsn_t[j];
		}

		if (dstd>0.0&&dstd<=stdmin) break;
	}

	if (std_ex) *std_ex=dstd_min;
	if (ave_ex) *ave_ex=dave_min;

	free(ibadsn_t);

	return badn_min;
}

