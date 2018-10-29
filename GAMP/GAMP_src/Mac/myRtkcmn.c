#include "gamp.h"

/* constants -----------------------------------------------------------------*/
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define MAX(x,y)    ((x)>=(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))


static double sign(double d, double d1)
{
	if (d1>0)      return fabs(d);
	else if (d1<0) return 0-fabs(d);
	else	           return 0;
}

//get site coordinates from sinex file
extern void getCoord_snx(char *filepath, char *sitename, double coord[4])
{
	FILE *fp=NULL;
	char ch[100],chtmp[5];
	int inregion=-1;

	if (strlen(filepath)<3) return;

	coord[0]=coord[1]=coord[2]=coord[3]=0.0;

	if ((fp=fopen(filepath,"r"))==NULL) {
		printf("*** ERROR: sinex file %s open failed!\n",filepath);
		return;
	}

	while (!feof(fp)) {
		fgets(ch,sizeof(ch),fp);
		if (strstr(ch,"+SOLUTION/ESTIMATE")) {
			inregion=1;
		}
		if (strstr(ch,"-SOLUTION/ESTIMATE")) {
			inregion=0;
			break;
		}
		if (inregion==1) {
			strncpy(chtmp,&ch[14],4); chtmp[4]='\0';
			if (strcasecmp(chtmp,sitename)==0&&strstr(ch,"STAX")) {
				coord[0]=str2num(ch,47,21);
			}
			if (strcasecmp(chtmp,sitename)==0&&strstr(ch, "STAY")) {
				coord[1]=str2num(ch,47,21);
			}
			if (strcasecmp(chtmp,sitename)==0&&strstr(ch,"STAZ")) {
				coord[2]=str2num(ch,47,21);
				break;
			}		
		}
	}
	if (coord[0]*coord[1]*coord[2]!=0.0) coord[3]=2.0;   //来自snx文件

	fclose(fp);
}

//get site coordinates from site.crd
extern int getCoord_i(char *filepath, char *sitename, double coord[])
{
	FILE *fp=NULL;
	char ch[1000],chtmp[5];

	if (strlen(filepath)<3) return 0;

	coord[0]=coord[1]=coord[2]=0.0;

	if ((fp=fopen(filepath,"r"))==NULL) {
		printf("*** ERROR: site.crd file NOT found!");
		return 0;
	}

	while (!feof(fp)) {
		fgets(ch,sizeof(ch),fp);
		xStrMid(chtmp,0,4,ch);
		if (strcasecmp(chtmp,sitename)==0) {
			//coord[0]=str2num(ch, 104, 14);
			//coord[1]=str2num(ch, 119, 14);
			//coord[2]=str2num(ch, 134, 14);
			sscanf(ch,"%s %lf %lf %lf",chtmp,coord,coord+1,coord+2);
			break;
		}
	}
	fclose(fp);

	if (coord[0]*coord[1]*coord[2]==0.0) {
		printf("*** WARNING: Not found the coordinates for site %s!\n",sitename);
		return 0;
	}

	return 1;
}

//C     input data
//C     ----------
//C     dmjd: modified julian date
//C     pgt  station coordinate
//C     ellipsoidal latitude in radians
//C     longitude in radians
//C     ellipsoidal height in m
//C
//C     output data
//C     -----------
//C     pres: pressure in hPa
//C     temp: temperature in Celsius
//C     undu: Geoid undulation in m (from a 9x9 EGM based model)
//////////////////////////////////////////////////////////////////////////
// Given:                                                                        
//    DMJD           d      Modified Julian Date                                 
//    DLAT           d      Latitude given in radians (North Latitude)           
//    DLON           d      Longitude given in radians (East Longitude)          
//    DHGT           d      Ellipsoidal height in meters                         
//                                                                                
// Returned:                                                                     
//    PRES           d      Pressure given in hPa                                
//    TEMP           d      Temperature in degrees Celsius                       
//    UNDU           d      Geoid undulation in meters (Note 1) 

/*+                                                                               
*  - - - - - - - - -                                                             
*   G P T                                                                        
*  - - - - - - - - -                                                             
*                                                                                
*  This routine is part of the International Earth Rotation and                  
*  Reference Systems Service (IERS) Conventions software collection.             
*                                                                                
*  This subroutine determines Global Pressure and Temperature (Boehm et al. 2007)
*  based on Spherical Harmonics up to degree and order 9.                        
*                                                                                
*  In general, Class 1, 2, and 3 models represent physical effects that          
*  act on geodetic parameters while canonical models provide lower-level         
*  representations or basic computations that are used by Class 1, 2, or         
*  3 models.                                                                     
*                                                                                
*  Status: Class 1 model	                                                       
*                                                                                
*     Class 1 models are those recommended to be used a priori in the            
*     reduction of raw space geodetic data in order to determine                 
*     geodetic parameter estimates.                                              
*     Class 2 models are those that eliminate an observational                   
*     singularity and are purely conventional in nature.                         
*     Class 3 models are those that are not required as either Class             
*     1 or 2.                                                                    
*     Canonical models are accepted as is and cannot be classified as a          
*     Class 1, 2, or 3 model.                                                    
*                                                                                
*  Given:                                                                        
*     DMJD           d      Modified Julian Date                                 
*     DLAT           d      Latitude given in radians (North Latitude)           
*     DLON           d      Longitude given in radians (East Longitude)          
*     DHGT           d      Ellipsoidal height in meters                         
*                                                                                
*  Returned:                                                                     
*     PRES           d      Pressure given in hPa                                
*     TEMP           d      Temperature in degrees Celsius                       
*     UNDU           d      Geoid undulation in meters (Note 1)                  
*                                                                                
*  Notes:                                                                        
*                                                                                
*  1) This is from a 9x9 Earth Gravitational Model (EGM).                        
*                                                                                
*  Test case:                                                                    
*     given input: DMJD = 55055D0                                                
*                  DLAT = 0.6708665767D0 radians (NRAO, Green Bank, WV)          
*                  DLON = -1.393397187D0 radians                                 
*                  DHGT = 812.546 meters                                         
*     expected output: PRES = 918.0710638757363995D0 hPa                         
*                      TEMP = 19.31914181012882992D0 degrees Celsius             
*                      UNDU = -42.19185643717770517D0 meters                     
*                                                                                
*  References:                                                                   
*                                                                                
*     Boehm, J., Heinkelmann, R. and Schuh, H., 2007, "Short Note: A             
*     Global model of pressure and temperature for geodetic applications",       
*     Journal of Geodesy, 81(10), pp. 679-683.                                   
*                                                                                
*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),                   
*     IERS Technical Note No. 36, BKG (2010)
*                                                                                
*  Revisions:                                                                    
*  2006 June 12 J. Boehm    Original code                                        
*  2006 June 16 J. Boehm    Accounted for geoid undulation                       
*  2006 August 14 O. Montenbruck Recursions for Legendre polynomials             
*  2009 February 13 B.E. Stetzler Added header and copyright                     
*  2009 March 30 B.E. Stetzler More modifications and defined twopi              
*  2009 March 31 B.E. Stetzler Provided test case                                
*  2009 July  29 B.E. Stetzler Capitalized all variables for FORTRAN 77          
*                              compatibility and corrected test case             
*                              latitude and longitude coordinates                
*  2011 October 18 B.E. Stetzler Corrected test case output mentioned            
*                                in the header                                   
*-----------------------------------------------------------------------         
*/
//ftp://tai.bipm.org/iers/convupdt/chapter9/GPT.F
///*
extern void getGPT(const double *pos,double dmjd, double *pres, double *temp, double *undu)
{
	double V[10][10],W[10][10];
	int I,N,M,NMAX,MMAX;                                                  
	double DOY,TEMP0,PRES0,APM,APA,ATM,ATA,HORT,X,Y,Z;
	double DLAT,DLON,DHGT,PRES,TEMP,UNDU;
	const double TWOPI=6.283185307179586476925287;

	double a_geoid[55] = {
		-5.6195e-001,-6.0794e-002,-2.0125e-001,-6.4180e-002,-3.6997e-002,
		+1.0098e+001,+1.6436e+001,+1.4065e+001,+1.9881e+000,+6.4414e-001,
		-4.7482e+000,-3.2290e+000,+5.0652e-001,+3.8279e-001,-2.6646e-002,
		+1.7224e+000,-2.7970e-001,+6.8177e-001,-9.6658e-002,-1.5113e-002,
		+2.9206e-003,-3.4621e+000,-3.8198e-001,+3.2306e-002,+6.9915e-003,
		-2.3068e-003,-1.3548e-003,+4.7324e-006,+2.3527e+000,+1.2985e+000,
		+2.1232e-001,+2.2571e-002,-3.7855e-003,+2.9449e-005,-1.6265e-004,
		+1.1711e-007,+1.6732e+000,+1.9858e-001,+2.3975e-002,-9.0013e-004,
		-2.2475e-003,-3.3095e-005,-1.2040e-005,+2.2010e-006,-1.0083e-006,
		+8.6297e-001,+5.8231e-001,+2.0545e-002,-7.8110e-003,-1.4085e-004,
		-8.8459e-006,+5.7256e-006,-1.5068e-006,+4.0095e-007,-2.4185e-008 };
		double b_geoid[55]={
			+0.0000e+000,+0.0000e+000,-6.5993e-002,+0.0000e+000,+6.5364e-002,
			-5.8320e+000,+0.0000e+000,+1.6961e+000,-1.3557e+000,+1.2694e+000,
			0.0000e+000,-2.9310e+000,+9.4805e-001,-7.6243e-002,+4.1076e-002,
			+0.0000e+000,-5.1808e-001,-3.4583e-001,-4.3632e-002,+2.2101e-003,
			-1.0663e-002,+0.0000e+000,+1.0927e-001,-2.9463e-001,+1.4371e-003,
			-1.1452e-002,-2.8156e-003,-3.5330e-004,+0.0000e+000,+4.4049e-001,
			+5.5653e-002,-2.0396e-002,-1.7312e-003,+3.5805e-005,+7.2682e-005,
			+2.2535e-006,+0.0000e+000,+1.9502e-002,+2.7919e-002,-8.1812e-003,
			+4.4540e-004,+8.8663e-005,+5.5596e-005,+2.4826e-006,+1.0279e-006,
			+0.0000e+000,+6.0529e-002,-3.5824e-002,-5.1367e-003,+3.0119e-005,
			-2.9911e-005,+1.9844e-005,-1.2349e-006,-7.6756e-009,+5.0100e-008
		};
		double ap_mean[55]= {
			+1.0108e+003,+8.4886e+000,+1.4799e+000,-1.3897e+001,+3.7516e-003,
			-1.4936e-001,+1.2232e+001,-7.6615e-001,-6.7699e-002,+8.1002e-003,
			-1.5874e+001,+3.6614e-001,-6.7807e-002,-3.6309e-003,+5.9966e-004,
			+4.8163e+000,-3.7363e-001,-7.2071e-002,+1.9998e-003,-6.2385e-004,
			-3.7916e-004,+4.7609e+000,-3.9534e-001,+8.6667e-003,+1.1569e-002,
			+1.1441e-003,-1.4193e-004,-8.5723e-005,+6.5008e-001,-5.0889e-001,
			-1.5754e-002,-2.8305e-003,+5.7458e-004,+3.2577e-005,-9.6052e-006,
			-2.7974e-006,+1.3530e+000,-2.7271e-001,-3.0276e-004,+3.6286e-003,
			-2.0398e-004,+1.5846e-005,-7.7787e-006,+1.1210e-006,+9.9020e-008,
			+5.5046e-001,-2.7312e-001,+3.2532e-003,-2.4277e-003,+1.1596e-004,
			+2.6421e-007,-1.3263e-006,+2.7322e-007,+1.4058e-007,+4.9414e-009
		};
		double bp_mean[55]= {
			+0.0000e+000,+0.0000e+000,-1.2878e+000,+0.0000e+000,+7.0444e-001,
			+3.3222e-001,+0.0000e+000,-2.9636e-001,+7.2248e-003,+7.9655e-003,
			+0.0000e+000,+1.0854e+000,+1.1145e-002,-3.6513e-002,+3.1527e-003,
			+0.0000e+000,-4.8434e-001,+5.2023e-002,-1.3091e-002,+1.8515e-003,
			+1.5422e-004,+0.0000e+000,+6.8298e-001,+2.5261e-003,-9.9703e-004,
			-1.0829e-003,+1.7688e-004,-3.1418e-005,+0.0000e+000,-3.7018e-001,
			+4.3234e-002,+7.2559e-003,+3.1516e-004,+2.0024e-005,-8.0581e-006,
			-2.3653e-006,+0.0000e+000,+1.0298e-001,-1.5086e-002,+5.6186e-003,
			+3.2613e-005,+4.0567e-005,-1.3925e-006,-3.6219e-007,-2.0176e-008,
			+0.0000e+000,-1.8364e-001,+1.8508e-002,+7.5016e-004,-9.6139e-005,
			-3.1995e-006,+1.3868e-007,-1.9486e-007,+3.0165e-010,-6.4376e-010
		};
		double ap_amp[55]= {
			-1.0444e-001,+1.6618e-001,-6.3974e-002,+1.0922e+000,+5.7472e-001,
			-3.0277e-001,-3.5087e+000,+7.1264e-003,-1.4030e-001,+3.7050e-002,
			+4.0208e-001,-3.0431e-001,-1.3292e-001,+4.6746e-003,-1.5902e-004,
			+2.8624e+000,-3.9315e-001,-6.4371e-002,+1.6444e-002,-2.3403e-003,
			+4.2127e-005,+1.9945e+000,-6.0907e-001,-3.5386e-002,-1.0910e-003,
			-1.2799e-004,+4.0970e-005,+2.2131e-005,-5.3292e-001,-2.9765e-001,
			-3.2877e-002,+1.7691e-003,+5.9692e-005,+3.1725e-005,+2.0741e-005,
			-3.7622e-007,+2.6372e+000,-3.1165e-001,+1.6439e-002,+2.1633e-004,
			+1.7485e-004,+2.1587e-005,+6.1064e-006,-1.3755e-008,-7.8748e-008,
			-5.9152e-001,-1.7676e-001,+8.1807e-003,+1.0445e-003,+2.3432e-004,
			+9.3421e-006,+2.8104e-006,-1.5788e-007,-3.0648e-008,+2.6421e-010
		};
		double bp_amp[55]= {
			+0.0000e+000,+0.0000e+000,+9.3340e-001,+0.0000e+000,+8.2346e-001,
			+2.2082e-001,+0.0000e+000,+9.6177e-001,-1.5650e-002,+1.2708e-003,
			+0.0000e+000,-3.9913e-001,+2.8020e-002,+2.8334e-002,+8.5980e-004,
			+0.0000e+000,+3.0545e-001,-2.1691e-002,+6.4067e-004,-3.6528e-005,
			-1.1166e-004,+0.0000e+000,-7.6974e-002,-1.8986e-002,+5.6896e-003,
			-2.4159e-004,-2.3033e-004,-9.6783e-006,+0.0000e+000,-1.0218e-001,
			-1.3916e-002,-4.1025e-003,-5.1340e-005,-7.0114e-005,-3.3152e-007,
			+1.6901e-006,+0.0000e+000,-1.2422e-002,+2.5072e-003,+1.1205e-003,
			-1.3034e-004,-2.3971e-005,-2.6622e-006,+5.7852e-007,+4.5847e-008,
			+0.0000e+000,+4.4777e-002,-3.0421e-003,+2.6062e-005,-7.2421e-005,
			+1.9119e-006,+3.9236e-007,+2.2390e-007,+2.9765e-009,-4.6452e-009
		};
		double at_mean[55]= {
			+1.6257e+001,+2.1224e+000,+9.2569e-001,-2.5974e+001,+1.4510e+000,
			+9.2468e-002,-5.3192e-001,+2.1094e-001,-6.9210e-002,-3.4060e-002,
			-4.6569e+000,+2.6385e-001,-3.6093e-002,+1.0198e-002,-1.8783e-003,
			+7.4983e-001,+1.1741e-001,+3.9940e-002,+5.1348e-003,+5.9111e-003,
			+8.6133e-006,+6.3057e-001,+1.5203e-001,+3.9702e-002,+4.6334e-003,
			+2.4406e-004,+1.5189e-004,+1.9581e-007,+5.4414e-001,+3.5722e-001,
			+5.2763e-002,+4.1147e-003,-2.7239e-004,-5.9957e-005,+1.6394e-006,
			-7.3045e-007,-2.9394e+000,+5.5579e-002,+1.8852e-002,+3.4272e-003,
			-2.3193e-005,-2.9349e-005,+3.6397e-007,+2.0490e-006,-6.4719e-008,
			-5.2225e-001,+2.0799e-001,+1.3477e-003,+3.1613e-004,-2.2285e-004,
			-1.8137e-005,-1.5177e-007,+6.1343e-007,+7.8566e-008,+1.0749e-009
		};
		double bt_mean[55]= {
			+0.0000e+000,+0.0000e+000,+1.0210e+000,+0.0000e+000,+6.0194e-001,
			+1.2292e-001,+0.0000e+000,-4.2184e-001,+1.8230e-001,+4.2329e-002,
			+0.0000e+000,+9.3312e-002,+9.5346e-002,-1.9724e-003,+5.8776e-003,
			+0.0000e+000,-2.0940e-001,+3.4199e-002,-5.7672e-003,-2.1590e-003,
			+5.6815e-004,+0.0000e+000,+2.2858e-001,+1.2283e-002,-9.3679e-003,
			-1.4233e-003,-1.5962e-004,+4.0160e-005,+0.0000e+000,+3.6353e-002,
			-9.4263e-004,-3.6762e-003,+5.8608e-005,-2.6391e-005,+3.2095e-006,
			-1.1605e-006,+0.0000e+000,+1.6306e-001,+1.3293e-002,-1.1395e-003,
			+5.1097e-005,+3.3977e-005,+7.6449e-006,-1.7602e-007,-7.6558e-008,
			+0.0000e+000,-4.5415e-002,-1.8027e-002,+3.6561e-004,-1.1274e-004,
			+1.3047e-005,+2.0001e-006,-1.5152e-007,-2.7807e-008,+7.7491e-009
		};
		double at_amp[55]= {
			-1.8654e+000,-9.0041e+000,-1.2974e-001,-3.6053e+000,+2.0284e-002,
			+2.1872e-001,-1.3015e+000,+4.0355e-001,+2.2216e-001,-4.0605e-003,
			+1.9623e+000,+4.2887e-001,+2.1437e-001,-1.0061e-002,-1.1368e-003,
			-6.9235e-002,+5.6758e-001,+1.1917e-001,-7.0765e-003,+3.0017e-004,
			+3.0601e-004,+1.6559e+000,+2.0722e-001,+6.0013e-002,+1.7023e-004,
			-9.2424e-004,+1.1269e-005,-6.9911e-006,-2.0886e+000,-6.7879e-002,
			-8.5922e-004,-1.6087e-003,-4.5549e-005,+3.3178e-005,-6.1715e-006,
			-1.4446e-006,-3.7210e-001,+1.5775e-001,-1.7827e-003,-4.4396e-004,
			+2.2844e-004,-1.1215e-005,-2.1120e-006,-9.6421e-007,-1.4170e-008,
			+7.8720e-001,-4.4238e-002,-1.5120e-003,-9.4119e-004,+4.0645e-006,
			-4.9253e-006,-1.8656e-006,-4.0736e-007,-4.9594e-008,+1.6134e-009
		};
		double bt_amp[55]= {
			+0.0000e+000,+0.0000e+000,-8.9895e-001,+0.0000e+000,-1.0790e+000,
			-1.2699e-001,+0.0000e+000,-5.9033e-001,+3.4865e-002,-3.2614e-002,
			+0.0000e+000,-2.4310e-002,+1.5607e-002,-2.9833e-002,-5.9048e-003,
			+0.0000e+000,+2.8383e-001,+4.0509e-002,-1.8834e-002,-1.2654e-003,
			-1.3794e-004,+0.0000e+000,+1.3306e-001,+3.4960e-002,-3.6799e-003,
			-3.5626e-004,+1.4814e-004,+3.7932e-006,+0.0000e+000,+2.0801e-001,
			+6.5640e-003,-3.4893e-003,-2.7395e-004,+7.4296e-005,-7.9927e-006,
			-1.0277e-006,+0.0000e+000,+3.6515e-002,-7.4319e-003,-6.2873e-004,
			-8.2461e-005,+3.1095e-005,-5.3860e-007,-1.2055e-007,-1.1517e-007,
			+0.0000e+000,+3.1404e-002,+1.5580e-002,-1.1428e-003,+3.3529e-005,
			+1.0387e-005,-1.9378e-006,-2.7327e-007,+7.5833e-009,-9.2323e-009
		};

		DLAT=pos[0];
		DLON=pos[1];
		DHGT=pos[2];

//      Reference day is 28 January 1980                                           
//      This is taken from Niell (1996) to be consistent (See References)          
//      For constant values use: doy = 91.3125
		DOY=dmjd-44239+1-28;

//		Define degree n and order m EGM  
		NMAX=MMAX=9;

//      Define unit vector
		X=cos(DLAT)*cos(DLON);
		Y=cos(DLAT)*sin(DLON);
		Z=sin(DLAT);

//      Legendre polynomials                                                       
		V[1-1][1-1] = 1.0;                                                             
		W[1-1][1-1] = 0.0;                                                            
		V[2-1][1-1] = Z * V[1-1][1-1];
		W[2-1][1-1] = 0.0;

		for (N=2;N<=NMAX;N++) {
			V[N+1-1][1-1] = ((2*N-1) * Z * V[N-1][1-1] - (N-1) * V[N-1-1][1-1]) / N ;             
			W[N+1-1][1-1] = 0.0;    
		}

		for( M=1;M<=NMAX;M++) {
			V[M+1-1][M+1-1] = (2*M-1) * (X*V[M-1][M-1] - Y*W[M-1][M-1]);                         
			W[M+1-1][M+1-1] = (2*M-1) * (X*W[M-1][M-1] + Y*V[M-1][M-1]);                            
			if (M < NMAX) {
				V[M+2-1][M+1-1] = (2*M+1) * Z * V[M+1-1][M+1-1];                          
				W[M+2-1][M+1-1] = (2*M+1) * Z * W[M+1-1][M+1-1];
			}                                                                 
			for (N=M+2;N<=NMAX;N++) {
				V[N+1-1][M+1-1] = ((2*N-1)*Z*V[N-1][M+1-1] - (N+M-1)*V[N-1-1][M+1-1]) / (N-M);
				W[N+1-1][M+1-1] = ((2*N-1)*Z*W[N-1][M+1-1] - (N+M-1)*W[N-1-1][M+1-1]) / (N-M);
			}
		}

//      Geoidal height                                                             
		UNDU = 0.0;                                                                
		I = 0;

		for (N=0;N<=NMAX;N++) {
			for (M=0;M<=N;M++) {
				I = I+1;
				UNDU = UNDU + (a_geoid[I-1]*V[N+1-1][M+1-1] + b_geoid[I-1]*W[N+1-1][M+1-1]);
			}
		}

//      orthometric height                                                         
		HORT = DHGT - UNDU;

//      Surface pressure on the geoid                                              
		APM = 0.0;
		APA = 0.0;
		I = 0;


		for (N=0;N<=NMAX;N++) {
			for (M=0;M<=N;M++) {
				I = I+1;
				APM = APM + ( ap_mean[I-1]*V[N+1-1][M+1-1] + bp_mean[I-1]*W[N+1-1][M+1-1] ) ;
				APA = APA + ( ap_amp[I-1] *V[N+1-1][M+1-1] + bp_amp[I-1] *W[N+1-1][M+1-1] ) ;           
			}
		}
		
		PRES0  = APM + APA*cos(DOY/365.25*TWOPI);

//      height correction for pressure                                             
		PRES = PRES0*pow(1.0-0.0000226*HORT, 5.225);

//      Surface temperature on the geoid                                           
		ATM = 0.0;                                                              
		ATA = 0.0;                                                                 
		I = 0;

		for (N=0;N<=NMAX;N++) {
			for (M=0;M<=N;M++) {                                               
				I = I+1;
				ATM = ATM + (at_mean[I-1]*V[N+1-1][M+1-1] + bt_mean[I-1]*W[N+1-1][M+1-1]);
				ATA = ATA + (at_amp[I-1] *V[N+1-1][M+1-1] + bt_amp[I-1] *W[N+1-1][M+1-1]);
			}
		}

		TEMP0 =  ATM + ATA*cos(DOY/365.25*TWOPI);                                

//      height correction for temperature 
		TEMP = TEMP0 - 0.0065*HORT;


		*pres=PRES;
		*temp=TEMP;
		*undu=UNDU;

//  Finished.                                                                      
}

//ftp://tai.bipm.org/iers/convupdt/chapter9/GMF.F
/*  GMF */
//    This subroutine determines the Global Mapping Functions GMF
//	Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006), 
//  Global Mapping Functions (GMF): A new empirical mapping function based on numerical weather model data,
//  Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
//    input data
//    ----------
//    dmjd: modified julian date
//    dlat: latitude in radians
//    dlon: longitude in radians
//    dhgt: height in m
//    zd:   zenith distance in radians

//    output data
//    -----------
//    gmfh: hydrostatic mapping function
//    gmfw: wet mapping function

//    Johannes Boehm, 2005 August 30
extern void tropmapf_gmf(gtime_t gt, const double blh[3], double elev, double *gmfh, double *gmfw)
{
	double dfac[20],P[10][10],aP[55],bP[55],t,dmjd;
	int i,j,ir,k,n,m,nmax,mmax;
	double doy,phh;
	double ah,bh,ch,aw,bw,cw;
	double ahm,aha,awm,awa;
	double c10h,c11h,c0h;
	double a_ht,b_ht,c_ht;
	double sine,beta,gamma,topcon;
	double hs_km,ht_corr,ht_corr_coef;
	double dlat,dlon,dhgt,sum,dt;

	//double PI = 3.1415926535897932384626433;
	double TWOPI = 6.283185307179586476925287;

	static double ah_mean[55] = {
		+1.2517e+02,	+8.503e-01,	+6.936e-02,	-6.760e+00, +1.771e-01,
		+1.130e-02,		+5.963e-01,	+1.808e-02, +2.801e-03, -1.414e-03,
		-1.212e+00,		+9.300e-02,	+3.683e-03, +1.095e-03, +4.671e-05,
		+3.959e-01,		-3.867e-02, +5.413e-03, -5.289e-04, +3.229e-04,
		+2.067e-05,		+3.000e-01, +2.031e-02, +5.900e-03, +4.573e-04,
		-7.619e-05,		+2.327e-06, +3.845e-06, +1.182e-01, +1.158e-02,
		+5.445e-03,		+6.219e-05, +4.204e-06, -2.093e-06, +1.540e-07,
		-4.280e-08,		-4.751e-01, -3.490e-02, +1.758e-03, +4.019e-04,
		-2.799e-06,		-1.287e-06, +5.468e-07, +7.580e-08, -6.300e-09,
		-1.160e-01,		+8.301e-03, +8.771e-04, +9.955e-05, -1.718e-06,
		-2.012e-06,		+1.170e-08, +1.790e-08, -1.300e-09, +1.000e-10
	};

	static double bh_mean[55] = {
		+0.000e+00,	+0.000e+00, +3.249e-02, +0.000e+00, +3.324e-02,
		+1.850e-02, +0.000e+00, -1.115e-01, +2.519e-02, +4.923e-03,
		+0.000e+00, +2.737e-02, +1.595e-02, -7.332e-04, +1.933e-04,
		+0.000e+00, -4.796e-02, +6.381e-03, -1.599e-04, -3.685e-04,
		+1.815e-05, +0.000e+00, +7.033e-02, +2.426e-03, -1.111e-03,
		-1.357e-04, -7.828e-06, +2.547e-06, +0.000e+00, +5.779e-03,
		+3.133e-03, -5.312e-04, -2.028e-05, +2.323e-07, -9.100e-08,
		-1.650e-08, +0.000e+00, +3.688e-02, -8.638e-04, -8.514e-05,
		-2.828e-05, +5.403e-07, +4.390e-07, +1.350e-08, +1.800e-09,
		+0.000e+00, -2.736e-02, -2.977e-04, +8.113e-05, +2.329e-07,
		+8.451e-07, +4.490e-08, -8.100e-09, -1.500e-09, +2.000e-10
	};

	static double ah_amp[55] = {
		-2.738e-01, -2.837e+00, +1.298e-02, -3.588e-01, +2.413e-02,
		+3.427e-02, -7.624e-01, +7.272e-02, +2.160e-02, -3.385e-03,
		+4.424e-01, +3.722e-02, +2.195e-02, -1.503e-03, +2.426e-04,
		+3.013e-01, +5.762e-02, +1.019e-02, -4.476e-04, +6.790e-05,
		+3.227e-05, +3.123e-01, -3.535e-02, +4.840e-03, +3.025e-06,
		-4.363e-05, +2.854e-07, -1.286e-06, -6.725e-01, -3.730e-02,
		+8.964e-04, +1.399e-04, -3.990e-06, +7.431e-06, -2.796e-07,
		-1.601e-07, +4.068e-02, -1.352e-02, +7.282e-04, +9.594e-05,
		+2.070e-06, -9.620e-08, -2.742e-07, -6.370e-08, -6.300e-09,
		+8.625e-02, -5.971e-03, +4.705e-04, +2.335e-05, +4.226e-06,
		+2.475e-07, -8.850e-08, -3.600e-08, -2.900e-09, +0.000e+00
	};

	static double bh_amp[55] = {
		+0.000e+00, +0.000e+00, -1.136e-01, +0.000e+00, -1.868e-01,
		-1.399e-02, +0.000e+00, -1.043e-01, +1.175e-02, -2.240e-03,
		+0.000e+00, -3.222e-02, +1.333e-02, -2.647e-03, -2.316e-05,
		+0.000e+00, +5.339e-02, +1.107e-02, -3.116e-03, -1.079e-04,
		-1.299e-05, +0.000e+00, +4.861e-03, +8.891e-03, -6.448e-04,
		-1.279e-05, +6.358e-06, -1.417e-07, +0.000e+00, +3.041e-02,
		+1.150e-03, -8.743e-04, -2.781e-05, +6.367e-07, -1.140e-08,
		-4.200e-08, +0.000e+00, -2.982e-02, -3.000e-03, +1.394e-05,
		-3.290e-05, -1.705e-07, +7.440e-08, +2.720e-08, -6.600e-09,
		+0.000e+00, +1.236e-02, -9.981e-04, -3.792e-05, -1.355e-05,
		+1.162e-06, -1.789e-07, +1.470e-08, -2.400e-09, -4.000e-10
	};

	static double aw_mean[55] = {
		+5.640e+01, +1.555e+00, -1.011e+00, -3.975e+00, +3.171e-02,
		+1.065e-01, +6.175e-01, +1.376e-01, +4.229e-02, +3.028e-03,
		+1.688e+00, -1.692e-01, +5.478e-02, +2.473e-02, +6.059e-04,
		+2.278e+00, +6.614e-03, -3.505e-04, -6.697e-03, +8.402e-04,
		+7.033e-04, -3.236e+00, +2.184e-01, -4.611e-02, -1.613e-02,
		-1.604e-03, +5.420e-05, +7.922e-05, -2.711e-01, -4.406e-01,
		-3.376e-02, -2.801e-03, -4.090e-04, -2.056e-05, +6.894e-06,
		+2.317e-06, +1.941e+00, -2.562e-01, +1.598e-02, +5.449e-03,
		+3.544e-04, +1.148e-05, +7.503e-06, -5.667e-07, -3.660e-08,
		+8.683e-01, -5.931e-02, -1.864e-03, -1.277e-04, +2.029e-04,
		+1.269e-05, +1.629e-06, +9.660e-08, -1.015e-07, -5.000e-10
	};

	static double bw_mean[55] = {
		+0.000e+00, +0.000e+00, +2.592e-01, +0.000e+00, +2.974e-02,
		-5.471e-01, +0.000e+00, -5.926e-01, -1.030e-01, -1.567e-02,
		+0.000e+00, +1.710e-01, +9.025e-02, +2.689e-02, +2.243e-03,
		+0.000e+00, +3.439e-01, +2.402e-02, +5.410e-03, +1.601e-03,
		+9.669e-05, +0.000e+00, +9.502e-02, -3.063e-02, -1.055e-03,
		-1.067e-04, -1.130e-04, +2.124e-05, +0.000e+00, -3.129e-01,
		+8.463e-03, +2.253e-04, +7.413e-05, -9.376e-05, -1.606e-06,
		+2.060e-06, +0.000e+00, +2.739e-01, +1.167e-03, -2.246e-05,
		-1.287e-04, -2.438e-05, -7.561e-07, +1.158e-06, +4.950e-08,
		+0.000e+00, -1.344e-01, +5.342e-03, +3.775e-04, -6.756e-05,
		-1.686e-06, -1.184e-06, +2.768e-07, +2.730e-08, +5.700e-09
	};

	static double aw_amp[55] = {
		+1.023e-01, -2.695e+00, +3.417e-01, -1.405e-01, +3.175e-01,
		+2.116e-01, +3.536e+00, -1.505e-01, -1.660e-02, +2.967e-02,
		+3.819e-01, -1.695e-01, -7.444e-02, +7.409e-03, -6.262e-03,
		-1.836e+00, -1.759e-02, -6.256e-02, -2.371e-03, +7.947e-04,
		+1.501e-04, -8.603e-01, -1.360e-01, -3.629e-02, -3.706e-03,
		-2.976e-04, +1.857e-05, +3.021e-05, +2.248e+00, -1.178e-01,
		+1.255e-02, +1.134e-03, -2.161e-04, -5.817e-06, +8.836e-07,
		-1.769e-07, +7.313e-01, -1.188e-01, +1.145e-02, +1.011e-03,
		+1.083e-04, +2.570e-06, -2.140e-06, -5.710e-08, +2.000e-08,
		-1.632e+00, -6.948e-03, -3.893e-03, +8.592e-04, +7.577e-05,
		+4.539e-06, -3.852e-07, -2.213e-07, -1.370e-08, +5.800e-09
	};

	static double bw_amp[55] = {
		+0.000e+00, +0.000e+00, -8.865e-02, +0.000e+00, -4.309e-01,
		+6.340e-02, +0.000e+00, +1.162e-01, +6.176e-02, -4.234e-03,
		+0.000e+00, +2.530e-01, +4.017e-02, -6.204e-03, +4.977e-03,
		+0.000e+00, -1.737e-01, -5.638e-03, +1.488e-04, +4.857e-04,
		-1.809e-04, +0.000e+00, -1.514e-01, -1.685e-02, +5.333e-03,
		-7.611e-05, +2.394e-05, +8.195e-06, +0.000e+00, +9.326e-02,
		-1.275e-02, -3.071e-04, +5.374e-05, -3.391e-05, -7.436e-06,
		+6.747e-07, +0.000e+00, -8.637e-02, -3.807e-03, -6.833e-04,
		-3.861e-05, -2.268e-05, +1.454e-06, +3.860e-07, -1.068e-07,
		+0.000e+00, -2.658e-02, -1.947e-03, +7.131e-04, -3.506e-05,
		+1.885e-07, +5.792e-07, +3.990e-08, +2.000e-08, -5.700e-09
	};

	mjd_t mjd;

	time2mjd(gt,&mjd);
	dmjd=mjd.day+(mjd.ds.sn+mjd.ds.tos)/86400.0;

	//reference day is 28 January
	//this is taken from Niell(1996) to be consistent
	//doy=dmjd-44239.0+1-28;
	doy=dmjd-44239.0-27;

	dlat=blh[0]; dlon=blh[1]; dhgt=blh[2];
	//parameter t
	t=sin(dlat);

	//degree n and order m
	nmax=9;
	mmax=9;

	//determine nmax!(faktorielle) moved by 1
	dfac[0]=1;
	for (i=1;i<=2*nmax+1;i++) 
		dfac[i]=dfac[i-1]*i;

	//determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
	for (i=0;i<=nmax;i++) {
		for (j=0;j<=MIN(i,mmax);j++) {
			ir=(int)((i-j)/2);
			sum=0.0;

			for (k=0;k<=ir;k++) {
				sum=sum+pow(-1.0,k)*dfac[2*i-2*k]/dfac[k]/dfac[i-k]/dfac[i-j-2*k]*pow(t,i-j-2*k);
			}

			//Legender functions moved by 1
			P[i][j]=1.0/pow(2.0,i)*sqrt(pow(1-t*t,j))*sum;
		}
	}

	//spherical harmonics
	i=0;
	for (n=0;n<=9;n++) {
		for (m=0;m<=n;m++) {
			i=i+1;
			dt=m*dlon;
			aP[i-1]=P[n][m]*cos(dt);
			bP[i-1]=P[n][m]*sin(dt);
		}
	}

	//hydrostatic
	bh=0.0029;
	c0h=0.062;

	if (dlat<0.0) {	//southern hemisphere
		phh=PI;
		c11h=0.007;
		c10h=0.002;
	}
	else {	//northern hemisphere
		phh=0;
		c11h=0.005;
		c10h=0.001;
	}

	ch=c0h+((cos(doy/365.25*TWOPI+phh)+1.0)*c11h/2.0+c10h)*(1.0-cos(dlat));

	ahm=0.0;
	aha=0.0;
	for (i=1;i<=55;i++) {
		ahm=ahm+(ah_mean[i-1]*aP[i-1]+bh_mean[i-1]*bP[i-1])*1.0e-5;
		aha=aha+(ah_amp[i-1] *aP[i-1]+bh_amp[i-1] *bP[i-1])*1.0e-5;
	}

	ah=ahm+aha*cos(doy/365.25*TWOPI);

	sine=sin(elev);
	beta=bh/(sine+ch);
	gamma=ah/(sine+beta);
	topcon=(1.0+ah/(1.0+bh/(1.0+ch)));
	*gmfh=topcon/(sine+gamma);

	//height correction for hydrostatic mapping function from Niell (1996)
	a_ht=2.53e-5;	//2.53 from http://maia.usno.navy.mil/conv2010/chapter9/GMF.F

	b_ht=5.49e-3;
	c_ht=1.14e-3;
	hs_km=dhgt/1000.0;

	beta=b_ht/(sine+c_ht);
	gamma=a_ht/(sine+beta);
	topcon=(1.0+a_ht/(1.0+b_ht/(1.0+c_ht)));
	ht_corr_coef=1.0/sine-topcon/(sine+gamma);
	ht_corr=ht_corr_coef*hs_km;
	*gmfh=*gmfh+ht_corr;

	//wet
	bw=0.00146;
	cw=0.04391;

	awm=0.0;
	awa=0.0;
	for (i=1;i<=55;i++) {
		awm=awm+(aw_mean[i-1]*aP[i-1]+bw_mean[i-1]*bP[i-1])*1e-5;
		awa=awa+(aw_amp[i-1] *aP[i-1]+bw_amp[i-1] *bP[i-1])*1e-5;
	}
	aw=awm+awa*cos(doy/365.25*TWOPI);

	beta=bw/(sine+cw);
	gamma=aw/(sine+beta);
	topcon=(1.0+aw/(1.0+bw/(1.0+cw)));
	*gmfw=topcon/(sine+gamma);
}

/*----------------------------------------------------------------------------
C
C       NAME	        ECLIPS (version Sep  2011)
C
C	PURPOSE 	DETECT ECLIPSING & YAW ROTATE ECLIP. SATELLITES
C                       (THE INPUT BODY-X UNIT VECTOR - SANTXYZ IS YAW-ROTATED
C                        BY PHI-YANGLE (THE ECL-NOMINAL) YAW ANGLE DIFFERENCE)  
C
C       COPYRIGHT       GEODETIC SURVEY DIVISION, 2011.
C                       ALL RIGHTS RESERVED.
C                       ALL TERMS AND CONDITIONS APPLY AS DETAILED IN
C                       " TERMS AND CONDITIONS FOR SOFTWARE " 
C
C       CONTACT         kouba@geod.nrcan.gc.ca
C
C       UPDATE HISTORY: Aug. 23, 2011:1996-2008 W. AVERAGES of JPL REPROCESSING
C                                YRATE SOLUTIONS FOR ALL II/IIA CODED IN DATA
C                                STATEMENT, THIS ENABLES REPROCESSING FROM 1996
C                       Sep 26, 2011: 1. Corrected bug causing Block IIF shadow 
C                               CROSSING MANEVURE WITH 0.06 DEG/SEC RATE EVEN 
C                               FOR BETADG > 8 DEG
C                                     2. CORRECTED/improved IIA RECOVERY logic
C
C	PARAMETERS	DESCRIPTION
C
C        IDIR		DIRECTION OF PROCESSING (1=FORWARD, -1=BACKWARD)
C        IPRN           SV PRN NUMBER (.LE.32 FOR GPS, .GT.32 FOR GLONASS)
C        TTAG           OBSERVATION EPOCH TIME TAG (EG SEC OF GPS WEEK)
C        SVBCOS         SV BETA ANGLE (BETWEEN SV AND SUN RADIUS) COSINE			对应论文中的E
C        ANOON          SV BETA ANGLE LIMIT (DEG) FOR A NOON TURN MANEVURE 
C        ANIGHT         SV BETA ANGLE LIMIT (DEG) FOR A NIGHT SHADOW CROSSING
C        NECLIPS        NUMBER OF ECLIPSING FOR THE PRN SATELLITE 
C        ECLSTM         ECLIPSE START TIME(EG SEC OF GPS WEEK)
C        ECLETM         ECLIPSE END TIME  ( "         "      )
C        IECLIPS        SV ECLIPSING (0=NO,1, 2=YES(1=night, 2=noon))
C        PI             = PI=3.1415926536D0
C        XSV(3)         SV X, Y, Z (m)(ITRF)
C        SANTXYZ        BODY-X UNIT VECTOR (ITRF)
C                       WARNING: THE IIA BODY-X ORIENTATION EXPECTED FOR THE IIR
C                        THE  BODY-X REVERSED FOR IIR (SEE BELOW) & RETURNED
C        VSVC           SV INERTIAL VELOCIRY VECTOR IN ITRF
C        BETA           90 DEG + THE SUN ANGLE WITH ORBITAL PLANE(IN RAD)
C        IBLK           SV BLOCK  1=I, 2=II, 3=IIA, IIR=(4, 5) IIF=6
C
C        INTERNAL PARAMETRS DESCRIPTION
C        YANGLE         THE NOMINAL YAW ANGLE
C        PHI            THE ECLIPSING YAW ANGLE            
C
C *********************************************************************
*/
static int eclips_(int IPRN, double SVBCOS, double ANIGHT, double BETA, double TTAG, 
	               double XSV[3], double SANTXYZ[3], double VSVC[3], int IBLK)
{
	int i,j,ii;
	int IECLIPS;
	int IDIR=1;
	double    TWOHR, HALFHR;
	double    ANOON;
	double    CNOON, CNIGHT;
	double    DTR, DTTAG;
	double    MURATE, YANGLE, DET, BETADG, PHI=0.0, SANTX, SANTY, v[3], r[3];
	double    YAWEND;
	double	  ECLSTM, ECLETM;
	int NOON, NIGHT;

	//  MAX YAW RATES OF CURRENT&PAST BLOCK II/IIA's,(AVER'D 1996-2008 JPL SOLUTIONS)  
	//  CHANGE IF REQUIRED OR INPUT IF ESTIMATED 
	//  PRN              01     02     03      04      05      06     07
	double YRATE[]= { .1211, .1339,  .123,  .1233,  .1180,  .1266, .1269,
		// 08     09     10     11      12      13      14     15
		.1033, .1278, .0978, 0.200,  0.199,  0.200, 0.0815, .1303,
		// PRN 16     17     18     19      20      21      22     23
		.0838, .1401, .1069,  .098,   .103, 0.1366,  .1025, .1140,
		// PRN 24     25     26     27      28      29      30     31
		.1089, .1001, .1227, .1194,  .1260,  .1228,  .1165, .0969,
		// PRN 32      33-64: GLONASS RATES (DILSSNER 2010)                            
		.1152,
		0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
		0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
		0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
		0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250
	};

	ECLSTM=ECLETM=-1e6;

	//  CHECK FOR BLOCK IIR AND FIX TO NOMINAL YAW RATE
	if( IPRN<=MAXPRNGPS && IBLK>=4 ) YRATE[IPRN-1]=0.2;

	// THE NEW GPS BLK IIF YAW RATES ( DILSSNER (2010) INSIDE GNSS)
	if( IPRN<=MAXPRNGPS && IBLK>=6 ) YRATE[IPRN-1]=0.11;

	IECLIPS=0;

	TWOHR = 7200.0;
	HALFHR= 1800.0;    
	DTR=D2R;

	// compute the noon beta angle limit (beta zero) FOR A NOON TURN from YRATEs
	// & THE ACTUAL SAT ORBIT ANGLE RATE (MURATE) (~0.00836 FOR GPS; ~ 0.00888 GLNS)
	MURATE= sqrt( ( pow(VSVC[1-1],2) + pow(VSVC[2-1],2) + pow(VSVC[3-1],2) ) / 
		( pow(XSV[1-1], 2) + pow(XSV[2-1], 2) + pow(XSV[3-1], 2) ) 
		) / DTR;
	ANOON=atan(MURATE/YRATE[IPRN-1])/DTR;

	CNOON=cos(ANOON*DTR);
	CNIGHT=cos(ANIGHT*DTR);

	//
	NOON=0;
	NIGHT=0;
	BETADG = BETA/DTR - 90.0; 

	if ( 40==IPRN && 14==PPP_Glo.ctNow[3] ) {
		ii=0;
	}

	//
	if ( IPRN>MAXPRNGPS && fabs(BETADG)<ANOON ) {
		// GLONASS NOON TURN MODE ACORDING TO DILSSNER 2010 
		YAWEND=75.0;
		//  ITERATION FOR YAWEND OF THE GLONASS  NOON TURN

		for (j=1;j<=3;j++) {
			YAWEND=fabs(   atan2( -tan(BETADG*DTR), sin(PI-DTR*MURATE*YAWEND/YRATE[IPRN-1]) ) / DTR 
				- atan2( -tan(BETADG*DTR), sin(PI+DTR*MURATE*YAWEND/YRATE[IPRN-1]) ) / DTR 
				)/2.0; 
		}

		// UPDATE ANOON, CNOON FOR NEW GLONASS NOON TURN LIMITS
		ANOON= MURATE*YAWEND/YRATE[IPRN-1];
		CNOON= cos(ANOON*DTR);
	}

	// BLK IIR'S
	if( IBLK==4 || IBLK==5 ) {
		CNIGHT=cos((ANOON+180.0)*DTR);
		for (j=1;j<=3;j++) {
			// BODY-X U VECTOR REVERSAL FOR IIR ONLY
			SANTXYZ[j-1]=-SANTXYZ[j-1];
		}
	}
	//
	if ( SVBCOS < CNIGHT ) 
		NIGHT=1;

	if ( SVBCOS > CNOON )
		NOON=1;

	//
	//     IF SV IN NIGHT SHADOW OR NOON TURN DURING FORWARD PASS
	//     STORE START AND END TIME OF YAW MANEUVRE (FOR THE BACKWARD RUN)
	//

	//acos: 0-pi
	// YAW ANLGE
	YANGLE= acos( (SANTXYZ[1-1]*VSVC[1-1] + SANTXYZ[2-1]*VSVC[2-1] + SANTXYZ[3-1]*VSVC[3-1] ) / 
		sqrt( pow(VSVC[1-1],2) + pow(VSVC[2-1],2) + pow(VSVC[3-1],2) )
		) / DTR;

	// IIR YANGLE has the same sign as beta, II/IIA has the opposite sign
	if( BETADG<0.0 && IBLK>=4 && IBLK<=5 )
		YANGLE=-YANGLE;
	if( BETADG>0.0 && IBLK!=4 && IBLK!=5 )
		YANGLE=-YANGLE;

	//
	if( (NIGHT || NOON) ) {
		DET=SQRT( pow(180.0-acos(SVBCOS)/DTR, 2) - pow(BETADG,2) );
		PHI = PI/2.0;
		// Check if already after a midnight or noon
		if ( NIGHT ) {
			if (IBLK==4 || IBLK==5) {
				if ( fabs(YANGLE)>90.0 )	DET=-DET;
				if ( DET!=0.0 )				PHI=atan2( tan(BETADG*DTR),-sin(-DET*DTR) )/DTR;
			}
			else {
				// BLK IIA & GLONASS TOO !
				if ( fabs(YANGLE)<90.0 )	DET=-DET;
				if ( DET!=0.0 )				PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR) )/DTR;
			}
		}
		if( NOON ) {
			DET=SQRT( pow(acos(SVBCOS)*180.0/PI,2) - pow(BETADG,2) );

			if( IBLK==4 || IBLK==5 ) {
				if ( fabs(YANGLE)<90.0 )	DET=-DET;
				if ( DET!=0.0 )				PHI=atan2(tan(BETADG*DTR), -sin(PI-DET*DTR))/DTR;
			}
			else {
				// BLK IIA & GLONASS !
				if( fabs(YANGLE)>90.0 )		DET=-DET;
				if( DET!=0.0 )				PHI=atan2(-tan(BETADG*DTR),sin(PI-DET*DTR))/DTR;
			}
		}



		// ONLY FORWARD
		//计算正午/午夜机动，地影恢复期时间段
		if (IDIR > 0 ) {
			//
			// INITIALIZE ECLIPSE START AND TIME TAG ARRAYS  
			//
			//1
			//if ( NECLIPS[IPRN-1] == 0 ) 
			{
				//NECLIPS[IPRN-1]=NECLIPS[IPRN-1]+1;
				ECLSTM=TTAG+DET/MURATE;
				// IIR MIDNIGHT/NOON TURN or II/IIA NOON TURN START
				// for IIR/GLONAS NIGHT (turn) only makes sense when BETADG < ANOON!
				// For IIA it gets here only when NOON is true and that happens  only when BETADG < ANOON!
				YAWEND=atan(MURATE/YRATE[IPRN-1])/DTR;

				if(((IBLK>3 && IBLK<=5) || NOON) && fabs(BETADG)<YAWEND ) {
					// GLONASS
					if ( IPRN > MAXPRNGPS ) {
						// GLONASS NOON TURN MODE ACORDING TO DILSSNER ET AL 2010 
						ECLSTM = ECLSTM - ANOON/MURATE;
						ECLETM = ECLSTM + 2.0*ANOON/MURATE;
					}
					else {
						// GPS IIA/IIR/IIF NOON OR IIR MIDNIGHT TURNs
						ECLSTM = ECLSTM -   fabs(BETADG)*sqrt(ANOON/fabs(BETADG)-1.0)/MURATE;
						ECLETM = ECLSTM + 2*fabs(BETADG)*sqrt(ANOON/fabs(BETADG)-1.0)/MURATE;
					}
				}

				// II/IIA SHADOW START & END TIMES
				if ( (IBLK<=3 || IBLK>5) && NIGHT ) {
					//if (ANIGHT<180)
					//	ANIGHT+=180;


					ECLSTM = ECLSTM -     SQRT( pow(ANIGHT-180.0,2) - pow(BETADG,2) )/MURATE;
					ECLETM = ECLSTM + 2.0*SQRT( pow(ANIGHT-180.0,2) - pow(BETADG,2) )/MURATE;
				}
				//
				// UPDATE SV COSINE AND TIME TAG ARRAYS
				// (TO BE USED DURING BACKWARDS RUN)
				//
				if ( (NIGHT && SVBCOS<CNIGHT) || (NOON && SVBCOS>CNOON) ) {
					DTTAG= fabs(TTAG-ECLSTM);
					//
					// ECLIPSE TIME IS MORE THAN 2 HOURS, THIS IS A NEW ECLIPSE!
					//
					if ( DTTAG>TWOHR ) {
						ECLSTM=TTAG+DET/MURATE;
						// IIR MIDNIGHT/NOON TURN  or II/IIA NOON TURN START
						// AND GLONASS NOON
						if ((IBLK>3 && IBLK<=5) || NOON) {
							// GLONASS

							if (IPRN>MAXPRNGPS) {
								// GLONASS NOON TURN MODE ACORDING TO DILSSNER ET AL 2010 
								ECLSTM = ECLSTM -     ANOON/MURATE;
								ECLETM = ECLSTM + 2.0*ANOON/MURATE;
							}
							else {
								// GPS TURNS ONLY
								ECLSTM = ECLSTM -  fabs(BETADG)*sqrt(ANOON/fabs(BETADG)-1.0)/MURATE;
								ECLSTM = ECLSTM +2*fabs(BETADG)*sqrt(ANOON/fabs(BETADG)-1.0)/MURATE;
							}
						}
					}

					//     II/IIA SHADOW START & END TIMES
					//   & GLONASS & IIF AS WELL !
					if ( (IBLK<=3||IBLK>5) && NIGHT) {
						//if (ANIGHT<180)
						//	ANIGHT+=180;

						ECLSTM = ECLSTM -     SQRT( pow(ANIGHT-180.0,2)-pow(BETADG,2) )/MURATE;
						ECLSTM = ECLSTM + 2.0*SQRT( pow(ANIGHT-180.0,2)-pow(BETADG,2) )/MURATE;
					}
				}
			}
			//  END OF FORWARD LOOP (IDIR = 1)
		}
	}
	//
	//     BOTH FWD (IDIR= 1) OR BWD (IDIR=-1)
	//     SET ECLIPSE FLAG (1=NIGHT SHADOW, 2=NOON TURN) 
	//
	if ( 1 ) {
		// CHECK IF IPRN IS ECLIPSING AND WHICH SEQ NO (I)
		i=0;


		for ( j=1;j<=1;j++ ) {
			if ( fabs(ECLETM+1.0e6)<=1.0e-8 && fabs(ECLSTM+1.0e6)<=1.0e-8 )
				continue;			

			if ( TTAG>=ECLSTM && TTAG<=(ECLETM+HALFHR) )
				i=j;
		}

		// CURRENTLY NOT ECLIPSING (i=0)
		if ( 0==i ) return IECLIPS;



		//判断此时时间是否在正午/午夜机动，地影恢复期
		if ( TTAG>=ECLSTM && TTAG<=(ECLETM+HALFHR) ) {
			// velocity & radius unit vectors V & R
			for ( j=1;j<=3;j++ ) {
				v[j-1]=VSVC[j-1]/SQRT( pow(VSVC[1-1],2)+pow(VSVC[2-1],2)+pow(VSVC[3-1],2) );
				r[j-1]=XSV[j-1] /SQRT( pow(XSV[1-1], 2)+pow(XSV[2-1], 2)+pow(XSV[3-1], 2) );
			}
			// ORBIT ANGLE MU AT ECLIPSE/TURN START
			DET= MURATE*(ECLETM-ECLSTM)/2.0;


			//！！！！！！！！！！！！计算此时具体的航偏角PHI，将名义姿态绕航偏角旋转
			if (SVBCOS < 0) {
				// SHADOW CROSSING
				// BLK IIA/IIF SHADOW CROSSING
				if ( IPRN<=MAXPRNGPS && (IBLK<=3||IBLK>5) ) {
					if ( TTAG<=ECLETM ) {
						// IIA NIGHT TURN
						if ( IBLK<=3 )  PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR + sign(YRATE[IPRN-1],0.50)*(TTAG-ECLSTM); 
						// IIF NIGHT TURN (DILSSNER  2010)
						if ( IBLK>5 )   PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR + sign(0.060, BETADG)*(TTAG-ECLSTM); 
					}
					else {
						// **** WARNING
						// IIA/IIF SHADOW EXIT RECOVERY: USING THE IIA DATA  DURING
						// THE IIA RECOVERY (UP TO 30 MIN) IS NOT RECOMMENDED!
						// **** WARNING
						// GPS IIA  AT SHADOW EXIT
						if ( IBLK<=3 )  PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR + sign(YRATE[IPRN-1],0.50)*(ECLETM-ECLSTM); 
						// GPS IIF AT SHADOW EXIT
						if ( IBLK>5 )   PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR + sign(0.060, BETADG)*(ECLETM-ECLSTM); 
						// YAWEND- HERE THE ACTUAL  YAW  AT THE SHADOW EXIT
						YAWEND= YANGLE- PHI; 
						YAWEND= fmod(YAWEND, 360.0);
						if ( fabs(YAWEND)>180.0) YAWEND= YAWEND-360.0*YAWEND/fabs(YAWEND);
						PHI=PHI + sign(YRATE[IPRN-1],YAWEND)*(TTAG-ECLETM); 
						// SANTX- THE CURRENT ANGLE DIFF, CONSISTENT WITH YAWEND
						SANTX= YANGLE-PHI;
						SANTX = fmod(SANTX , 360.0);
						if ( fabs(SANTX)>180.0) SANTX = SANTX -360.0* SANTX /fabs(SANTX );
						// STOP! THE NOMINAL YAW (YANGLE) REACHED!
						if ( fabs(SANTX)>fabs(YAWEND) ) return IECLIPS;
						if ( YAWEND!=0.0 && ((SANTX)/YAWEND)<0.0) return IECLIPS;
						// SET PHI <-180,+180>
						PHI= fmod(PHI, 360.0);
						if ( fabs(PHI)>180.0) PHI= PHI-360.0*PHI/fabs(PHI);
					}
				}

				// GLONASS
				if( IPRN>MAXPRNGPS ) {
					// GLONASS/GPS  NIGHT TURN (DILSSNER AT AL 2010 )
					if ( TTAG>ECLETM ) return IECLIPS;
					YAWEND=YRATE[IPRN-1]; 
					PHI=atan2(-tan(BETADG*DTR), sin(-DET*DTR))/DTR + sign(YAWEND,BETADG)*(TTAG-ECLSTM); 
					// YAWEND -YAW ANGLE AT THE (GLONASS) SHADOW EXIT
					YAWEND=atan2(-tan(BETADG*DTR), sin( DET*DTR))/DTR;

					if ((YAWEND/PHI)>=1.0 || (PHI/YAWEND)<0.0)   
						PHI = YAWEND;
				}

				if ( IPRN<=MAXPRNGPS && IBLK>5 ) 
					// GPS BLK IIF NIGHT YAW RATE(DILSSNER 2010):
					if ( fabs(BETADG)>8.0) return IECLIPS;

				if ( IBLK>3 && IBLK<=5) {
					// BLK II R SHADOW (MIDNIGHT TURN) CROSSING
					PHI=atan2( tan(BETADG*DTR),-sin(-DET*DTR))/DTR + sign(YRATE[IPRN-1],BETADG)*(TTAG-ECLSTM); 
					if( (PHI/YANGLE)>=1.0 || (PHI/YANGLE)<0.0) return IECLIPS;
				}
				//             write(*,*)"R",IPRN-32,TTAG,YANGLE, PHI,DET,
				//    & BETADG, ECLETM(IPRN,I),I
				IECLIPS=1;
			}
			else {
				// NOON TURNS 
				PHI=atan2(-tan(BETADG*DTR),sin(PI-DET*DTR))/DTR -sign(YRATE[IPRN-1],BETADG)*(TTAG-ECLSTM);
				if ( IBLK>3 && IBLK<=5 ) {
					// BLK IIR NOON TURNS ONLY
					PHI=atan2( tan(BETADG*DTR),-sin(PI-DET*DTR))/DTR -sign(YRATE[IPRN-1],BETADG)*(TTAG-ECLSTM);
					// IIR END TURN CHECK
					if ( (YANGLE/PHI)>=1.0 || (PHI/YANGLE)<0.0) return IECLIPS;
				}
				else {
					// GLONASS END TURN CHECK
					if ( IPRN>MAXPRNGPS && TTAG>ECLETM ) return IECLIPS;
					// IIA OR IIF END TURN CHECK
					if ( IPRN<=MAXPRNGPS && ((PHI/YANGLE)>=1.0 || (PHI/YANGLE)<0.0)) return IECLIPS;
				}
				//             write(*,*)"S",IPRN-32,TTAG,YANGLE, PHI,DET,
				//    & BETADG, ECLSTM(IPRN,I)
				IECLIPS=2;
			}
			// ROTATE X-VECTOR TO ECLIPSING YAW ANGLE PHI 
			// ECLIPSING (II/IIA) NOT TO BE USED  A HALF HR AFTER SHADOW !
			SANTX=(cos((PHI-YANGLE)*DTR)*(v[2-1]-v[3-1]*r[2-1]/r[3-1])-cos(PHI*DTR)*(SANTXYZ[2-1]-SANTXYZ[3-1]*r[2-1])/r[3-1])/(SANTXYZ[1-1]*v[2-1]-SANTXYZ[2-1]*v[1-1])
				+((SANTXYZ[2-1]*v[3-1]-SANTXYZ[3-1]*v[2-1])*r[1-1]+(SANTXYZ[3-1]*v[1-1]-SANTXYZ[1-1]*v[3-1])*r[2-1])/r[3-1];
			SANTY = (cos(PHI*DTR) - (v[1-1]-v[3-1]*r[1-1]/r[3-1])*SANTX)/(v[2-1]-v[3-1]*r[2-1]/r[3-1]);
			// THE BODY-X UNIT VECTOR ROTATED BY (PHI-YANGLE) RETURNED
			SANTXYZ[1-1]= SANTX;
			SANTXYZ[2-1]= SANTY;
			SANTXYZ[3-1]= (-r[1-1]*SANTX-r[2-1]*SANTY)/r[3-1];
		}
	}

	return IECLIPS;
}

extern int calEclips(int prn, double *satp, const double *satv, double *sunp, 
	                 double TTAG, double SANTXYZ[3], const nav_t *nav)
{
	double SVBCOS,BETA=0.0,eSunP[3],eSatP[3],eSatV[3],vec[3],ANIGHT,satv_[3];
	double angleLmt;
	const char *type;
	int IBLK=-1;		//IBLK SV BLOCK  1=I, 2=II, 3=IIA, IIR=(4, 5) IIF=6

	if (prn>MAXPRNGPS+MAXPRNGLO) return 0;

	type=nav->pcvs[prn-1].type;

	if (type) {
		if (strstr(type, "BLOCK I           "))	IBLK=1;
		else if (strstr(type, "BLOCK II     "))	IBLK=2;
		else if (strstr(type, "BLOCK IIA    "))	IBLK=3;
		else if (strstr(type, "BLOCK IIR-B  "))	IBLK=4;
		else if (strstr(type, "BLOCK IIR-M  "))	IBLK=5;
		else if (strstr(type, "BLOCK IIF    "))	IBLK=6;
	}

	if (prn>MAXPRNGPS) IBLK=6;

	satv_[0]=satv[0]-OMGE*satp[1];
	satv_[1]=satv[1]+OMGE*satp[0];
	satv_[2]=satv[2];

	normv3(satp,eSatP);
	normv3(satv_,eSatV);
	normv3(sunp,eSunP);

	SVBCOS=dot(eSatP,eSunP,3);

	cross3(eSatP,eSatV,vec);

	BETA=dot(vec,eSunP,3);
	BETA=acos(BETA);
	BETA=-BETA+PI;

	angleLmt=76.116;  //acos(RE_WGS84/26580000.0)*R2D;

	ANIGHT=90+angleLmt-1.5;

	return eclips_(prn,SVBCOS,ANIGHT,BETA,TTAG,satp,SANTXYZ,satv_,IBLK);
}

extern int lsqPlus(const double *A, const double *y, const int nx, const int nv, double *x, double *Q)
{
	int i,j,k,info=0;
	int *ix;
	double *A_,*x_,*Q_;

	ix=imat(nx,1);

	for (i=k=0;i<nx;i++) {
		for (j=0;j<nv;j++) {
			if (fabs(A[j*nx+i])>1.0e-10) {
				ix[k++]=i;
				j=1000;
			}
		}
	}

	A_=mat(k*nv,1); x_=mat(k,1); Q_=mat(k*k,1);

	for (j=0;j<k;j++) {
		for (i=0;i<nv;i++) {
			A_[i*k+j]=A[i*nx+ix[j]];
		}

		x_[j]=x[ix[j]];
	}

	/* least square estimation */
	info=lsq(A_,y,k,nv,x_,Q_);

	for (i=0;i<nx*nx;i++) Q[i]=0.0;

	for (i=0;i<k;i++) {
		x[ix[i]]=x_[i];

		for (j=0;j<k;j++)
			Q[ix[i]+ix[j]*nx]=Q_[i+j*k];
	}

	free(ix); free(A_); free(x_); free(Q_);

	return info;
}

//to compute sagnac effect correction
extern double sagnac(const double *rs, const double *rr)
{
	return OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}

//convert satellite PCO to distance
//extern double svpco2dist(gtime_t time, const double *rr, const double *rs, int sat, const nav_t *nav)
//{
//	const pcv_t *pcv=nav->pcvs+sat-1;
//	double rk[3],ri[3],rj[3],r[3],rrho[3],pco[3],rsun[3],gmst,erpv[5]={0};
//	double svant[3],svPCcorr=0.0,ac=0.0,nadir=0.0,elev=0.0;
//	int i,j=0,k=1,sys;
//
//	/* sun position in ecef */
//	sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);
//
//	//Unitary vector from satellite to Earth mass center (ECEF)
//	for (i=0;i<3;i++) r[i]=-rs[i];
//	if (!normv3(r,rk)) return 0.0;
//
//	//Unitary vector from Earth mass center to Sun (ECEF)
//	if (!normv3(rsun,ri)) return 0.0;
//
//	//rj = rk x ri: Rotation axis of solar panels (ECEF)
//	cross3(rk,ri,rj);
//
//	//Redefine ri: ri = rj x rk (ECEF)
//	cross3(rj,rk,ri);
//	//Let's convert ri to an unitary vector. (ECEF)
//	if (!normv3(ri,ri)) return 0.0;
//
//	//Compute unitary vector vector from satellite to RECEIVER
//	for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
//	if (!normv3(r,rrho)) return 0.0;
//
//	ac=dot(rrho,rk,3);
//	ac=ac<-1.0?-1.0:(ac>1.0?1.0:ac);
//	nadir=acos(ac)*R2D;
//	elev=90.0-nadir;
//
//	sys=satsys(sat,NULL);
//	if (sys==SYS_GPS) {
//		j=0;
//		k=1;
//	}
//	else if (sys==SYS_GLO) {
//		j=0+NFREQ;
//		k=1+NFREQ;
//	}
//	else if (sys==SYS_CMP) {
//		j=0+2*NFREQ;
//		k=1+2*NFREQ;
//	}
//	else if (sys==SYS_GAL) {
//		j=0+3*NFREQ;
//		k=1+3*NFREQ;
//	}
//	for (i=0;i<3;i++) pco[i]=pcv->off[j][0];
//	for (i=0;i<3;i++) svant[i]=ri[i]*pco[0]+rj[i]*pco[1]+rk[i]*pco[2];
//	svPCcorr=dot(rrho,svant,3);
//
//	return svPCcorr;
//}
