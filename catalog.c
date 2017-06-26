/**
 * @file    catalog.c
 * @brief   High contrast imaging target list
 * 
 * To compile:
 * gcc catalog.c -lm
 * 
 * READS FOLLOWING FILES :
 *   SUPERBLINK.dat           : SUPERBLINK catalog
 *   CNS3-Gliese_catalog.dat  : Gliese catalog
 *   Gliese2MASS_clean.dat    : 2MASS Gliese
 *   RECONS100.dat            : RECONS catalog
 *   hip_main.dat             : HIPPARCOS catalog (main)
 *   hip_main2.dat            : HIPPARCOS catalog 
 *   hip_doubles.dat          : HIPPARCOS double stars
 * 
 *   miscdata.dat             : custom entries
 * 
 *   corleak.dat              : CORONAGRAPH LEAK PROFILE FOR 0.01 l/D RADIUS STAR, VALUE EVERY TENTH OF L/D
 *   conf_spectralband.txt    : spectral band (single character, for eg, "H"
 *   coroData_0.dat           : coherent light contrast
 *   coroData_1.dat	          : incoherent light contrast
 * 
 * @author  O. Guyon
 * @date    26 Jun 2017
 *
 * @bug No known bugs. 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define SBUFFERSIZE 1000



int simulMode = 1;                        /**< write colors and magn from EEM data */

// TELESCOPE AND CORONAGRAPH PARAMETERS

int VERBOSE = 0;

double  telDiam = 30.0;                    /**< Telescope diameter, [m] */
double  Efficiency = 0.5;                  /**< Efficiency (excludes coronagraph concept throughput */
double spectralR = 100.0;                  /**< Spectral resolution */
   
char spectralBand = 'H';                   /**< spectral band */
double requ_SNR = 10.0;                    /**< Required SNR at spectral resolution */
double requ_etime = 1.0*24.0*3600.0;       /**< Exposure time within which SNR must be met */
double ContrastQlim = 1000.0;              /**< maximum calibration factor on background, used to compute required RAW contrast needed to meet SNR */


// CORONAGRAPH 
double IWAld = 1.2;                        /**< Inner Working Angle */
double FWHMld = 1.11;                      /**< Full Width at Half Max post-coronagraph. Planet image is wider due to Lyor stop and/or apodization */


// Throughput is for incohere background (zodi and exozodi)
// AiryThroughput is for planet light
double Throughput = 0.8;                   /**< Coronagraph concept throughput */
double AiryThroughput = 0.7;               /**< Airy Throughput */


// Instrumental RAW contrast (excluding angular size) levels are specified inside and outside IWAld1:
// if r < IWAld        -> no detection
// IWAld < r < IWAld1  -> raw contrast = RAWcontrast1
// r > IWAld1          -> raw contrast = RAWcontrast2 
// 
// Angular size leak is read from file corleak.dat and added to RAWcontrast1 or RAWcontrast2 and written into RAWcontrast
double RAWcontrast;                        /**< Variable holding total incoherent contrast */              
double RAWcontrast1 = 1.0e-6;
double RAWcontrast2 = 1.0e-7;
double IWAld1 = 2.0;                      /**< IWA for higher contrast */


double contrastFactor = 100.0;            /**< Maximum calibration gain between raw contrast and detection (excludes photon noise) */


float lambdaB = 0.436e-6;
float lambdaV = 0.545e-6;
float lambdaR = 0.638e-6;
float lambdaI = 0.797e-6;
float lambdaJ = 1.22e-6;
float lambdaH = 1.63e-6;
float lambdaK = 2.19e-6;

double zodibckgB = 22.0 + 0.642; 
double zodibckgV = 22.0; // mag per squ arcsec in V [ref] 
double zodibckgR = 22.0; // ?? 
double zodibckgI = 22.0 - 0.706; 
double zodibckgJ = 22.0 - 1.18; 
double zodibckgH = 22.0 - 1.473; 
double zodibckgK = 22.0 - 1.545; 

double BVsun = 0.642;
double VIsun = 0.706;
double VJsun = 1.18;
double VHsun = 1.473;
double VKsun = 1.545;

#define DISTMAX 30.0 // max distance in pc

long NBSOURCES;

// A Modern Mean Stellar Color and Effective Temperatures (Teff) - EEM
// Only valid for dwarfs (type V)
float stfarray_x[67]; // Spectral type
float stfarray_Teff[67]; // Teff
float stfarray_logT[67];
float stfarray_BCv[67];
float stfarray_Mv[67];
float stfarray_logL[67];
float stfarray_UB[67];
float stfarray_BV[67];
float stfarray_VI[67];
float stfarray_VK[67];
float stfarray_JH[67];
float stfarray_HK[67];
float stfarray_KW[67];
float stfarray_Msun[67];
float stfarray_lgAge[67];
float stfarray_by[67];
float stfarray_MJ[67];
float stfarray_MK[67];
float stfarray_Mbol[67];
float stfarray_KL[67]; 



float simulBCv;  
float simulx;
float simulUB;
float simulBV;
float simulVI;
float simulVJ;
float simulVH;
float simulVK;
float simulJH;
float simulHK;


float simulMsun;
 
float valx, errx;
float valUB, errUB;
float valBV, errBV;
float valVI, errVI;
float valVJ, errVJ;
float valVK, errVK;
float valJH, errJH;
float valHK, errHK;
float valTeff, errTeff;

double RA1950, DEC1950, RA2000, DEC2000; // [rad]





long ContrastDATA_NBpt = 100;
float ContrastDATA_step = 0.2;
float *ContrastCOH;
float *ContrastINCOH;
float SSIZE_INCOH = 1.0;











typedef struct
{
  int FLAG_SUPERBLINK;
  int FLAG_2MASSCNS3;
  int FLAG_CNS3;
  int FLAG_HIPS;
  int FLAG_HIPD;
  int FLAG_RECONS;

  char name[200];
  char name1[200];
  char name2[200];
  char name3[200];

  float RAdeg;
  float DECdeg;

  float dist; // pc
  float disterr; // pc

  float VJ; // V-J
  float BTmag; // Tycho B mag
  float VTmag; // Tycho V mag
  float BJmag; // USNO B mag
  float RFmag; // USNO R mag
  float INmag; // USNO I mag

  float Vemag; // SUPERBLINK effective Vmag

  float Jmag; // 2MASS
  float Hmag; // 2MASS
  float Kmag; // 2MASS
  float SpectralIndex; // from SUPERBLINK .. -1 if unknown

  float VI;
  float VK;
  float VH;
  float JH;
  float HK;

  float TeffM; // Teff measurement if exists

  char SUPERBLINKCNS3name1[200];
  char SUPERBLINKCNS3name2[200];
  char SUPERBLINKCNS3name3[200];
  int SUPERBLINKHIPindex;

  char CNS3name1[200];
  char CNS3name2[200];
  char CNS3name3[200];
  float SpType; // from CNS3 ... -1 if unknown
  float mV; // from CNS3
  float BV; // from CNS3
  float UB; // from CNS3
  float RI; // from CNS3
  float CNS3dist; // from CNS3
  float MV; // from CNS3 if available

  long HIPindex;
  float HIPHpmag; 
  float HIPVmag; // HIPPARCOS Vmag
  float HIPSpType;
  float HIPBT;
  float HIPVT;
  float HIPBV;
  float HIPVI;

  char RECONSname1[200];
  char RECONSname2[200];
  char RECONSname3[200];
  float RECONSmV;
  float RECONSMV;
  float RECONSdist;
  float RECONSSpType;
  float RECONSmass;

  float Teff; // derived from colors and Teff estimate
  float Lmag; // L band mag
  float simulx; // simulated spectral type
  float Blum; // Bolometric Luminosity
  float HZSeparcsec;
  float HZContrast;
  float Msun; // star mass 
} STAR;

STAR *star;
long NBstar;



typedef struct
{
  char name[200];
  int nameA;
  int nameB;
  int nameC;
  char name1[200];
  char name2[200];
  char name3[200];

  float RAdeg;
  float DECdeg;
 
  float SpType; // -1 in not standard OBAFGKM, 0 for WD
  float mV;
  float BV;
  float UB;
  float RI;
  float dist;
  float MV;

} CNS3STAR;

CNS3STAR *CNS3star;
long NBstarCNS3;


typedef struct
{
  char name[200];
  char name1[200];
  char name2[200];
  char name3[200];

  float RAdeg;
  float DECdeg;

  float Jmag;
  float Hmag;
  float Kmag;

} CNS3_2MASS_STAR;

CNS3_2MASS_STAR *CNS3_2MASS_star;
long NBstar_CNS3_2MASS;


typedef struct
{
  long index;

  float Vmag;
  float Hpmag;

  float RAdeg;
  float DECdeg;

  float dist;
  float BT;
  float VT;
  float BV;
  float VI;

  float SpType;

} HIPSTAR;

HIPSTAR *HIPstar;
long NBstarHIP;



typedef struct
{
  int index;
  int component;

  float Hpmag;

  float RAdeg;
  float DECdeg;

  float dist;
  float BT;
  float VT;
  float rho;
} HIPDSTAR;

HIPDSTAR *HIPDstar;
long NBstarHIPD;


typedef struct
{
  char name[200];
  float RAdeg;
  float DECdeg;
  float dist;  
} WDSTAR;

WDSTAR *WDstar;
long NBstarWD;


typedef struct
{
  char name1[200];
  char name2[200];
  char name3[200];
  float RAdeg;
  float DECdeg;
  float dist;
  float SpType;
  float mV;
  float MV;
  float mass;
} RECONSSTAR;

RECONSSTAR *RECONSstar;
long NBstarRECONS;







int read_ContrastData()
{
	long NBpt;
	FILE *fp0;
	FILE *fp1;
	float ssize0 = 0.0001; // l/D radius
	float ssize1 = 0.0316227766;
	long i;
	float v00, v01, v02;
	float v10, v11, v12;
	char fname[200];


	SSIZE_INCOH = ssize1;
	NBpt = ContrastDATA_NBpt;

	ContrastCOH = (float*) malloc(sizeof(float)*NBpt);
	ContrastINCOH = (float*) malloc(sizeof(float)*NBpt);

	sprintf(fname, "coroData_0.dat");
	fp0 = fopen(fname, "r");
	if(fp0==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
	

	sprintf(fname, "coroData_1.dat");	
	fp1 = fopen(fname, "r");
	if(fp1==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
		
	for(i=0;i<NBpt;i++)
		{
			fscanf(fp0, "%f %f %f", &v00, &v01, &v02);
			fscanf(fp1, "%f %f %f", &v10, &v11, &v12);
			ContrastCOH[i] = v02;
			ContrastINCOH[i] = v12;
		}
	fclose(fp0);
	fclose(fp1);

	return(0);
}



int read_EEMtable()
{
  FILE *fp;
  char stf_string[20];
  int i;
  char word[200];
	char fname[200];


  // NOTE: ALL UNDEFINED VALUES HAVE BEEN REPLACED BY -100
  // NOTE: K-L COLORS ADDED BY HAND IN FILE FROM ALLEN 
	sprintf(fname, "EEM_dwarf_UBVIJHK_colors_Teff.dat1");

  fp = fopen(fname, "r");
  if(fp==NULL)
  {
	  printf("ERROR: File \"%s\" missing\n", fname);
	exit(0);
  }
  
  for(i=0;i<67;i++)
    {
      fscanf(fp,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f\n", stf_string, &stfarray_Teff[i], &stfarray_logT[i], &stfarray_BCv[i], &stfarray_Mv[i], &stfarray_logL[i], &stfarray_UB[i], &stfarray_BV[i], &stfarray_VI[i], &stfarray_VK[i], &stfarray_JH[i], &stfarray_HK[i], &stfarray_KW[i], &stfarray_Msun[i], &stfarray_lgAge[i], &stfarray_by[i], stf_string, &stfarray_MJ[i], &stfarray_MK[i], &stfarray_Mbol[i], &stfarray_KL[i]);
      word[0] = stf_string[1];
      word[1] = '\0';
      switch ( stf_string[0] ) {
      case 'O':
	stfarray_x[i] = 1.0 + 0.1*atoi(word);
	break;
      case 'B':
	stfarray_x[i] = 2.0 + 0.1*atoi(word);
	break;
      case 'A':
	stfarray_x[i] = 3.0 + 0.1*atoi(word);
	break;
      case 'F':
	stfarray_x[i] = 4.0 + 0.1*atoi(word);
	break;
      case 'G':
	stfarray_x[i] = 5.0 + 0.1*atoi(word);
	break;
      case 'K':
	stfarray_x[i] = 6.0 + 0.1*atoi(word);
	break;
      case 'M':
	stfarray_x[i] = 7.0 + 0.1*atoi(word);
	break;
      default:
	printf("ERROR in read_EEMtable: %d  \'%c\' \"%s\"\n", i, stf_string[0], stf_string);
	exit(0);
	break;
      }
    }

  fclose(fp);

  /*  for(i=0;i<67;i++)
    printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", stfarray_x[i], stfarray_Teff[i], stfarray_logT[i], stfarray_BCv[i], stfarray_Mv[i], stfarray_logL[i], stfarray_UB[i], stfarray_BV[i], stfarray_VI[i], stfarray_VK[i], stfarray_JH[i], stfarray_HK[i], stfarray_KW[i], stfarray_Msun[i], stfarray_lgAge[i], stfarray_by[i], stfarray_MJ[i], stfarray_MK[i], stfarray_Mbol[i], stfarray_KL[i]);

    exit(0);*/

  return(0);
}


int SimulMS_colors(float Teff, int dispR)
{
  long i0, i1;
  float ifrac;

  i1 = 1;
  while(stfarray_Teff[i1]>Teff)
    {
      i1++;
      if(i1==67)
	{
	  printf("ERROR in SimulMS_colors\n");
	  exit(0);
	}
    }
  i0 = i1-1;
  ifrac = (Teff-stfarray_Teff[i0])/(stfarray_Teff[i1]-stfarray_Teff[i0]);

  simulBCv = (1.0-ifrac)*stfarray_BCv[i0] + ifrac*stfarray_BCv[i1];
  
  simulx = (1.0-ifrac)*stfarray_x[i0] + ifrac*stfarray_x[i1];
  simulUB = (1.0-ifrac)*stfarray_UB[i0] + ifrac*stfarray_UB[i1];
  simulBV = (1.0-ifrac)*stfarray_BV[i0] + ifrac*stfarray_BV[i1];
  simulVI = (1.0-ifrac)*stfarray_VI[i0] + ifrac*stfarray_VI[i1];
  simulVK = (1.0-ifrac)*stfarray_VK[i0] + ifrac*stfarray_VK[i1];
  simulJH = (1.0-ifrac)*stfarray_JH[i0] + ifrac*stfarray_JH[i1];
  simulHK = (1.0-ifrac)*stfarray_HK[i0] + ifrac*stfarray_HK[i1];
  
  simulVJ = simulVK - simulJH - simulHK;
  simulVH = simulVK - simulHK;
  
  simulMsun = (1.0-ifrac)*stfarray_Msun[i0] + ifrac*stfarray_Msun[i1];


  if(dispR==1) // print result
    {
      printf("x  : %f %ld %f %ld %f\n", ifrac, i0, stfarray_x[i0], i1, stfarray_x[i1] );
      printf("BV : %f %ld %f %ld %f\n", ifrac, i0, stfarray_BV[i0], i1, stfarray_BV[i1] );

      printf("x     -> %f\n", simulx);
      printf("UB    -> %f\n", simulUB);
      printf("BV    -> %f\n", simulBV);
      printf("VI    -> %f\n", simulVI);
      printf("VK    -> %f\n", simulVK);
      printf("JH    -> %f\n", simulJH);
      printf("HK    -> %f\n", simulHK);
    
      printf("Msun  -> %f\n", simulMsun);
      exit(0);
    }
  return(0);
}








float Teff_solve()
{
  float Teff, Teffopt;
  float val, valbest;
  float tmp;
  FILE *fp;


  fp = fopen("TeffFit.txt", "w");
  Teffopt = 2700.0;
  valbest = 1.0e12;

  for(Teff=2700.0; Teff<30000.0;Teff*=1.01)
    {
      val = 0.0;
      SimulMS_colors(Teff, 0);
      
      tmp = (valx-simulx)/errx;
      val += tmp*tmp;
      
      tmp = (valUB-simulUB)/errUB;
      val += tmp*tmp;
      
      tmp = (valBV-simulBV)/errBV;
      val += tmp*tmp;
       
      tmp = (valVI-simulVI)/errVI;
      val += tmp*tmp;
      
      tmp = (valVK-simulVK)/errVK;
      val += tmp*tmp;
    
      tmp = (valJH-simulJH)/errJH;
      val += tmp*tmp;
      
      tmp = (valHK-simulHK)/errHK;
      val += tmp*tmp;
      
      tmp = (valTeff-Teff)/errTeff;
      val += tmp*tmp;

      if (val<valbest)
		{
			valbest = val;
			Teffopt = Teff;
		}
      fprintf(fp,"%f %f\n", Teff, val);
    }
  fclose(fp);

  
  SimulMS_colors(Teffopt, 0);
  if(VERBOSE==1) // print result
    {
      printf("OPTIMAL TEFF = %f\n", Teffopt);
      printf("x     %f %f    -> %f\n", valx, errx, simulx);
      printf("UB    %f %f    -> %f\n", valUB, errUB, simulUB);
      printf("BV    %f %f    -> %f\n", valBV, errBV, simulBV);
      printf("VI    %f %f    -> %f\n", valVI, errVI, simulVI);
      printf("VK    %f %f    -> %f\n", valVK, errVK, simulVK);
      printf("JH    %f %f    -> %f\n", valJH, errJH, simulJH);
      printf("HK    %f %f    -> %f\n", valHK, errHK, simulHK);
      printf("Teff  %f %f    -> %f\n", valTeff, errTeff, Teffopt);
    }

  return(Teffopt);
}






/*
  x = (stf-stfarray_x[i0])/(stfarray_x[i1]-stfarray_x[i0]);
  value = (1.0-x)*stfarray_v[i0] + x*stfarray_v[i1];

  for(i=0;i<15;i++)
    stfarray_Blum[i] = stfarray_Blum[i] / pow(stfarray_v[i]/10000.0,4.0);

  Blum =  (1.0-x)*stfarray_Blum[i0] + x*stfarray_Blum[i1];
  Blum = Blum * pow(value/10000.0,4.0);

  VK = (1.0-x)*stfarray_VK[i0] + x*stfarray_VK[i1];
  HK = (1.0-x)*stfarray_HK[i0] + x*stfarray_HK[i1];
  KL = (1.0-x)*stfarray_KL[i0] + x*stfarray_KL[i1];
*/




int read_SUPERBLINK()
{
  FILE *fp;
  char line[SBUFFERSIZE];
  long i, offset, k;
  char word[200];
	char fname[200];

  int RAhr, RAmin;
  float RAsec;
  int DECdeg, DECmin;
  float DECsec;

  float sign;

  printf("READING SUPERBLINK ... ");


//  fp = fopen("SUPERBLINK.dat","r");
  
  
	sprintf(fname, "SUPERBLINK.dat");	
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  
  k = 0;
  while((fgets(line,SBUFFERSIZE,fp)!=NULL)&&(k<8889))
    {
      if(line[0]=='|')
	{
	  offset = 57;
	  for(i=0;i<8;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';

	  //printf("par = \"%s\" %f\n", word, atof(word));
	  if(atof(word)>0.000001) // USE PARALLAX IF AVAILABLE
	    {
	      star[k].dist = 1.0/atof(word);
	      offset = 66;
	      for(i=0;i<14;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      //printf("parerr = \"%s\" %f\n", word, atof(word));
	      star[k].disterr = star[k].dist - (1.0/(1.0/star[k].dist + atof(word)));
	    }
	  else // USE PHOTOMETRIC PARALLAX
	    {
	      offset = 81;
	      for(i=0;i<13;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      //printf("parphot = \"%s\" %f\n", word, atof(word));
	      star[k].dist = 1.0/atof(word);
	      
	      offset = 95;
	      for(i=0;i<19;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      //     printf("parphoterr = \"%s\" %f\n", word, atof(word));
	      star[k].disterr = star[k].dist - (1.0/(1.0/star[k].dist + atof(word)));
	    }
	  
	  //	  printf("%ld %f\n", k, star[k].dist);

	  if(star[k].dist<2.0*DISTMAX)
	    {

	      offset = 1;
	      for(i=0;i<17;i++)
		star[k].name1[i] = line[i+offset];
	      star[k].name1[i] = '\0';
	      

	      
	      offset = 19;
	      for(i=0;i<2;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      RAhr = atoi(word);
	      //printf("RAhr = \"%s\" %d\n", word, atoi(word));
	      
	      offset = 22;
	      for(i=0;i<2;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      RAmin = atoi(word);
	      //printf("RAmin = \"%s\" %d\n", word, atoi(word));
	      
	      offset = 25;
	      for(i=0;i<5;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      RAsec = atof(word);
	      // printf("RAsec = \"%s\" %f\n", word, atof(word));
	      
	      star[k].RAdeg = 360.0*(1.0*RAhr + 1.0*RAmin/60.0 + 1.0*RAsec/3600.0)/24.0;
	      
	      
	      
	      
	      
	      if(line[31]=='+')
		sign = 1.0;
	      else
		sign = -1.0;
	      offset = 32;
	      for(i=0;i<2;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      DECdeg = atoi(word);
	      // printf("DECdeg = \"%s\" %d\n", word, atoi(word));
	      
	      offset = 35;
	      for(i=0;i<2;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      DECmin = atoi(word);
	      //  printf("DECmin = \"%s\" %d\n", word, atoi(word));
	      
	      offset = 38;
	      for(i=0;i<4;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      DECsec = atof(word);
	      //	  printf("DECsec = \"%s\" %f\n", word, atof(word));
	      
	      star[k].DECdeg = sign*(1.0*DECdeg + 1.0*DECmin/60.0 + 1.0*DECsec/3600.0);
	      
	      

	      
	      
	      offset = 115;
	      for(i=0;i<5;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"     ")==0)
		star[k].Jmag = -100.0;
	      else
		star[k].Jmag = atof(word);
	      
	      offset = 122;
	      for(i=0;i<7;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"       ")==0)
		star[k].VJ = -100.0;
	      else
		star[k].VJ = atof(word);
	      
	      offset = 284;
	      for(i=0;i<5;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"     ")==0)
		star[k].Hmag = -100.0;
	      else
		star[k].Hmag = atof(word);

	      
	      offset = 291;
	      for(i=0;i<5;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"     ")==0)
		star[k].Kmag = -100.0;
	      else
		star[k].Kmag = atof(word);
	      
	      
	      offset = 249;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		star[k].BTmag = -100.0;
	      else
		star[k].BTmag = atof(word);
	      
	      offset = 256;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		star[k].VTmag = -100.0;
	      else
		star[k].VTmag = atof(word);
	      
	      offset = 263;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		star[k].BJmag = -100.0;
	      else
		star[k].BJmag = atof(word);
	      
	      offset = 270;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		star[k].RFmag = -100.0;
	      else
		star[k].RFmag = atof(word);
	      
	      offset = 277;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		star[k].INmag = -100.0;
	      else
		star[k].INmag = atof(word);
	      
	      offset = 297;
	      for(i=0;i<8;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      star[k].Vemag = atof(word);
	      if(strcmp(word,"        ")==0)
		star[k].Vemag = -100.0;
	      else
		star[k].Vemag = atof(word);

	      
	      offset = 157;
	      for(i=0;i<13;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      sscanf(word,"%s %s %s", star[k].SUPERBLINKCNS3name1, star[k].SUPERBLINKCNS3name2, star[k].SUPERBLINKCNS3name3);

	      
	      word[0] = line[307];
	      word[1] = '\0';
	      if(line[306]=='m')
		star[k].SpectralIndex = 7.0 + 0.1*atoi(word);
	      else if(line[306]=='k')
		star[k].SpectralIndex = 6.0 + 0.1*atoi(word);
	      else 
		{
		  printf("ERROR: unknown spectral type char\n");
		  exit(0);
		}
	      
	      //	  printf("%05ld %17s %2.10g %3.10g %3.8g %3.8g %f %f %f %f %f %f %f %f %f %f\n",k, star[k].name1, star[k].RAdeg, star[k].DECdeg, star[k].dist, star[k].disterr, star[k].BTmag, star[k].BJmag, star[k].VTmag, star[k].Vemag, star[k].RFmag, star[k].INmag, star[k].Jmag, star[k].Hmag, star[k].Kmag, star[k].SpectralIndex);
	      
	      star[k].FLAG_SUPERBLINK = 1;
	      k++;
	    }
	}
    }
  NBstar = k;
  fclose(fp);
  
  printf("%ld sources\n", NBstar);

  return(0);
}





// Using Murray 1989 matrix 28
int B1950t0J2000()
{
  double x0, y0, z0;
  double x1, y1, z1;

  double Rxx =  0.9999256794956877;
  double Rxy = -0.0111814832204662;
  double Rxz = -0.0048590038153592;
  double Ryx =  0.0111814832391717;
  double Ryy =  0.9999374848933135;
  double Ryz = -0.0000271625947142;
  double Rzx =  0.0048590037723143;
  double Rzy = -0.0000271702937440;
  double Rzz =  0.9999881946023742;
  
  x0 = cos(RA1950)*cos(DEC1950);
  y0 = sin(RA1950)*cos(DEC1950);
  z0 = sin(DEC1950);

  x1 = Rxx * x0 + Rxy * y0 + Rxz * z0;
  y1 = Ryx * x0 + Ryy * y0 + Ryz * z0;
  z1 = Rzx * x0 + Rzy * y0 + Rzz * z0;
  
  RA2000 = atan2(y1,x1);
  DEC2000 = atan2(z1,sqrt(x1*x1+y1*y1));
  
  return(0);
}






int read_CNS3()
{
  FILE *fp;
  
  int i, j, offset, iend;
  char word[200];
  char word1[200];
  char name1[20];
  int name1A, name1B, name1C;

  char line[SBUFFERSIZE];
  int RAhr, RAmin, RAsec;
  int DECdeg;
  float DECmin;
  float RAdegf, DECdegf;
  float sign;
  long k;
  long cnt0, cnt1;
  float val, val1, val2; 

	char fname[200];

  cnt0 = 0;
  cnt1 = 0;
//  fp = fopen("CNS3-Gliese_catalog.dat","r");
  
  
	sprintf(fname, "CNS3-Gliese_catalog.dat");	
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  
  printf("READING CNS3 ...  ");

  k = 0;
  while(fgets(line,SBUFFERSIZE,fp)!=NULL)
    {

      offset = 109;
      for(i=0;i<5;i++)
	word[i] = line[i+offset];
      word[5] = '\0';
      val = 1000.0/atof(word);

      if(val>0.1)
	{
	  CNS3star[k].dist = 1000.0/atof(word);
      
	  for(i=0;i<10;i++)
	    word[i] = line[i];
	  word[10] = '\0';
	  //   printf("\"%s\" ",word);
	  
	  offset = 0;
	  j = 0;
	  for(i=0;i<2;i++)
	    {
	      if(line[i+offset]!=' ')
		{
		  CNS3star[k].name1[j] = line[i+offset];
		  j++;
		}
	    }
	  
	  offset = 2;
	  j = 0;
	  for(i=0;i<6;i++)
	    {
	      if(line[i+offset]!=' ')
		{
		  CNS3star[k].name2[j] = line[i+offset];
		  j++;
		}
	    }
	  
	  offset = 8;
	  j = 0;
	  for(i=0;i<2;i++)
	    {
	      if(line[i+offset]!=' ')
		{
		  CNS3star[k].name3[j] = line[i+offset];
		  j++;
		}
	    }
	  
	  //      printf("NAME: %s %s %s\n", CNS3star[k].name1,CNS3star[k].name2,CNS3star[k].name3);
	  
	  
	  name1A = 0;
	  name1B = 0;
	  name1C = 0;
	  for(i=0;i<11;i++)
	    {
	      CNS3star[k].name[i] = word[i];
	      if(word[i]=='A')
		{
		  CNS3star[k].name[i] = ' ';
		  CNS3star[k].nameA = 1;
		}
	      if(name1[i]=='B')
		{
		  CNS3star[k].name[i] = ' ';
		  CNS3star[k].nameB = 1;
		}
	      if(name1[i]=='C')
		{
		  CNS3star[k].name[i] = ' ';
		  CNS3star[k].nameC = 1;
		}
	    }
      	

	  offset = 12;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAhr = atoi(word);
	  //printf("RAhr = \"%s\" %d\n", word, atoi(word));
	  
	  offset = 15;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAmin = atoi(word);
	  //printf("RAmin = \"%s\" %d\n", word, atoi(word));
	  
	  offset = 18;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAsec = atoi(word);
	  //printf("RAsec = \"%s\" %d\n", word, atoi(word));
	  RAdegf = 360.0*(1.0*RAhr + 1.0*RAmin/60.0 + 1.0*RAsec/3600.0)/24.0;
	  
	  
	  switch(line[21]){
	  case '+':
	    sign = 1.0;
	    break;
	  case '-':
	    sign = -1.0;
	    break;
	  default:
	    printf("Sign char not recognized: %c\n", line[21]);
	    exit(0);
	    break;
	  }

	  offset = 22;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECdeg = atoi(word);
	  //printf("DECdeg = \"%s\" %d\n", word, atoi(word));
	  
	  offset = 25;
	  for(i=0;i<4;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECmin = atof(word);
	  // printf("DECmin = \"%s\" %f\n", word, atof(word));
	  DECdegf = sign*(1.0*DECdeg + 1.0*DECmin/60.0);
	  
	  //    printf("%ld %f %f ", k, RAdegf, DECdegf);
	  
	  RA1950 = RAdegf/180.0*M_PI;
	  DEC1950 = DECdegf/180.0*M_PI;
	  B1950t0J2000();
	  
	  CNS3star[k].RAdeg = RA2000/M_PI*180.0;
	  CNS3star[k].DECdeg = DEC2000/M_PI*180.0;
	  
	  
	  // SPECTRAL TYPE
	  offset = 55;
	  for(i=0;i<10;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  
	  CNS3star[k].SpType = -100.0; // unknown
	  if(word[0]=='D')
	    CNS3star[k].SpType = 0.0;
	  if((word[0]==' ')||(word[0]=='d'))
	    {
	      switch (word[1]) {
	      case 'O':
		CNS3star[k].SpType = 1.0;
		break;
	      case 'B':
		CNS3star[k].SpType = 2.0;
		break;
	      case 'A':
		CNS3star[k].SpType = 3.0;
		break;
	      case 'F':
		CNS3star[k].SpType = 4.0;
		break;
	      case 'G':
		CNS3star[k].SpType = 5.0;
		break;
	      case 'K':
		CNS3star[k].SpType = 6.0;
		break;
	      case 'M':
		CNS3star[k].SpType = 7.0;
		break;
	      }
	    }
	  
	  if(CNS3star[k].SpType > 0.01) // not a white dwarf
	    {
	      
	      if(isdigit(word[2])!=0)
		{
		  word1[0] = word[2];
		  word1[1] = '\0';
		  iend = 3;
		  if(word[3]=='.')
		    {
		      word1[1] = '.';
		      word1[2] = word[4];
		      word1[3] = '\0';
		      iend = 5;
		    }
		  val1 = atof(word1);
		  val2 = val1;
		  
		  if(word[iend]=='-')
		    {
		      word1[0] = word[iend+1];
		      word1[1] = '\0';
		      if(word[iend+2]=='.')
			{
			  word1[1] = '.';
			  word1[2] = word[iend+3];
			  word1[3] = '\0';
			}
		      val2 = atof(word1);
		    }
		  CNS3star[k].SpType += 0.1*(val1+val2)*0.5;
		}
	    }
	  
	  
	  
	  
	  offset = 67;
	  for(i=0;i<6;i++)
	    word[i] = line[i+offset];
	  word[6] = '\0';
	  if(strcmp(word,"      ")==0)
	    CNS3star[k].mV = -100.0;
	  else
	    CNS3star[k].mV = atof(word);

	  //	  if((CNS3star[k].mV<0.1)&&(CNS3star[k].mV>-0.1))
	  //printf("%s %s %s *********** \"%s\"  -> %f\n", CNS3star[k].name1, CNS3star[k].name2, CNS3star[k].name3, word, CNS3star[k].mV);


      
	  offset = 75;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[5] = '\0';
	  if((word[0]=='+')||(word[0]=='-'))
	    CNS3star[k].BV = atof(word);
	  else
	    CNS3star[k].BV = -100.0;
	  
	  
	  offset = 82;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[5] = '\0';
	  if((word[0]=='+')||(word[0]=='-'))
	    CNS3star[k].UB = atof(word);
	  else
	    CNS3star[k].UB = -100.0;
	  
	  
	  
	  offset = 89;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[5] = '\0';
	  if((word[0]=='+')||(word[0]=='-'))
	    CNS3star[k].RI = atof(word);
	  else
	    CNS3star[k].RI = -100.0;      
	  
	  
	  offset = 121;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[5] = '\0';
	  if(strcmp(word,"     ")==0)
	    CNS3star[k].MV = -100.0;
	  else
	    CNS3star[k].MV = atof(word);

	  
	  if((CNS3star[k].nameB!=1)&&(CNS3star[k].nameC!=1))	  
	    {
	      k++;
	    }
	}
    }
 
  fclose(fp);
    
  NBstarCNS3 = k;
  printf("%ld sources\n", NBstarCNS3);
  //  exit(0);
  return(0);
}



int read_CNS3_2MASS()
{
  // http://heasarc.gsfc.nasa.gov/W3Browse/star-catalog/gliese2mas.html
  
  FILE *fp;
  long k;
  char line[SBUFFERSIZE];
  long i, j, offset;
  char word[200];
	char fname[200];

  int RAhr, RAmin, DECdeg, DECmin;
  float RAsec, DECsec;
  float sign;
  int namecnt;

  printf("READING CNS3_2MASS ... ");

//  fp = fopen("Gliese2MASS_clean.dat", "r");  
	
	sprintf(fname, "Gliese2MASS_clean.dat");	
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}

  k = 0;
  while(fgets(line,SBUFFERSIZE,fp)!=NULL)
    {
      if(line[0]=='|')
	{
	  offset = 1;
	  for(i=0;i<12;i++)
	    CNS3_2MASS_star[k].name[i] = line[i+offset];
	  CNS3_2MASS_star[k].name[12] = '\0';

	  namecnt = 0;
	  i = 0;
	  j = 0;
	  while((i<12)&&(namecnt<3))
	    {
	      if(CNS3_2MASS_star[k].name[i] == ' ')
		{
		  switch(namecnt){
		  case 0:
		    CNS3_2MASS_star[k].name1[j] = '\0';
		    break;
		  case 1:
		    CNS3_2MASS_star[k].name2[j] = '\0';
		    break;
		  case 2:
		    CNS3_2MASS_star[k].name3[j] = '\0';
		    break;
		  }
		  namecnt++;
		  j = 0;
		}
	      else
		{
		  switch(namecnt){
		  case 0:
		    CNS3_2MASS_star[k].name1[j] = CNS3_2MASS_star[k].name[i];
		    break;
		  case 1:
		    CNS3_2MASS_star[k].name2[j] = CNS3_2MASS_star[k].name[i];
		    break;
		  case 2:
		    CNS3_2MASS_star[k].name3[j] = CNS3_2MASS_star[k].name[i];
		    break;
		  default:
		    namecnt++;
		    break;
		  }
		  j++;
		}
	      i++;
	    }
	  
	  //  printf("%s -> %s %s %s\n", CNS3_2MASS_star[k].name, CNS3_2MASS_star[k].name1, CNS3_2MASS_star[k].name2, CNS3_2MASS_star[k].name3);

	  offset = 32;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAhr = atoi(word);
	  //printf("RAhr = \"%s\" %d\n", word, atoi(word));
    
	  offset = 35;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAmin = atoi(word);
	  //printf("RAmin = \"%s\" %d\n", word, atoi(word));
      
	  offset = 38;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAsec = atof(word);
	  //printf("RAsec = \"%s\" %d\n", word, atoi(word));
	  CNS3_2MASS_star[k].RAdeg = 360.0*(1.0*RAhr + 1.0*RAmin/60.0 + 1.0*RAsec/3600.0)/24.0;
      
	  	 
	  if(line[44]=='+')
	    sign = 1.0;
	  else
	    sign = -1.0;
	  offset = 45;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECdeg = atoi(word);
	  //printf("DECdeg = \"%s\" %d\n", word, atoi(word));
    
	  offset = 48;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECmin = atof(word);
	  // printf("DECmin = \"%s\" %f\n", word, atof(word));
    
	  offset = 51;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECsec = atof(word);
	  //printf("RAsec = \"%s\" %d\n", word, atoi(word));

	  CNS3_2MASS_star[k].DECdeg = sign*(1.0*DECdeg + 1.0*DECmin/60.0 + 1.0*DECsec/3600);

	  //    printf("%ld %f %f ", k, RAdegf, DECdegf);
    

    
    
	  offset = 80;
	  for(i=0;i<6;i++)
	    word[i] = line[i+offset];
	  word[6] = '\0';
	  if(strcmp(word,"      ")==0)
	    CNS3_2MASS_star[k].Jmag = -100.0;
	  else
	    CNS3_2MASS_star[k].Jmag = atof(word);
    
	  offset = 87;
	  for(i=0;i<6;i++)
	    word[i] = line[i+offset];
	  word[6] = '\0';
	  if(strcmp(word,"      ")==0)
	    CNS3_2MASS_star[k].Hmag = -100.0;
	  else
	    CNS3_2MASS_star[k].Hmag = atof(word);
    
	  offset = 94;
	  for(i=0;i<6;i++)
	    word[i] = line[i+offset];
	  word[6] = '\0';
	  if(strcmp(word,"      ")==0)
	    CNS3_2MASS_star[k].Kmag = -100.0;
	  else
	    CNS3_2MASS_star[k].Kmag = atof(word);
    
	  //if((CNS3_2MASS_star[k].name3[0]!='B')&&(CNS3_2MASS_star[k].name3[0]!='C'))
	  //{
	  star[k].FLAG_2MASSCNS3 = 1;
	  k++;
	      //  }
	}
    }
  fclose(fp);
  NBstar_CNS3_2MASS = k;

  printf("%ld sources\n", NBstar_CNS3_2MASS);

  return(0);
}



int read_HIP()
{
  FILE *fp;
  FILE *fp2;
  long k, cnt;
  char word[200];
  char word1[200];
  long i, iend, offset;
  char line[SBUFFERSIZE];
  char line2[SBUFFERSIZE];
  double val, val1, val2;
	char fname[200];

  long index2;
  long cnt2;
  int OK2;
  float BV2, VI2, Hpmag2, dist2;
  

  printf("READING HIPPARCOS ... ");
  fflush(stdout);

  sprintf(fname, "hip_main.dat");
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  


  sprintf(fname, "hip_main2.dat"); // I/311
	fp2 = fopen(fname, "r");
	if(fp2==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  


  k = 0;

  cnt = 0;
  while(cnt<101)
    {
      fgets(line2,SBUFFERSIZE,fp2);

      offset = 29;
      for(i=0;i<6;i++)
	word[i] = line2[i+offset];
      word[i] = '\0';
      index2 = atoi(word);
      //      printf("index2 = %ld\n", index2);
      cnt++;
    }
  //  printf("%s\n",line2);

  //  exit(0);

  while(fgets(line,SBUFFERSIZE,fp)!=NULL)
    {
      if((line[0]=='H')&&(line[1]='|'))
	{


	  offset = 8;
	  for(i=0;i<6;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  HIPstar[k].index = atoi(word);
	  
	  cnt2 = 0;
	  while((index2!=HIPstar[k].index)&&(cnt2<2))
	    {
	      if(fgets(line2,SBUFFERSIZE,fp2)==NULL)
		{
		  cnt2++;
		  fclose(fp2);
		  
//		  fp2 = fopen("hip_main2.dat", "r");
		    sprintf(fname, "hip_main2.dat"); 
			fp2 = fopen(fname, "r");
			if(fp2==NULL)
			{
				printf("ERROR: File \"%s\" missing\n", fname);
				exit(0);
			}
  
		  
		  cnt = 0;
		  while(cnt<101)
		    {
		      fgets(line2,SBUFFERSIZE,fp2);
		      
		      offset = 29;
		      for(i=0;i<6;i++)
			word[i] = line2[i+offset];
		      word[i] = '\0';
		      index2 = atoi(word);
		      //      printf("index2 = %ld\n", index2);
		      cnt++;
		    }
		}
	      offset = 29;
	      for(i=0;i<6;i++)
		word[i] = line2[i+offset];
	      word[i] = '\0';
	      index2 = atoi(word);	  
	      
	    }
	  //	  printf("--------- %ld %ld -----------\n", HIPstar[k].index, index2);
	  if(index2!=HIPstar[k].index)
	    {
	      OK2 = 0;
	      //	      printf("WARNING: HIP %ld not found\n", HIPstar[k].index);
	    }
	  else
	    {	      	      
	      offset = 86;
	      for(i=0;i<7;i++)
		word[i] = line2[i+offset];
	      word[i] = '\0';
	      dist2 = 1000.0/atof(word);
	      if(val<0.0)
		dist2 = 10000.0;
	     
	      offset = 153;
	      for(i=0;i<7;i++)
		word[i] = line2[i+offset];
	      word[i] = '\0';
	      Hpmag2 = atof(word);
	      
	      offset = 178;
	      for(i=0;i<6;i++)
		word[i] = line2[i+offset];
	      word[i] = '\0';
	      BV2 = atof(word);
	      
	      offset = 192;
	      for(i=0;i<6;i++)
		word[i] = line2[i+offset];
	      word[i] = '\0';
	      VI2 = atof(word);
	      
	      OK2 = 1;
	    }


	  if(((HIPstar[k].index == 114110)||(HIPstar[k].index == 82725)||(HIPstar[k].index == 114176))&&(OK2==0)) // artefacts
	    {
	    }
	  else
	    {
	      offset = 79;
	      for(i=0;i<7;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      val = 1000.0/atof(word);
	      if(val<0.0)
		val = 10000.0;
	      
	      
	      /*	      if(HIPstar[k].index==50)
		{
		  printf("read_HIP: HIP 15689 is here - par: \"%s\" dist = %f\n", word, val);
		  exit(0);
		  }*/
	      
	      
	      
	      // if(val<DISTMAX)
	      //{
	      HIPstar[k].dist = dist2; //val;
	      //	      printf("%f -> %f\n", val, dist2);
	      
	      offset = 8;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPstar[k].index = atoi(word);
	      
	      offset = 41;
	      for(i=0;i<5;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPstar[k].Vmag = atof(word);
	      

	      offset = 51;
	      for(i=0;i<12;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPstar[k].RAdeg = atof(word);
	      
	      offset = 64;
	      for(i=0;i<12;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPstar[k].DECdeg = atof(word);
	      
	    

	      offset = 274;
	      for(i=0;i<7;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPstar[k].Hpmag = Hpmag2; //atof(word);
	      

	      
	      offset = 217;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		HIPstar[k].BT = -100.0;
	      else
		HIPstar[k].BT = atof(word);
	      

	      offset = 230;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		HIPstar[k].VT = -100.0;
	      else
		HIPstar[k].VT = atof(word);
	      
 
	      offset = 245;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		HIPstar[k].BV = -100.0;
	      else
		HIPstar[k].BV = atof(word);
	      if(fabs(BV2)>0.0001)
		HIPstar[k].BV = BV2;


	      offset = 260;
	      for(i=0;i<4;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"    ")==0)
		HIPstar[k].VI = -100.0;
	      else
		HIPstar[k].VI = atof(word);
	      /*    if(HIPstar[k].index == 67593)
		{
		  printf("---------- %f %f -------- \n", HIPstar[k].VI , VI2);
		  exit(0);
		  }*/		
	      if(fabs(VI2)>0.0001)
		HIPstar[k].VI = VI2;
	      
	      
	      HIPstar[k].SpType = -1; // unknown
	       offset = 435;
	      for(i=0;i<12;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	     
	      switch (word[0]) {
	      case 'O':
		HIPstar[k].SpType = 1.0;
		break;
	      case 'B':
		HIPstar[k].SpType = 2.0;
		break;
	      case 'A':
		HIPstar[k].SpType = 3.0;
		break;
	      case 'F':
		HIPstar[k].SpType = 4.0;
		break;
	      case 'G':
		HIPstar[k].SpType = 5.0;
		break;
	      case 'K':
		HIPstar[k].SpType = 6.0;
		break;
	      case 'M':
		HIPstar[k].SpType = 7.0;
		break;
	      }
	       
	      
	      if(HIPstar[k].SpType > 0.0)
		{
		  
		  if(isdigit(word[1])!=0)
		    {
		      word1[0] = word[1];
		      word1[1] = '\0';
		      iend = 3;
		      if(word[2]=='.')
			{
			  word1[1] = '.';
			  word1[2] = word[3];
			  word1[3] = '\0';
			  iend = 4;
			}
		      val1 = atof(word1);
		      val2 = val1;
		      
		      if(word[iend]=='-')
			{
			  word1[0] = word[iend+1];
			  word1[1] = '\0';
			  if(word[iend+2]=='.')
			    {
			      word1[1] = '.';
			      word1[2] = word[iend+3];
			      word1[3] = '\0';
			    }
			  val2 = atof(word1);
			}
		      HIPstar[k].SpType += 0.1*(val1+val2)*0.5;
		    }
		}
	      
	      k ++;
	      //  }
	    }
	}      
    }
  NBstarHIP = k;    

  printf(" %ld stars\n", NBstarHIP);
  fflush(stdout);

  fclose(fp);
  fclose(fp2);

  return(0);
}



int read_HIPD()
{
  FILE *fp;
  long k;
  char word[200];
  long i, offset;
  char line[SBUFFERSIZE];
  double val;
	char fname[200];
  long lcnt = 0;


  sprintf(fname, "hip_doubles.dat");
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  


  k = 0;
  while((fgets(line,SBUFFERSIZE,fp)!=NULL))
    {
      //      printf("%ld line : \"%s\"\n", lcnt, line);
      if((lcnt>121)&&(line[0]!='#'))
	{
	  //	   printf("%ld line : \"%s\"\n", lcnt, line);
	   offset = 143;
	  for(i=0;i<7;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  val = 1000.0/atof(word);
	  if(val<0.0)
	    val = 10000.0;
	  //printf("dist = %f pc\n", val);
	  
	  //	  if((val<DISTMAX)&&(val>0.1))
	  //  {
	      HIPDstar[k].dist = val;
	      //	      printf("---------------------------------------------------------------------\n");
	      //  printf("\"%s\"\n",line);
	      // printf("%ld %ld dist = %f pc\n", lcnt, k, val);
	      
	      
	      offset = 68;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPDstar[k].index = atoi(word);
 
	      switch (line[66]) {
	      case 'A':
		HIPDstar[k].component = 0;
		break;
	      case 'S':
		HIPDstar[k].component = 1;
		break;
	      case 'P':
		HIPDstar[k].component = 8;
		break;
	      case 'E':
		HIPDstar[k].component = 8;
		break;
	      case 'R':
		HIPDstar[k].component = 8;
		break;
	      case 'K':
		HIPDstar[k].component = 8;
		break;
	      case 'F':
		HIPDstar[k].component = 8;
		break;
	      case 'G':
		HIPDstar[k].component = 8;
		break;
	      case 'H':
		HIPDstar[k].component = 8;
		break;
	      case 'N':
		HIPDstar[k].component = 8;
		break;
	      case 'B':
		HIPDstar[k].component = 1;
		break;
	      case 'C':
		HIPDstar[k].component = 2;
		break;
	      case 'D':
		HIPDstar[k].component = 3;
		break;
	      default:
		printf("ERROR: %ld: char 66 not recognized: %c\n", k, line[66]);
		printf("%s\n", line);
		exit(0);
		break;
	      }
		
	      offset = 75;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPDstar[k].Hpmag = atof(word);
	      

	      offset = 117;
	      for(i=0;i<12;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPDstar[k].RAdeg = atof(word);
	      
	      offset = 130;
	      for(i=0;i<12;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      HIPDstar[k].DECdeg = atof(word);
	      
	    
	      
	      offset = 89;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		HIPDstar[k].BT = -100.0;
	      else
		HIPDstar[k].BT = atof(word);
	      

	      offset = 103;
	      for(i=0;i<6;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		HIPDstar[k].VT = -100.0;
	      else
		HIPDstar[k].VT = atof(word);
	      
 
	      offset = 215;
	      for(i=0;i<7;i++)
		word[i] = line[i+offset];
	      word[i] = '\0';
	      if(strcmp(word,"      ")==0)
		HIPDstar[k].rho = -100.0;
	      else
		HIPDstar[k].rho = atof(word);
	      
	  
	      //	      printf("%d %d %f %f %f %f %f %f\n", HIPDstar[k].index, HIPDstar[k].component, HIPDstar[k].RAdeg, HIPDstar[k].DECdeg, HIPDstar[k].Hpmag, HIPDstar[k].BT, HIPDstar[k].VT, HIPDstar[k].rho);

     	      k ++;
	      //    }
	}      
      lcnt ++;
    }
  NBstarHIPD = k;    

  printf("HIPPARCOS doubles: %ld stars within %f pc\n", NBstarHIPD, DISTMAX);

  fclose(fp);


  return(0);
}




int read_RECONS()
{
  FILE *fp;
  char line[200];
  char fname[200];
  int offset;
  int i, iend;
  long k;
  char word[200];
  char word1[200];
  float sign;
  int RAhr, RAmin;
  float RAsec;
  int DECdeg, DECmin, DECsec;
  int selOK;
  /*
    typedef struct
    {
    char name1[200];
    char name2[200];
    char name3[200];
    float RAdeg;
    float DECdeg;
    float SpType;
    float mV;
    float MV;
    float mass;
    } RECONSSTAR;
  */

  printf("READING RECONS ... ");
  
//  fp = fopen("RECONS100.dat","r");
  sprintf(fname, "RECONS100.dat"); 
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  


  k = 0;
  while((fgets(line,SBUFFERSIZE,fp)!=NULL))
    {					      
      selOK = 1;
      if(line[3]=='.')
	{
	  offset = 5;
	  for(i=0;i<16;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  sscanf(word,"%s %s %s", RECONSstar[k].name1, RECONSstar[k].name2,  RECONSstar[k].name3 );


	  
	  offset = 34;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAhr = atoi(word);
		      
	  offset = 37;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAmin = atoi(word);
	  
	  offset = 40;
	  for(i=0;i<4;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RAsec = atof(word);
	  
	  RECONSstar[k].RAdeg = 360.0*(1.0*RAhr + 1.0*RAmin/60.0 + 1.0*RAsec/3600.0)/24.0;
	      
	      
	  if(line[45]=='+')
	    sign = 1.0;
	  else
	    sign = -1.0;
	  offset = 46;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECdeg = atoi(word);
	  
	  offset = 49;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECmin = atoi(word);
	  
	  offset = 52;
	  for(i=0;i<2;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  DECsec = atof(word);
	  	      
	  RECONSstar[k].DECdeg = sign*(1.0*DECdeg + 1.0*DECmin/60.0 + 1.0*DECsec/3600.0);
	  

	  offset = 75;
	  for(i=0;i<7;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RECONSstar[k].dist = 1.0/atof(word);


	  // SPECTRAL TYPE
	  offset = 98;
	  for(i=0;i<8;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  
	  RECONSstar[k].SpType = -100.0; // unknown
	  if(word[0]=='D')
	    RECONSstar[k].SpType = 0.0; // White dwarf
	  if((word[0]==' '))
	    {
	      switch (word[1]) {
	      case 'O':
		RECONSstar[k].SpType = 1.0;
		break;
	      case 'B':
		RECONSstar[k].SpType = 2.0;
		break;
	      case 'A':
		RECONSstar[k].SpType = 3.0;
		break;
	      case 'F':
		RECONSstar[k].SpType = 4.0;
		break;
	      case 'G':
		RECONSstar[k].SpType = 5.0;
		break;
	      case 'K':
		RECONSstar[k].SpType = 6.0;
		break;
	      case 'M':
		RECONSstar[k].SpType = 7.0;
		break;
	      default:
		//	printf("NOT SELECTED\n");
		selOK = 0;
		break;
	      }
	    }
	  
	  if(RECONSstar[k].SpType > 0.1) // not a white dwarf
	    {
	      if(isdigit(word[2])!=0)
		{
		  word1[0] = word[2];
		  word1[1] = '\0';
		  iend = 3;
		  if(word[3]=='.')
		    {
		      word1[1] = '.';
		      word1[2] = word[4];
		      word1[3] = '\0';
		      iend = 5;
		    }
		  RECONSstar[k].SpType += 0.1*atof(word1);
		}	      	  
	    }
	
	  
	  offset = 111;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RECONSstar[k].mV = atof(word);

	  offset = 119;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RECONSstar[k].MV = atof(word);

	  offset = 126;
	  for(i=0;i<5;i++)
	    word[i] = line[i+offset];
	  word[i] = '\0';
	  RECONSstar[k].mass = atof(word);


	  /*	  if(selOK == 1)
	    {
	      printf("%s", line);
	      printf("%s %s %s  %f %f %f %f %f %f\n", RECONSstar[k].name1, RECONSstar[k].name2,  RECONSstar[k].name3, RECONSstar[k].RAdeg, RECONSstar[k].DECdeg, RECONSstar[k].SpType, RECONSstar[k].mV, RECONSstar[k].MV, RECONSstar[k].mass );
	      printf("\n");
	      }*/
	  if(selOK == 1)
	    k++;
	}
    }
  fclose(fp);

  NBstarRECONS = k;    
  printf("%ld sources\n", NBstarRECONS);
  return(0);
}




int MatchCatalogs()
{
  // CATALOG 1 (MAIN) : SUPERBLINK
  // CATALOG 2 : CNS3
  // CATALOG 3 : Gliese 2MASS
  
  FILE *fp;
  FILE *fp1;
  
  long k1, k2, k3, k2m;
  double dist, distmin;
  long k1min;

  long cntOK = 0;
  int m1, m2, m3;
  int mOK;
  long mOKcnt = 0;
  int mmOK, k2mm;
  
  long CNS3match;
  long HIPmatch;
  long HIPDmatch;
  long RECONSmatch;

  long NBstar0;

  long NBadd;
  long l2;

  long cnt1, cnt2;
  int matchOK;


  printf("SUPERBLINK ONLY : %ld\n", NBstar);

  //
  // ADD Gliese 2MASS info to SUPERBLINK
  //
  cnt1 = 0;
  cnt2 = 0;
  NBstar0 = NBstar;
  fp1 = fopen("MatchLog.txt","w");
  for(k3=0;k3<NBstar_CNS3_2MASS;k3++)
    {
      RA2000 = CNS3_2MASS_star[k3].RAdeg/180.0*M_PI;
      DEC2000 = CNS3_2MASS_star[k3].DECdeg/180.0*M_PI;

      // start from CNS3_2MASS, find corresponding entry in CNS3
      k2mm = -1;
      k2 = 0;
      mOK = 1;
      while((k2<NBstarCNS3)&&(mOK!=0))
	{
	  m1 = strcmp(CNS3_2MASS_star[k3].name1, CNS3star[k2].name1);
	  m2 = strcmp(CNS3_2MASS_star[k3].name2, CNS3star[k2].name2);
	  m3 = strcmp(CNS3_2MASS_star[k3].name3, CNS3star[k2].name3);

	  mOK = m2*m2;
	  if((mOK==0)&&(CNS3_2MASS_star[k3].name3[0]==CNS3star[k2].name3[0]))
	    mOK = 0;
	  if((mOK==0)&&(CNS3_2MASS_star[k3].name3[0]=='A')&&(CNS3star[k2].name3[0]=='\0'))
	    mOK = 0;
	  if((mOK==0)&&(CNS3_2MASS_star[k3].name3[0]=='\0')&&(CNS3star[k2].name3[0]=='A'))
	    mOK = 0;
	  if((mOK==0)&&(CNS3_2MASS_star[k3].name3[0]=='A')&&(CNS3star[k2].name3[0]==' '))
	    mOK = 0;
	  if((mOK==0)&&(CNS3_2MASS_star[k3].name3[0]==' ')&&(CNS3star[k2].name3[0]=='A'))
	    mOK = 0;
	  
	  mmOK = m2*m2;
	  if(mmOK==0)
	    k2mm = k2;
	  if(mOK==0)
	    k2m = k2;
	  k2++;	  	
	}
      if((mOK!=0)&&(k2mm>-1))
	{
	  mOK = 0;
	  k2m = k2mm;
	  CNS3match = 1;
	}

      if(mOK!=0)
	{
	  //	  printf("No Match found for CNS3_2MASS %s %s %s in CNS3 -> unknown distance\n", CNS3_2MASS_star[k3].name1, CNS3_2MASS_star[k3].name2, CNS3_2MASS_star[k3].name3);
	  CNS3match = 0;
	  mOKcnt ++;
	}


      // k3: index in CNS3_2MASS
      // k2m: if CNS3match == 1, entry in CNS3
      
      
      // look for corresponding entry in SUPERBLINK

      // Try to match names
      k1min = 0;
      matchOK = 0;
      for(k1=0;k1<NBstar0;k1++)
	{
	  //	  printf("TEST : %d \"%s\" \"%s\" \"%s\"  <->  %ld \"%s\" \"%s\" \"%s\"\n", k3, CNS3star[k3].name1, CNS3star[k3].name2, CNS3star[k3].name3, k1, star[k1].SUPERBLINKCNS3name1,  star[k1].SUPERBLINKCNS3name2,  star[k1].SUPERBLINKCNS3name3);
	  //fflush(stdout);
	  m2 = strcmp(star[k1].SUPERBLINKCNS3name2, CNS3_2MASS_star[k3].name2);
	  m3 = strcmp(star[k1].SUPERBLINKCNS3name3, CNS3_2MASS_star[k3].name3);
	  l2 = strlen(star[k1].SUPERBLINKCNS3name2);
	  if(star[k1].SUPERBLINKCNS3name3[0] == CNS3_2MASS_star[k3].name3[0])
	    m3 = 0;
	  if((star[k1].SUPERBLINKCNS3name3[0] == 'A')&&(CNS3_2MASS_star[k3].name3[0]=='\0'))
	    m3 = 0;
	  if((star[k1].SUPERBLINKCNS3name3[0] == '\0')&&(CNS3_2MASS_star[k3].name3[0]=='A'))
	    m3 = 0;
	  if((star[k1].SUPERBLINKCNS3name3[0] == 'A')&&(CNS3_2MASS_star[k3].name3[0]==' '))
	    m3 = 0;
	  if((star[k1].SUPERBLINKCNS3name3[0] == ' ')&&(CNS3_2MASS_star[k3].name3[0]=='A'))
	    m3 = 0;
	  
	 
	  if((m2*m2+m3*m3==0)&&(l2>0))
	    {
	      k1min = k1;
	      matchOK = 1;


	      printf("NAME MATCHED : %ld  %s %s %s  <->  %ld  %s %s %s\n", k3, CNS3_2MASS_star[k3].name1, CNS3_2MASS_star[k3].name2, CNS3_2MASS_star[k3].name3, k1, star[k1].SUPERBLINKCNS3name1,  star[k1].SUPERBLINKCNS3name2,  star[k1].SUPERBLINKCNS3name3);
	      cnt1 ++;
	    }
	}

      if(matchOK==0) // no match in name, try coordinate matching
	{
	  distmin = 100.0;
	  for(k1=0;k1<NBstar0;k1++)
	    {
	      dist = acos( sin(DEC2000)*sin(star[k1].DECdeg/180.0*M_PI)+cos(DEC2000)*cos(star[k1].DECdeg/180.0*M_PI)*cos(RA2000-star[k1].RAdeg/180.0*M_PI) );
	      dist *= 180.0/M_PI; // dist in degree
	      
	      if((dist<distmin)&&(m3==0)) // stars are close and same component identifier
		if(fabs(star[k1].Jmag-CNS3_2MASS_star[k3].Jmag)<0.1) // if SUPERBLINK Jmag and 2MASS Jmag are same
		  {
		    distmin = dist;
		    k1min = k1;
		  }	      
	    }
	  if(distmin < 0.05) // within 3 acrmin
	    {	     
	      printf("COORD MATCHED : %ld  %s %s %s  %f %f <->  %ld  %s %s %s  %f %f    -> %f %f\n", k3, CNS3_2MASS_star[k3].name1, CNS3_2MASS_star[k3].name2, CNS3_2MASS_star[k3].name3, CNS3_2MASS_star[k3].RAdeg, CNS3_2MASS_star[k3].DECdeg, k1min, star[k1min].SUPERBLINKCNS3name1,  star[k1min].SUPERBLINKCNS3name2,  star[k1min].SUPERBLINKCNS3name3, star[k1min].RAdeg, star[k1min].DECdeg, distmin, star[k1min].Jmag-CNS3_2MASS_star[k3].Jmag);
	      matchOK = 1;
	    }
	}
      

    
      if(matchOK==1) // MATCH FOUND - DO NOT ADD NEW TARGET - update fields
	{
	  cnt2 ++;
	  k1 = k1min;
	  //	  printf("COORD MATCHED : %ld  %s %s %s  <->  %ld  %s %s %s\n", k3, CNS3_2MASS_star[k3].name1, CNS3_2MASS_star[k3].name2, CNS3_2MASS_star[k3].name3, k1, star[k1].SUPERBLINKCNS3name1,  star[k1].SUPERBLINKCNS3name2,  star[k1].SUPERBLINKCNS3name3);

	  star[k1].FLAG_2MASSCNS3 = 1;
	  strcpy(star[k1].name1, CNS3_2MASS_star[k3].name1);
	  strcpy(star[k1].name2, CNS3_2MASS_star[k3].name2);
	  strcpy(star[k1].name3, CNS3_2MASS_star[k3].name3);
	  
	  star[k1].RAdeg = CNS3_2MASS_star[k3].RAdeg;
	  star[k1].DECdeg = CNS3_2MASS_star[k3].DECdeg;
	  	      
	  star[k1].Jmag = CNS3_2MASS_star[k3].Jmag; // 2MASS
	  star[k1].Hmag = CNS3_2MASS_star[k3].Hmag; // 2MASS
	  star[k1].Kmag = CNS3_2MASS_star[k3].Kmag; // 2MASS	  
	  if(CNS3match==1)
	    {
	      star[k1].FLAG_CNS3 = 1;
	      star[k1].dist = CNS3star[k2m].dist;
	      
	      //	      printf("MATCH : %d \"%s\" \"%s\" \"%s\"  <->  %ld %s %s %s [%f]\n", k2m, CNS3star[k2m].name1, CNS3star[k2m].name2, CNS3star[k2m].name3, k1, star[k1].SUPERBLINKCNS3name1,  star[k1].SUPERBLINKCNS3name2,  star[k1].SUPERBLINKCNS3name3, mmin);
	      
	      strcpy(star[k1].CNS3name1, CNS3star[k2m].name1);
	      strcpy(star[k1].CNS3name2, CNS3star[k2m].name2);
	      strcpy(star[k1].CNS3name3, CNS3star[k2m].name3);

	      star[k1].SpType = CNS3star[k2m].SpType; // from CNS3 ... -1 if unknown
	      star[k1].mV = CNS3star[k2m].mV; // from CNS3
	      star[k1].BV = CNS3star[k2m].BV; // from CNS3
	      star[k1].UB = CNS3star[k2m].UB; // from CNS3
	      star[k1].RI = CNS3star[k2m].RI; // from CNS3
	      star[k1].CNS3dist = CNS3star[k2m].dist; // from CNS3
	      star[k1].MV = CNS3star[k2m].MV; // from CNS3
	    } 
	  cntOK ++;
	}
      else // NO MATCH FOUND IN SUPERBLINK - ADD NEW TARGET
	{
	  k1 = NBstar;
	  star[k1].FLAG_2MASSCNS3 = 1;
	  strcpy(star[k1].name1, CNS3_2MASS_star[k3].name1);
	  strcpy(star[k1].name2, CNS3_2MASS_star[k3].name2);
	  strcpy(star[k1].name3, CNS3_2MASS_star[k3].name3);

	  star[k1].RAdeg = CNS3_2MASS_star[k3].RAdeg;
	  star[k1].DECdeg = CNS3_2MASS_star[k3].DECdeg;
	  	      
	  star[k1].Jmag = CNS3_2MASS_star[k3].Jmag; // 2MASS
	  star[k1].Hmag = CNS3_2MASS_star[k3].Hmag; // 2MASS
	  star[k1].Kmag = CNS3_2MASS_star[k3].Kmag; // 2MASS
	      
	  if(CNS3match==1)
	    {
	      star[k1].FLAG_CNS3 = 1;
	      star[k1].dist = CNS3star[k2m].dist;
	      
	      strcpy(star[k1].CNS3name1, CNS3star[k2m].name1);
	      strcpy(star[k1].CNS3name2, CNS3star[k2m].name2);
	      strcpy(star[k1].CNS3name3, CNS3star[k2m].name3);

	      star[k1].SpType = CNS3star[k2m].SpType; // from CNS3 ... -1 if unknown
	      star[k1].mV = CNS3star[k2m].mV; // from CNS3
	      star[k1].BV = CNS3star[k2m].BV; // from CNS3
	      star[k1].UB = CNS3star[k2m].UB; // from CNS3
	      star[k1].RI = CNS3star[k2m].RI; // from CNS3
	      star[k1].CNS3dist = CNS3star[k2m].dist; // from CNS3
	      star[k1].MV = CNS3star[k2m].MV; // from CNS3
	    }
	  NBstar++;
	  // fprintf(fp1,"%ld %ld %g %g %g\n", k2, k1min, distmin, CNS3star[k2m].dist, star[k1min].dist);
	}
    }
  fclose(fp1);
  printf("2MASS/CNS3 targets found in SUPERBLINK : %ld\n", cnt1);

  printf("2MASS not matched with original CNS3: %ld / %ld\n", mOKcnt, NBstar_CNS3_2MASS);

  printf("2MASS_CNS3 [%ld] + SUPERBLINK [%ld] -> %ld\n", NBstar_CNS3_2MASS, NBstar0, NBstar);

  printf("%ld out of %ld CNS3+2MASS targets were already in SUPERBLINK\n", cntOK, NBstar_CNS3_2MASS);



  // ADDING HIPPARCOS ENTRIES
  printf("\nAdding Hipparcos singles ...\n");
  fp = fopen("HIPmatch.txt","w");
  HIPmatch = 0;
  NBstar0 = NBstar;
  for(k2=0;k2<NBstarHIP;k2++) // scan HIPPARCOS catalog
    {
      printf("\r[%ld / %ld]   ", k2, NBstarHIP);
      fflush(stdout);
      RA2000 = HIPstar[k2].RAdeg/180.0*M_PI;
      DEC2000 = HIPstar[k2].DECdeg/180.0*M_PI;
      
      distmin = 1000.0; 
      for(k1=0;k1<NBstar0;k1++)
	{
	  dist = acos( sin(DEC2000)*sin(star[k1].DECdeg/180.0*M_PI)+cos(DEC2000)*cos(star[k1].DECdeg/180.0*M_PI)*cos(RA2000-star[k1].RAdeg/180.0*M_PI) );
	  dist *= 180.0/M_PI;
	  if(dist<distmin)
	    {
	      distmin = dist;
	      k1min = k1;
	    }	  
	}
      /*      if(HIPstar[k2].index==29271)
	{
	  printf("HIP 29271 is here\n");
	  printf("distmin = %f %ld\n", distmin, k1min);
	  exit(0);
	  }*/
      
      if(distmin<0.03) // UPDATE EXISTING ENTRY
	{
	  star[k1min].RAdeg = HIPstar[k2].RAdeg;
	  star[k1min].DECdeg = HIPstar[k2].DECdeg;
	  star[k1min].HIPindex = HIPstar[k2].index;
	  star[k1min].dist = HIPstar[k2].dist; // HIPPARCOS TAKES PRIORITY ON PARALLAX
	  star[k1min].HIPVmag = HIPstar[k2].Vmag;
	  star[k1min].HIPHpmag = HIPstar[k2].Hpmag;
	  star[k1min].HIPSpType = HIPstar[k2].SpType;
	  star[k1min].HIPBT = HIPstar[k2].BT;
	  star[k1min].HIPVT = HIPstar[k2].VT;
	  star[k1min].HIPBV = HIPstar[k2].BV;
	  star[k1min].HIPVI = HIPstar[k2].VI;
	  HIPmatch++;
	}
      else // ADD NEW ENTRY
	{
	  k1min = NBstar;
	  star[k1min].RAdeg = HIPstar[k2].RAdeg;
	  star[k1min].DECdeg = HIPstar[k2].DECdeg;
	  star[k1min].HIPindex = HIPstar[k2].index;
	  star[k1min].dist = HIPstar[k2].dist;
	  star[k1min].HIPVmag = HIPstar[k2].Vmag;
	  star[k1min].HIPHpmag = HIPstar[k2].Hpmag;
	  star[k1min].HIPSpType = HIPstar[k2].SpType;
	  star[k1min].HIPBT = HIPstar[k2].BT;
	  star[k1min].HIPVT = HIPstar[k2].VT;
	  star[k1min].HIPBV = HIPstar[k2].BV;
	  star[k1min].HIPVI = HIPstar[k2].VI;
	  NBstar++;
	}
      star[k1min].FLAG_HIPS = 1;


      fprintf(fp,"%ld %g %ld %g %g %g %g %g %g %g\n", k2, distmin, HIPstar[k2].index, HIPstar[k2].RAdeg, HIPstar[k2].DECdeg, HIPstar[k2].Vmag, HIPstar[k2].dist, HIPstar[k2].BV, star[k1].Vemag, star[k1].mV);
      
      
    }
  printf("%ld / %ld HIPPARCOS sources matched, %ld others added -> %ld\n", HIPmatch, NBstarHIP, NBstarHIP-HIPmatch, NBstar);
  fclose(fp);

  // ADDING HIPPARCOS DOUBLES
  printf("\nAdding Hipparcos doubles  ...\n");
  HIPDmatch = 0;
  NBadd = 0;
  for(k2=0; k2<NBstarHIPD; k2++) // scan HIPPARCOS MULTIPLES catalog
    {
      printf("\r[%ld/%ld]   ", k2, NBstarHIPD);
      fflush(stdout);
      if(HIPDstar[k2].component == 0) // only add primary
	{
	  RA2000 = HIPDstar[k2].RAdeg/180.0*M_PI;
	  DEC2000 = HIPDstar[k2].DECdeg/180.0*M_PI;
	  
	  distmin = 1000.0; 
	  for(k1=0;k1<NBstar0;k1++)
	    {
	      dist = acos( sin(DEC2000)*sin(star[k1].DECdeg/180.0*M_PI)+cos(DEC2000)*cos(star[k1].DECdeg/180.0*M_PI)*cos(RA2000-star[k1].RAdeg/180.0*M_PI) );
	      dist *= 180.0/M_PI;
	      if(dist<distmin)
		{
		  distmin = dist;
		  k1min = k1;
		}	  
	    }
	  
	  if(distmin<0.03) // UPDATE EXISTING ENTRY
	    {
	      star[k1min].RAdeg = HIPDstar[k2].RAdeg;
	      star[k1min].DECdeg = HIPDstar[k2].DECdeg;
	      star[k1min].HIPindex = HIPDstar[k2].index;
	      star[k1min].dist = HIPDstar[k2].dist; // HIPPARCOS TAKES PRIORITY ON PARALLAX
	      star[k1min].HIPHpmag = HIPDstar[k2].Hpmag;
	      star[k1min].HIPBT = HIPDstar[k2].BT;
	      star[k1min].HIPVT = HIPDstar[k2].VT;
	      HIPDmatch++;
	    }
	  else // ADD NEW ENTRY
	    {
	      k1min = NBstar;
	      star[k1min].RAdeg = HIPDstar[k2].RAdeg;
	      star[k1min].DECdeg = HIPDstar[k2].DECdeg;
	      star[k1min].HIPindex = HIPDstar[k2].index;
	      star[k1min].dist = HIPDstar[k2].dist;
	      star[k1min].HIPHpmag = HIPDstar[k2].Hpmag;
	      star[k1min].HIPBT = HIPDstar[k2].BT;
	      star[k1min].HIPVT = HIPDstar[k2].VT;
	      NBstar++;	      
	      NBadd++;
	    }
	  star[k1min].FLAG_HIPD = 1;
	}
    }
  printf("%ld / %ld HIPPARCOS double sources primary matched, %ld other primaries added -> %ld\n", HIPDmatch, NBstarHIPD, NBadd, NBstar);



  // ADDING RECONS
  printf("\nAdding RECONS sources ...\n");
  NBstar0 = NBstar;
  RECONSmatch = 0;
  NBadd = 0;
  printf("\n");
  for(k2=0; k2<NBstarRECONS; k2++) // scan RECONS catalog
    {
      printf("\r[%ld/%ld]   ", k2, NBstarRECONS);
      fflush(stdout);
     
      RA2000 = RECONSstar[k2].RAdeg/180.0*M_PI;
      DEC2000 = RECONSstar[k2].DECdeg/180.0*M_PI;
      //      printf("  %f %f\n", RA2000, DEC2000);
      // fflush(stdout);
      
      distmin = 1000.0; 
      for(k1=0;k1<NBstar0;k1++)
	{
	  dist = acos( sin(DEC2000)*sin(star[k1].DECdeg/180.0*M_PI)+cos(DEC2000)*cos(star[k1].DECdeg/180.0*M_PI)*cos(RA2000-star[k1].RAdeg/180.0*M_PI) );
	  dist *= 180.0/M_PI;
	  
	  if(star[k1].FLAG_2MASSCNS3 == 1)
	    {
	      m3 = strcmp(star[k1].name3, RECONSstar[k2].name3);    
	      if(star[k1].name3[0]==RECONSstar[k2].name3[0])
		m3 = 0;
	      if((star[k1].name3[0]=='\0')&&(RECONSstar[k2].name3[0]=='A'))
		m3 = 0;
	      if((star[k1].name3[0]=='A')&&(RECONSstar[k2].name3[0]=='\0'))
		m3 = 0;
	      if((star[k1].name3[0]==' ')&&(RECONSstar[k2].name3[0]=='A'))
		m3 = 0;
	      if((star[k1].name3[0]=='A')&&(RECONSstar[k2].name3[0]==' '))
		m3 = 0;	      
	    }
	  else	  
	    m3 = 0;
	  
	  if((dist<distmin)&&(m3==0))
	    {
	      distmin = dist;
	      k1min = k1;
	    }	  
	}
      //   printf("DISTMIN = %f (%ld)\n", distmin, k1min);
      //fflush(stdout);


      if(distmin<0.03) // UPDATE EXISTING ENTRY
	{
	  //  printf("UPDATE ENTRY %ld (%ld)\n", k1, k2);
	  // fflush(stdout);
	  strcpy(star[k1].RECONSname1, RECONSstar[k2].name1);
	  strcpy(star[k1].RECONSname2, RECONSstar[k2].name2);
	  strcpy(star[k1].RECONSname3, RECONSstar[k2].name3);
	  star[k1min].RECONSmV = RECONSstar[k2].mV;
	  star[k1min].RECONSMV = RECONSstar[k2].MV;
	  star[k1min].RECONSdist = RECONSstar[k2].dist;
	  star[k1min].RECONSSpType = RECONSstar[k2].SpType;
	  star[k1min].RECONSmass = RECONSstar[k2].mass;
	  RECONSmatch++;
	}
      else // ADD NEW ENTRY
	{
	  printf("\nADD RECONS ENTRY %ld (%ld) - %s %s %s\n", k1, k2, RECONSstar[k2].name1, RECONSstar[k2].name2, RECONSstar[k2].name3);
	  fflush(stdout);
	  k1min = NBstar;
	  strcpy(star[k1min].RECONSname1, RECONSstar[k2].name1);
	  strcpy(star[k1min].RECONSname2, RECONSstar[k2].name2);
	  strcpy(star[k1min].RECONSname3, RECONSstar[k2].name3);
	  star[k1min].RECONSmV = RECONSstar[k2].mV;
	  star[k1min].RECONSMV = RECONSstar[k2].MV;
	  star[k1min].RECONSdist = RECONSstar[k2].dist;
	  star[k1min].RECONSSpType = RECONSstar[k2].SpType;
	  star[k1min].RECONSmass = RECONSstar[k2].mass;


	  strcpy(star[k1min].name1, RECONSstar[k2].name1);
	  strcpy(star[k1min].name2, RECONSstar[k2].name2);
	  strcpy(star[k1min].name3, RECONSstar[k2].name3);
	  star[k1min].RAdeg = RECONSstar[k2].RAdeg;
	  star[k1min].DECdeg = RECONSstar[k2].DECdeg;

	  NBstar++;	      
	  NBadd++;
	}
      star[k1min].FLAG_RECONS = 1;       
    }
  printf("%ld / %ld RECONS sources primary matched, %ld other added -> %ld\n", RECONSmatch, NBstarRECONS, NBadd, NBstar);
 
  return(0);
}






int Read_MiscData()
{
  FILE *fp;
  long k;
  long j; // entry number
  int n;
  char line[500];
  char word1[200];
  char word2[200];
	char fname[200];

  char name[200];
  int name_entry = 0;
  float Teff;
  int Teff_entry = 0;
  float Blum;
  int Blum_entry = 0;

  int OKline = 0;
  int match = 0;

  printf("READING ADDITIONAL MISC DATA ... \n");
  
  k = 8617;
  printf("name %ld = \"%s\"\n", k, star[k].name);

//  fp = fopen("miscdata.dat","r");
    sprintf(fname, "miscdata.dat"); 
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  
  
  
  j = 0;
  while((fgets(line,SBUFFERSIZE,fp)!=NULL))
    {					      
      if(line[0]!='#')
		{
			n = sscanf(line, "%s %s", word1, word2);
			if(n<1)
				{
					if(OKline==1)
						{
							k = 0;
							while((match==0)&&(k<NBstar))
							{		     
							if(strncmp(star[k].name, name, strlen(name))==0)
								{
									printf("Match found: %ld \"%s\" \"%s\"\n", k, star[k].name, name);
									match = 1;

									if(Teff_entry==1)
									{
										printf("=== ENTERING Teff=%f for %s ===\n", Teff, name);
										star[k].TeffM = Teff;
									}
									if(Blum_entry==1)
									{
										printf("=== ENTERING Blum=%f for %s ===\n", Teff, name);
										star[k].Blum = Blum;
									}

								}		    
							k++;
							}		
						if(match==0)
							printf("%s: No match found\n", name);
						match = 0;
						}
					OKline = 0;
				}
			else
				{
					if(OKline == 0)
						{
							name_entry = 0;
							Teff_entry = 0;		  
							Blum_entry = 0;		  
							j++;
						}
					OKline = 1;
					if(strcmp(word1,"name")==0)
					{
						sprintf(name, "%s", word2);
						//  printf("%ld name = %s\n", j, name);
						name_entry = 1;
					}
					if(strcmp(word1,"Teff")==0)
					{
						Teff = atof(word2);
						//printf("%ld Teff = %f\n", j, Teff);
						Teff_entry = 1;
					}
					if(strcmp(word1, "Blum")==0)
					{
						Blum = atof(word2);
						Blum_entry = 1;
					}
				}
	  

		}
    }
 
  fclose(fp);

  printf("misc data: measurements for %ld sources integrated\n", j);

  // exit(0);

  return(0);
}






















int main()
{
  FILE *fp;
  FILE *fpjs; // for javascript output
  long k;
  float Mbol;
  long tmpl;
  char line[5000];
  float tmpRAdeg;
  float tmpDECdeg, tmpdist, tmpVJ, tmpSpectralIndex;
  long tmpHIPindex;
  float tpmmV, tmpHIPHpmag, tmpJmag, tmpHmag, tmpKmag, tmpMV, tmpBV, tmpHIPVI, tmpVI, tmpTeff, tmpsimulx, tmpBlum, tmpHZSeparcsec, tmpHZContrast, tmpMsun, tmpVK, tmpVH;
  int tmpFLAG_SUPERBLINK, tmpFLAG_2MASSCNS3, tmpFLAG_CNS3, tmpFLAG_HIPS, tmpFLAG_HIPD, tmpFLAG_RECONS;
  char tmpname[200];
  //  char tmpname2[200];
  //  char tmpname3[200];

  float zeropt;
  float zeroptB, zeroptV, zeroptR, zeroptI, zeroptJ, zeroptH, zeroptK;
  float lambda, dlambda, colorcorr;

  long i;
  float zodibckg, ezodibckg;
  float Sphot0, Sphot, Pphot, Zphot, EZphot, SNR;
  float starDiam, starDiamld;

  int WRITEJS = 1;
  long NBJSentries = 0;
  long JScnt = 0;
  float JSdistlimit = 10.0; // distance limit for javascript catalog

  int OKstar = 1;

  // coronagraph leak due to stellar angular size
  float corleakarray[100];
  float corleakarrays[100];
  float v1, v2;
  long j, jmin, jmax;


	double requ_Contrast; // minimum raw contrast to reach required SNR
	double ContrastQ; // factor between source contrast and requ contrast
	double rSNR;

	char string[100];
	double starld;
	double AchievedContrast;
	long corodata_index;

	char fname[200];



	// read spectral band
	fp = fopen("conf_spectralband.txt", "r");
	if(fp!=NULL)
	{
		printf("Reading spectral band\n");
		fscanf(fp, "%s", string);
		fclose(fp);
		spectralBand = string[0];
	}

	printf("SPECTRAL BAND : %c\n", spectralBand);





  read_EEMtable();
	
	read_ContrastData();





  // LEAK PROFILE FOR 0.01 l/D RADIUS STAR
  // FORMAT IS EVERY TENTH OF L/D
//  fp = fopen("corleak.dat","r");
    sprintf(fname, "corleak.dat"); 
	fp = fopen(fname, "r");
	if(fp==NULL)
	{
		printf("ERROR: File \"%s\" missing\n", fname);
		exit(0);
	}
  
  
  for(i=0;i<100;i++)
    {
      fscanf(fp,"%f %f\n", &v1, &v2);
      corleakarray[i] = v2;
    }
  fclose(fp);
  // SMOOTH BY ~0.5 l/D
  for(i=0;i<100;i++)
    {
      jmin = i-5;
      jmax = i+5;
      if(jmin<0)
	jmin = 0;
      if(jmax>99)
	jmax = 99;
      corleakarrays[i] = 0;
      for(j=jmin;j<jmax;j++)
	corleakarrays[i] += corleakarray[j];
      corleakarrays[i] /= (jmax-jmin);
    }

  for(i=0;i<100;i++)
    printf("%f %g\n", 0.1*i, corleakarrays[i]);
  //  exit(0);

  star = (STAR*) malloc(sizeof(STAR)*200000);
  for (k=0;k<200000;k++)
    {
      star[k].FLAG_SUPERBLINK = 0;
      star[k].FLAG_2MASSCNS3 = 0;
      star[k].FLAG_CNS3 = 0;
      star[k].FLAG_HIPS = 0;
      star[k].FLAG_HIPD = 0;
      star[k].FLAG_RECONS = 0;

      strcpy(star[k].name,"\0");
      strcpy(star[k].name1,"\0");
      strcpy(star[k].name2,"\0");
      strcpy(star[k].name3,"\0");
      
      star[k].RAdeg = -100.0;
      star[k].DECdeg = -100.0;
      
      star[k].dist = -100.0; // pc
      star[k].disterr = -100.0; // pc

      star[k].VJ = -100.0; // V-J
      
      star[k].BTmag = -100.0; // Tycho B mag
      star[k].VTmag = -100.0; // Tycho V mag
      star[k].BJmag = -100.0; // USNO B mag
      star[k].RFmag = -100.0; // USNO R mag
      star[k].INmag = -100.0; // USNO I mag
      
      star[k].Vemag = -100.0; // SUPERBLINK effective Vmag

      star[k].Jmag = -100.0; // 2MASS
      star[k].Hmag = -100.0; // 2MASS
      star[k].Kmag = -100.0; // 2MASS
      star[k].SpectralIndex = -100.0; // from SUPERBLINK .. -1 if unknown

      star[k].VI = -100.0;
      star[k].VK = -100.0;
      star[k].JH = -100.0;
      star[k].HK = -100.0;
      
      star[k].TeffM = -100.0;

      strcpy(star[k].CNS3name1,"\0");
      strcpy(star[k].CNS3name2,"\0");
      strcpy(star[k].CNS3name3,"\0");
      strcpy(star[k].SUPERBLINKCNS3name1,"\0");
      strcpy(star[k].SUPERBLINKCNS3name2,"\0");
      strcpy(star[k].SUPERBLINKCNS3name3,"\0");
      star[k].SpType = -100.0; // from CNS3 ... -1 if unknown
      star[k].mV = -100.0; // from CNS3
      star[k].BV = -100.0; // from CNS3
      star[k].UB = -100.0; // from CNS3
      star[k].RI = -100.0; // from CNS3
      star[k].CNS3dist = -100.0; // from CNS3
      star[k].MV = -100.0; // from CNS3
      
      star[k].HIPindex = -100;
      star[k].HIPVmag = -100.0; // HIPPARCOS Vmag
      star[k].HIPHpmag = -100.0;
      star[k].HIPSpType = -100.0;
      star[k].HIPBT = -100.0;
      star[k].HIPVT = -100.0;
      star[k].HIPBV = -100.0;
      star[k].HIPVI = -100.0;

      strcpy(star[k].RECONSname1, "\0");
      strcpy(star[k].RECONSname2, "\0");
      strcpy(star[k].RECONSname3, "\0");
      star[k].RECONSmV = -100.0;
      star[k].RECONSMV = -100.0;
      star[k].RECONSdist = -100.0;
      star[k].RECONSSpType = -100.0;
      star[k].RECONSmass = -100.0;

      star[k].simulx = -100.0;
      star[k].Teff = -100.0;
      star[k].Lmag = -100.0;
      star[k].Blum = -100.0; // Bolometric Luminosity
      star[k].HZSeparcsec = -100.0;
      star[k].HZContrast = -100.0;
      star[k].Msun = -100.0;
    }


  // READ STARLIST
  if((fp = fopen("StarList.txt","r"))!=NULL)
    {
      k = 0;
      //      while(fscanf(fp,"%ld %f %f %f %f %f %ld %f %f\n", &tmpl, &tmpRAdeg, &tmpDECdeg, &tmpdist, &tmpVJ, &tmpSpectralIndex, &tmpHIPindex, &tpmmV, &tmpHIPHpmag)>5)
      
      while(fgets ( line, sizeof line, fp ) != NULL)
	{
	  sscanf(line,"%ld %f %f %f %f %f %ld %f %f %f %f %f %f %f %f %f %f %f %f %f %g %g %d %d %d %d %d %d %s %f %f\n", &tmpl, &tmpRAdeg, &tmpDECdeg, &tmpdist, &tmpVJ, &tmpSpectralIndex, &tmpHIPindex, &tpmmV, &tmpHIPHpmag, &tmpJmag, &tmpHmag, &tmpKmag, &tmpMV, &tmpBV, &tmpHIPVI, &tmpVI, &tmpTeff, &tmpsimulx, &tmpBlum, &tmpHZSeparcsec, &tmpHZContrast, &tmpMsun, &tmpFLAG_SUPERBLINK, &tmpFLAG_2MASSCNS3, &tmpFLAG_CNS3, &tmpFLAG_HIPS, &tmpFLAG_HIPD, &tmpFLAG_RECONS, tmpname, &tmpVK, &tmpVH);
	//	{
	  
	  star[k].RAdeg = tmpRAdeg;
	  star[k].DECdeg = tmpDECdeg;
	  star[k].dist = tmpdist; 
	  star[k].VJ = tmpVJ;
	  star[k].SpectralIndex = tmpSpectralIndex;
	  star[k].HIPindex = tmpHIPindex;
	  star[k].mV = tpmmV;
	  star[k].HIPHpmag = tmpHIPHpmag;
	  star[k].Jmag = tmpJmag;
	  star[k].Hmag = tmpHmag;
	  star[k].Kmag = tmpKmag;
	  star[k].MV = tmpMV;
	  star[k].BV = tmpBV;
	  star[k].HIPVI = tmpHIPVI;
	  star[k].VI = tmpVI;
	  star[k].Teff = tmpTeff;
	  star[k].simulx = tmpsimulx;
	  star[k].Blum = tmpBlum;
	  star[k].HZSeparcsec = tmpHZSeparcsec;
	  star[k].HZContrast = tmpHZContrast;
	  star[k].Msun = tmpMsun;
	  star[k].FLAG_SUPERBLINK = tmpFLAG_SUPERBLINK;
	  star[k].FLAG_2MASSCNS3 = tmpFLAG_2MASSCNS3;
	  star[k].FLAG_CNS3 = tmpFLAG_CNS3;
	  star[k].FLAG_HIPS = tmpFLAG_HIPS;
	  star[k].FLAG_HIPD = tmpFLAG_HIPD;
	  star[k].FLAG_RECONS = tmpFLAG_RECONS;
	  strcpy(star[k].name, tmpname);
	  star[k].VK = tmpVK;
	  star[k].VH = tmpVH;
	  
	  //	  printf("--------- %ld %ld\n", k, tmpl);
	  k++;
	}
      NBstar = k;
      printf("%ld entries read\n", NBstar);
      fclose(fp);
    }
  else
    {
      RECONSstar = (RECONSSTAR*) malloc(sizeof(RECONSSTAR)*300);
      read_RECONS();
      
      read_SUPERBLINK();
      
      CNS3star = (CNS3STAR*) malloc(sizeof(CNS3STAR)*5000);
      read_CNS3();
      
      CNS3_2MASS_star = (CNS3_2MASS_STAR*) malloc(sizeof(CNS3_2MASS_STAR)*5000);
      read_CNS3_2MASS();
      
      HIPstar = (HIPSTAR*) malloc(sizeof(HIPSTAR)*120000);
      read_HIP();
      
      HIPDstar = (HIPDSTAR*) malloc(sizeof(HIPDSTAR)*30000);
      read_HIPD();
      
      // MERGE CATALOGS TO CREATE LIST OF ENTRIES
      MatchCatalogs();
  


      // CREATE SINGLE NAME
      
     for(k=0;k<NBstar;k++)
	{
	  sprintf(star[k].name, "%s%s%s", star[k].name1, star[k].name2, star[k].name3);
	  for(i=0;i<strlen(star[k].name);i++)
	    if(star[k].name[i]==' ')
	      star[k].name[i] = '_';
	  i = strlen(star[k].name)-1;
	  while((star[k].name[i] == '_')&&(i>0))
	    {
	      star[k].name[i] = ' ';
	      i--;
	    }
	}

      // read additional data
      Read_MiscData();



      // COMBINE INFO FOR EACH ENTRY
      fp = fopen("StarList.txt","w");
      for(k=0;k<NBstar;k++)
	{
	  OKstar = 1; 


	  /*	  sprintf(star[k].name, "%s%s%s", star[k].name1, star[k].name2, star[k].name3);
	  for(i=0;i<strlen(star[k].name);i++)
	    if(star[k].name[i]==' ')
	      star[k].name[i] = '_';
	  i = strlen(star[k].name)-1;
	  while((star[k].name[i] == '_')&&(i>0))
	    {
	      star[k].name[i] = ' ';
	      i--;
	    }
	  */


	  // distance 
	  // match routine did use SUPERBLINK (trig then phot), then CSN3, then HIP
	  if(star[k].FLAG_RECONS==1) // if in RECONS, use RECONS distance
	    star[k].dist = star[k].RECONSdist;
	  
	  if(star[k].dist<1.0)
	    {
	      //	      printf("COULD NOT ESTIMATE distance FOR STAR %ld -> DISCARDING ENTRY\n", k);
	      OKstar = 0;
	    }


	  if((star[k].dist<1.4)&&(OKstar==1))
	    {
	      VERBOSE = 1;
	      printf("-------- STAR %ld ---- dist = %f ----\n", k, star[k].dist);
	    }
	  else
	    VERBOSE = 0;
	 


	  // Spectral type
	  // match routine did use SUPERBLINK 
	  //if((star[k].dist<DISTMAX)&&(star[k].dist>0.1))
	  // {
	      if(star[k].SpectralIndex<0)
		star[k].SpectralIndex = star[k].SpType; // from CNS3
	      if((star[k].SpectralIndex<0)&&((star[k].FLAG_HIPS==1)||(star[k].FLAG_HIPD==1)))
		star[k].SpectralIndex = star[k].HIPSpType; // from HIPPARCOS
	      if(star[k].RECONSSpType>-10.0) // use RECONS spectral type if available
		star[k].SpectralIndex = star[k].RECONSSpType;
	      
	      // mV & MV
	      // match routine uses CNS3
	      if(star[k].mV<-10.0)
		{	      
		  if(star[k].Vemag > -10.0) // SUPERBLINK
		    star[k].mV = star[k].Vemag;
		  if(star[k].VTmag > -10.0) // SUPERBLINK
		    star[k].mV = star[k].VTmag;
		  if(star[k].HIPVmag > -10.0) // HIP
		    star[k].mV = star[k].HIPVmag;	      
		}
	      if(star[k].RECONSmV > -10.0) // use RECONS if available
		star[k].mV = star[k].RECONSmV;
	      
	      if((star[k].MV<-10.0)&&(star[k].mV>-10.0))
		star[k].MV = star[k].mV - 5.0*(log10(star[k].dist)-1.0);
	      
	      // use HIP mag as backup
	      if((star[k].MV<-10.0)&&(star[k].HIPHpmag>-10.0))
		star[k].MV = star[k].HIPHpmag - 5.0*(log10(star[k].dist)-1.0);
	  
	      if(OKstar==1)
		if(star[k].MV<-10.0)
		  {
		    //		    printf("COULD NOT ESTIMATE MV FOR STAR %ld -> DISCARDING ENTRY\n", k);
		    OKstar = 0;
		  }

	      // B-V
	      if(star[k].HIPBV > -10.0)
		star[k].BV = star[k].HIPBV;
	      
	      // V-I	  
	      if(star[k].HIPVI > -10.0)
		star[k].VI = star[k].HIPVI;
	      
	      // V-J
	      if(star[k].VJ<-10.0) // compute V-J color if missing
		{
		  if((star[k].Jmag>-10.0)&&(star[k].mV>-10.0))
		    star[k].VJ = star[k].mV-star[k].Jmag;
		}

	      // V-K
	      if(star[k].VK<-10.0) // compute V-K color if missing
		{
		  if((star[k].Kmag>-10.0)&&(star[k].mV>-10.0))
		    star[k].VK = star[k].mV-star[k].Kmag;
		}

	      // J-H
	      if(star[k].JH<-10.0) // compute J-H color if missing
		{
		  if((star[k].Hmag>-10.0)&&(star[k].Jmag>-10.0))
		    star[k].JH = star[k].Jmag - star[k].Hmag;		
		}

	      // H-K
	      if(star[k].HK<-10.0) // compute H-K color if missing
		{
		  if((star[k].Hmag>-10.0)&&(star[k].Kmag>-10.0))
		    star[k].HK = star[k].Hmag - star[k].Kmag;		
		}



	      // fill in the missing data by comparing photometry with MS photometry table
	      if(star[k].SpectralIndex>-1.0)
		{
		  valx = star[k].SpectralIndex;
		  errx = 0.5;
		}
	      else
		{
		  valx = 5.0;
		  errx = 1.0e10;
		}
	  
	      if(star[k].UB>-10.0)
		{
		  valUB = star[k].UB;
		  errUB = 0.1;
		}
	      else
		{
		  valUB = 0.0;
		  errUB = 1.0e10;
		}

	      if(star[k].BV>-10.0)
		{
		  valBV = star[k].BV;
		  errBV = 0.1;
		}
	      else
		{
		  valBV = 0.0;
		  errBV = 1.0e10;
		}
	  
	      if(star[k].VI>-10.0)
		{
		  valVI = star[k].VI;
		  errVI = 0.1;
		}
	      else
		{
		  valVI = 0.0;
		  errVI = 1.0e10;
		}
	  
	      if(star[k].VJ>-10.0)
		{
		  valVJ = star[k].VJ;
		  errVJ = 0.1;
		}
	      else
		{
		  valVJ = 0.0;
		  errVJ = 1.0e10;
		}
	  
	      if(star[k].VK>-10.0)
		{
		  valVK = star[k].VK;
		  errVK = 0.1;
		}
	      else
		{
		  valVK = 0.0;
		  errVK = 1.0e10;
		}
	  
	  
	      if(star[k].JH>-10.0)
		{
		  valJH = star[k].JH;
		  errJH = 0.1;
		}
	      else
		{
		  valJH = 0.0;
		  errJH = 1.0e10;
		}
	  
	  
	      if(star[k].HK>-10.0)
		{
		  valHK = star[k].HK;
		  errHK = 0.1;
		}
	      else
		{
		  valHK = 0.0;
		  errHK = 1.0e10;
		}
	  
	      if(star[k].TeffM>10.0)
		{
		  valTeff = star[k].TeffM;
		  errTeff = 0.01;
		}
	      else
		{
		  valTeff = 5000.0;
		  errTeff = 100000000.0;		
		}


	      star[k].Teff = Teff_solve();
	      
	      
	      if(star[k].VI<-50.0)
			star[k].VI = simulVI;
	      
	      star[k].simulx = simulx;

	      // compute bolometric luminosity
	      Mbol = star[k].MV + simulBCv - 0.076;
	      //	  BCSun = -0.076;
	      if(star[k].Blum < 0.0)
			star[k].Blum = pow(2.51188643,-(Mbol-4.83));

		if(simulMode==1) // replace all color measurements by EEM values
		{
			// input is MV, dist

			// colors
			star[k].UB = simulUB;
			star[k].BV = simulBV;
			star[k].VI = simulVI;
			star[k].VJ = simulVJ;
			star[k].VH = simulVH;
			star[k].VK = simulVK;
			star[k].JH = simulJH;
			star[k].HK = simulHK;
			
			// magnitudes
			star[k].mV = star[k].MV + 5.0*(log10(star[k].dist)-1.0);
			star[k].Jmag = star[k].mV - star[k].VJ;
			star[k].Hmag = star[k].mV - star[k].VH;
			star[k].Kmag = star[k].mV - star[k].VK;			
		}

	      star[k].HZSeparcsec = sqrt(star[k].Blum)/star[k].dist;
	      star[k].HZContrast = 1.74e-10/star[k].Blum; 
	      //	  Earth-like planet, albedo = 0.3, phase = pi/2

	      if(star[k].RECONSmass>0.0)
		star[k].Msun = star[k].RECONSmass;
	      else if (star[k].SpectralIndex>1.0)
		star[k].Msun = simulMsun;
	      else if (fabs(star[k].SpectralIndex)<0.1)
		star[k].Msun = 0.6; // default value for white dwarfs



	      if(OKstar==1)
		{
		  // WRITE RESULTS
		  
		  // col 1: index
		  // col 2: RA deg
		  // col 3: DEC deg
		  // col 4: dist
		  // col 5: V-J
		  // col 6: spectral index
		  // col 7: HIP index
		  // col 8: V mag
		  // col 9: HIP mag
		  // col 10: J mag
		  // col 11: H mag
		  // col 12: K mag
		  // col 13: MV
		  // col 14: B-V
		  // col 15: HIP V-I
		  // col 16: V-I
		  // col 17: Teff estimate
		  // col 18: spectral type estimate
		  // col 19: Blum
		  // col 20: HZ sep arcsec
		  // col 21: HZ contrast (for Earth)
		  // col 22: stellar mass      
		  // col 23: FLAG (SUPERBLINK)
		  // col 24: FLAG (2MASSCNS3)
		  // col 25: FLAG (CNS3)
		  // col 26: FLAG (HIPS)
		  // col 27: FLAG (HIPD)
		  // col 28: FLAG (RECONS)
		  // col 29: name
		  // col 30: V-K
		  // col 31: V-H

		  
		  
		  
		  fprintf(fp,"%ld %f %f %f %f %f %ld %f %f %f %f %f %f %f %f %f %f %f %f %f %g %g %d %d %d %d %d %d %s %f %f\n", k, star[k].RAdeg, star[k].DECdeg, star[k].dist, star[k].VJ, star[k].SpectralIndex, star[k].HIPindex, star[k].mV, star[k].HIPHpmag, star[k].Jmag, star[k].Hmag, star[k].Kmag, star[k].MV, star[k].BV, star[k].HIPVI, star[k].VI, star[k].Teff, star[k].simulx, star[k].Blum, star[k].HZSeparcsec, star[k].HZContrast, star[k].Msun, star[k].FLAG_SUPERBLINK, star[k].FLAG_2MASSCNS3, star[k].FLAG_CNS3, star[k].FLAG_HIPS, star[k].FLAG_HIPD, star[k].FLAG_RECONS, star[k].name, star[k].VK, star[k].VH);
		}
	}
      fclose(fp);
    }

  if(1==1)
    {

      zeroptB = 1.3846e11;
      zeroptV = 9.9690e10;
      zeroptR = 7.2384e10;
      zeroptI = 4.5825e10;
      zeroptJ = 1.9422e10;
      zeroptH = 9.4440e9;
      zeroptK = 4.3829e9;
      

      fp = fopen("ResultSNR.txt","w");
      
      if(WRITEJS==1)
	{
	  fpjs = fopen("catalog.js", "w");
	  for(k=0;k<NBstar;k++)
	    if((star[k].dist < JSdistlimit)&&(star[k].dist>1.0))
	      JScnt++;
	  NBJSentries = JScnt+1;
	  printf("%ld entries in JS catalog\n", NBJSentries);
	}

      fprintf(fp, "# Telescope diameter [m]       : %.2f \n", telDiam);
      fprintf(fp, "# Optical system Efficiency    : %.4f\n", Efficiency);
      fprintf(fp, "# Spectral resolution          : %.2f\n", spectralR);
      fprintf(fp, "# spectralBand                 : %c\n", spectralBand);
      
      fprintf(fp, "# Coronagraph IWA [l/D]        : %f\n", IWAld);
      fprintf(fp, "# Coronagraph image FWHM [l/D] : %f\n", FWHMld); 
      fprintf(fp, "# Coronagraph throughput       : %f\n", Throughput); 
      fprintf(fp, "# Coronagraph Airy throughput  : %f\n", AiryThroughput);
      fprintf(fp, "# Coronagraph raw contrast 1   : %f\n", RAWcontrast1);
      fprintf(fp, "# high contrast angle limit    : %f\n", IWAld1);
      fprintf(fp, "# Coronagraph raw contrast 2   : %f\n", RAWcontrast2);
      fprintf(fp, "# Coronagraph contrast factor  : %f\n", contrastFactor);




      fprintf(fp, "\n");

      if(WRITEJS==1)
	{
	  JScnt = 0;
      fprintf(fpjs, "var catalogNBtargets = %ld;\n", NBJSentries);
      fprintf(fpjs, "var catalogTelDiam = %f;\n", telDiam);

      fprintf(fpjs, "var catalogRAdeg = new Array();\n");
      fprintf(fpjs, "var catalogDECdeg = new Array();\n");
      fprintf(fpjs, "var catalogdist = new Array();\n");
      fprintf(fpjs, "var catalogVJ = new Array();\n");
      fprintf(fpjs, "var catalogSprectalIndex = new Array();\n");
      fprintf(fpjs, "var catalogHIPindex = new Array();\n");
      fprintf(fpjs, "var catalogVmag = new Array();\n");
      fprintf(fpjs, "var catalogHIPmag = new Array();\n");
      fprintf(fpjs, "var catalogJmag = new Array();\n");
      fprintf(fpjs, "var catalogHmag = new Array();\n");
      fprintf(fpjs, "var catalogKmag = new Array();\n");
      fprintf(fpjs, "var catalogMV = new Array();\n");
      fprintf(fpjs, "var catalogBV = new Array();\n");
      fprintf(fpjs, "var catalogHIP_VI = new Array();\n");
      fprintf(fpjs, "var catalogVI = new Array();\n");
      fprintf(fpjs, "var catalogTeff = new Array();\n");
      fprintf(fpjs, "var catalogSpectralType = new Array();\n");
      fprintf(fpjs, "var catalogBlum = new Array();\n");
      fprintf(fpjs, "var catalogHZarcsec = new Array();\n");
      fprintf(fpjs, "var catalogHZlogcontrast = new Array();\n");
      fprintf(fpjs, "var catalogStellarMass = new Array();\n");
      fprintf(fpjs, "var catalogFLAG_SUPERBLINK = new Array();\n");
      fprintf(fpjs, "var catalogFLAG_2MASSCNS3 = new Array();\n");
      fprintf(fpjs, "var catalogFLAG_CNS3 = new Array();\n");
      fprintf(fpjs, "var catalogFLAG_HIPS = new Array();\n");
      fprintf(fpjs, "var catalogFLAG_HIPD = new Array();\n");
      fprintf(fpjs, "var catalogFLAG_RECONS = new Array();\n");
      fprintf(fpjs, "var catalogPlanetPhotonRate = new Array();\n");
      fprintf(fpjs, "var catalogStarLeakPhotonRate = new Array();\n");
      fprintf(fpjs, "var catalogZodiPhotonRate = new Array();\n");
      fprintf(fpjs, "var catalogExoZodiPhotonRate = new Array();\n");
      fprintf(fpjs, "var catalogSNR1hr = new Array();\n");
      fprintf(fpjs, "var catalogName = new Array();\n");
	}


	  fprintf(fp, "\n");
	  fprintf(fp, "\n");
	  
	fprintf(fp, "# col 1: index\n");
	fprintf(fp, "# col 2: RA deg\n");
	fprintf(fp, "# col 3: DEC deg \n");
	fprintf(fp, "# col 4: dist\n");
    fprintf(fp, "# col 5: V-J\n");
    fprintf(fp, "# col 6: spectral index\n");
    fprintf(fp, "# col 7: HIP index\n");
    fprintf(fp, "# col 8: V mag\n");
    fprintf(fp, "# col 9: HIP mag\n");
    fprintf(fp, "# col 10: J mag\n");
    fprintf(fp, "# col 11: H mag\n");
    fprintf(fp, "# col 12: K mag\n");
    fprintf(fp, "# col 13: MV\n");
    fprintf(fp, "# col 14: B-V\n");
    fprintf(fp, "# col 15: HIP V-I\n");
    fprintf(fp, "# col 16: V-I\n");
    fprintf(fp, "# col 17: Teff estimate\n");
    fprintf(fp, "# col 18: spectral type estimate\n");
    fprintf(fp, "# col 19: Blum\n");
    fprintf(fp, "# col 20: HZ sep arcsec\n");
    fprintf(fp, "# col 21: HZ contrast (for Earth)\n");
    fprintf(fp, "# col 22: stellar mass      \n");
    fprintf(fp, "# col 23: FLAG (SUPERBLINK)\n");
    fprintf(fp, "# col 24: FLAG (2MASSCNS3)\n");
    fprintf(fp, "# col 25: FLAG (CNS3)\n");
    fprintf(fp, "# col 26: FLAG (HIPS)\n");
    fprintf(fp, "# col 27: FLAG (HIPD)\n");
    fprintf(fp, "# col 28: FLAG (RECONS)\n");
    fprintf(fp, "# col 29: Star diameter [Sun]\n");
    fprintf(fp, "# col 30: Star photon rate\n");
    fprintf(fp, "# col 31: Planet photon rate (per sec)\n");
    fprintf(fp, "# col 32: Star leak photon rate (per sec)\n");
    fprintf(fp, "# col 33: Zodiacal light photon rate (per sec)\n");
    fprintf(fp, "# col 34: ExoZodi light photon rate (per sec)\n");
    fprintf(fp, "# col 35: SNR (for 1hr exposure)\n");
    fprintf(fp, "# col 36: name\n");
    fprintf(fp, "# col 37: requ RAW contrast\n");
    fprintf(fp, "# col 38: achieved SNR assuming background calib limit\n");
	fprintf(fp, "# col 39: star rad l/D\n");
	fprintf(fp, "# col 40: contrast achieved\n");
	fprintf(fp, "# col 41: HZ sep l/D\n");


	  fprintf(fp, "\n");
	  fprintf(fp, "\n");



	for(k=0;k<NBstar;k++)
	{
		switch ( spectralBand ) {

			case 'B':
			lambda = lambdaB;
			zeropt = zeroptB;
			zodibckg = zodibckgB;
			ezodibckg = zodibckg + (star[k].BV - BVsun);
			colorcorr = -star[k].BV;
			break;

			case 'V':
			lambda = lambdaV;
			zeropt = zeroptV;
			zodibckg = zodibckgV;
			ezodibckg = zodibckg;
			colorcorr = 0.0;
			break;

/*			case 'R':
			lambda = lambdaR;
			zeropt = zeroptR;
			zodibckg = zodibckgR;
			ezodibckg = zodibckgK - (star[k].VK - VKsun);
			colorcorr = star[k].VR;
			break;
	*/		
			case 'I':
			lambda = lambdaI;
			zeropt = zeroptI;
			zodibckg = zodibckgI;
			ezodibckg = zodibckg - (star[k].VI - VIsun);
			colorcorr = star[k].VI;
			break;
			
			case 'J':
			lambda = lambdaJ;
			zeropt = zeroptJ;
			zodibckg = zodibckgJ;
			ezodibckg = zodibckg - (star[k].VJ - VJsun);
			colorcorr = star[k].VJ;
			break;

			case 'H':
			lambda = lambdaH;
			zeropt = zeroptH;
			zodibckg = zodibckgH;
			ezodibckg = zodibckg - ( star[k].VH - VHsun);
			colorcorr = star[k].VH;
			break;

			case 'K':
			lambda = lambdaK;
			zeropt = zeroptK;
			zodibckg = zodibckgK;
			ezodibckg = zodibckg - (star[k].VK - VKsun);
			colorcorr = star[k].VK;
			break;
			
			
			default:
			printf("ERROR: spectral band %c not supported\n", spectralBand);
			exit(0);
			break;
      }


	



	  dlambda = lambda/spectralR; // bandpass

	  // stellar diam
	  // bol lum = diam^2 T^4 -> diam = sqrt(Bol)/T^2
	  starDiam = sqrt(star[k].Blum)/pow(star[k].Teff/5778.0,2.0);
	  starDiamld = starDiam/star[k].dist*0.0093 / (lambda/telDiam/M_PI*180.0*3600.0);

	  if(star[k].HZSeparcsec/3600.0/180.0*M_PI < IWAld1*lambda/telDiam)
	    RAWcontrast = RAWcontrast1;
	  if(star[k].HZSeparcsec/3600.0/180.0*M_PI > IWAld1*lambda/telDiam)
	    RAWcontrast = RAWcontrast2;

	  

	  // Total star photon into coronagraph
	  if(star[k].mV>-10.0)
	    Sphot0 = pow(2.51188643152,-(star[k].mV-colorcorr))*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6);
	  else
	    Sphot0 = pow(2.51188643152,-((star[k].MV-colorcorr) + 5.0*(log10(star[k].dist)-1.0)))*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6);
	  Sphot0 *= Efficiency;

	  // planet photon
	  if(star[k].mV>-10.0)
	    Pphot = pow(2.51188643152,-(star[k].mV-colorcorr))*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6)*star[k].HZContrast;
	  else
	    Pphot = pow(2.51188643152,-((star[k].MV-colorcorr) + 5.0*(log10(star[k].dist)-1.0)))*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6)*star[k].HZContrast;
	  Pphot *= Efficiency*AiryThroughput;
	  
	  
	  // Star photon after coronagraph is added to coronagraph angular leak
	  i = (long) (10.0*star[k].HZSeparcsec/3600.0/180.0*M_PI/(lambda/telDiam)+0.5);
	  if(i>99)
	    i = 99;
	  RAWcontrast += corleakarrays[i]*pow(starDiamld/0.02,2.0);
	  
	  if((star[k].dist>0.1)&&(star[k].dist<5.0))
	    printf("dist = %f   STAR SIZE = %f Sun = %f l/D -> RAWcontrast at %f = %g x %g = %g\n", star[k].dist, starDiam, starDiamld, star[k].HZSeparcsec/3600.0/180.0*M_PI/(lambda/telDiam), corleakarrays[i], pow(starDiamld/0.02,2.0), RAWcontrast);

	  if(star[k].mV>-10.0)
	    Sphot = pow(2.51188643152,-(star[k].mV-colorcorr))*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6)*RAWcontrast;
	  else
	    Sphot = pow(2.51188643152,-((star[k].MV-colorcorr) + 5.0*(log10(star[k].dist)-1.0)))*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6)*RAWcontrast;
	  Sphot *= Efficiency*AiryThroughput;
	  
	  // Zodi photon
	  Zphot = pow(2.51188643152,-zodibckg)*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6)*pow(((lambda/telDiam)/M_PI*180.0*3600.0)*FWHMld,2.0);
	  Zphot *= Efficiency*Throughput;
	  
	  // exozodi photon
		EZphot = pow(2.51188643152,-ezodibckg)*zeropt*(telDiam*telDiam/4.0)*M_PI*(dlambda*1.0e6)*pow(((lambda/telDiam)/M_PI*180.0*3600.0)*FWHMld,2.0);
	  EZphot *= 3.0;
	  
	  
	  
	  // SNR 
	  SNR = Pphot/sqrt(Pphot+Sphot+Zphot+EZphot);
	  if(star[k].HZSeparcsec/3600.0/180.0*M_PI<IWAld*lambda/telDiam)
	    SNR = 0.0;
	  
	  if(star[k].HZContrast<RAWcontrast/contrastFactor)
	    SNR = 0.0;
	  
	  
	  
		starld = 0.5*starDiam/star[k].dist*10.0; // mas 
		starld *= M_PI/(1000.0*3600.0*180.0);
		starld /= (lambda/telDiam);
		
		corodata_index = (long) (1.0*star[k].HZSeparcsec*M_PI/(3600.0*180.0) / (lambda/telDiam) / ContrastDATA_step + 0.5);
		if(corodata_index > ContrastDATA_NBpt-2)
			corodata_index = ContrastDATA_NBpt-1;
			
		AchievedContrast = ContrastCOH[corodata_index] + ContrastINCOH[corodata_index] * pow(starld/SSIZE_INCOH, 2.0);
		
//float *ContrastCOH;
//float *ContrastINCOH;
//float SSIZE_INCOH = 1.0;

	  // ContrastAchieved = ;
		
	  
	  // WRITE RESULTS
	  
	  



      for(i=0;i<strlen(star[k].name);i++)
	if(star[k].name[i]==' ')
	  star[k].name[i] = '_';

      i = strlen(star[k].name)-1;
       while((star[k].name[i] == '_')&&(i>0))
	{
	  star[k].name[i] = ' ';
	  i--;
	}
	  
	requ_Contrast = 1.0/Sphot0 * ( Pphot*Pphot*requ_etime/(requ_SNR*requ_SNR) - (Pphot + Zphot + EZphot) );
	
	ContrastQ = requ_Contrast / ( Pphot/Sphot0 );
	if(ContrastQ>ContrastQlim)
		requ_Contrast = Pphot/Sphot0 * ContrastQlim;
	
	rSNR = ( Pphot*requ_etime ) / sqrt( (Pphot + requ_Contrast*Sphot0 + Zphot + EZphot)*requ_etime );	  
	
	if(star[k].HZSeparcsec/3600.0/180.0*M_PI<IWAld*lambda/telDiam)
	    {
			rSNR = 0.0;
			requ_Contrast = -1.0;
		}


       fprintf(fp,"%8ld %12.8f %+12.8f %7.3f %+7.3f %3.1f %7ld %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f     %8.2f %5.2f %12.6f %8.6f %16g %6.3f %d %d %d %d %d %d  %7.3f    %12g %12g %12g %12g %12g %12g %16s %12g %12g  %6.4f %12g %10.6f\n", k, star[k].RAdeg, star[k].DECdeg, star[k].dist, star[k].VJ, star[k].SpectralIndex, star[k].HIPindex, star[k].mV, star[k].HIPHpmag, star[k].Jmag, star[k].Hmag, star[k].Kmag, star[k].MV, star[k].BV, star[k].HIPVI, star[k].VI, star[k].Teff, star[k].simulx, star[k].Blum, star[k].HZSeparcsec, star[k].HZContrast, star[k].Msun, star[k].FLAG_SUPERBLINK, star[k].FLAG_2MASSCNS3, star[k].FLAG_CNS3, star[k].FLAG_HIPS, star[k].FLAG_HIPD, star[k].FLAG_RECONS, starDiam, Sphot0, Pphot, Sphot, Zphot, EZphot, SNR*60.0, star[k].name, requ_Contrast, rSNR, starld, AchievedContrast, 1.0*star[k].HZSeparcsec*M_PI/(3600.0*180.0) / (lambda/telDiam));
       // %g %g %g %g %g %s
       // Pphot, Sphot, Zphot, EZphot, SNR*60.0, star[k].name);

       
 	 
	      
      if(WRITEJS==1)
	if((star[k].dist < JSdistlimit)&&(star[k].dist>1.0))
	 {
	   fprintf(fpjs, "\n");
	   fprintf(fpjs, "catalogRAdeg[%ld] = %f;\n", JScnt, star[k].RAdeg);
	   fprintf(fpjs, "catalogDECdeg[%ld] = %f;\n", JScnt, star[k].DECdeg);
	   fprintf(fpjs, "catalogdist[%ld] = %f;\n", JScnt, star[k].dist);
	   fprintf(fpjs, "catalogVJ[%ld] = %f;\n", JScnt, star[k].VJ);
	   fprintf(fpjs, "catalogSprectalIndex[%ld] = %f;\n", JScnt, star[k].SpectralIndex);
	   fprintf(fpjs, "catalogHIPindex[%ld] = %ld;\n", JScnt, star[k].HIPindex);
	   fprintf(fpjs, "catalogVmag[%ld] = %f;\n", JScnt, star[k].mV);
	   fprintf(fpjs, "catalogHIPmag[%ld] = %f;\n", JScnt, star[k].HIPHpmag);
	   fprintf(fpjs, "catalogJmag[%ld] = %f;\n", JScnt, star[k].Jmag);
	   fprintf(fpjs, "catalogHmag[%ld] = %f;\n", JScnt, star[k].Hmag);
	   fprintf(fpjs, "catalogKmag[%ld] = %f;\n", JScnt, star[k].Kmag);
	   fprintf(fpjs, "catalogMV[%ld] = %f;\n", JScnt, star[k].MV);
	   fprintf(fpjs, "catalogBV[%ld] = %f;\n", JScnt, star[k].BV);
	   fprintf(fpjs, "catalogHIP_VI[%ld] = %f;\n", JScnt, star[k].HIPVI);
	   fprintf(fpjs, "catalogVI[%ld] = %f;\n", JScnt, star[k].VI);
	   fprintf(fpjs, "catalogTeff[%ld] = %f;\n", JScnt, star[k].Teff);
	   fprintf(fpjs, "catalogSpectralType[%ld] = %f;\n", JScnt, star[k].simulx);
	   fprintf(fpjs, "catalogBlum[%ld] = %f;\n", JScnt, star[k].Blum);
	   fprintf(fpjs, "catalogHZarcsec[%ld] = %f;\n", JScnt, star[k].HZSeparcsec);
	   fprintf(fpjs, "catalogHZlogcontrast[%ld] = %f;\n", JScnt, log10(star[k].HZContrast));
	   fprintf(fpjs, "catalogStellarMass[%ld] = %f;\n", JScnt, star[k].Msun);
	   
	   if(star[k].FLAG_SUPERBLINK==1)
	     fprintf(fpjs, "catalogFLAG_SUPERBLINK[%ld] = true;\n", JScnt);
	   else
	     fprintf(fpjs, "catalogFLAG_SUPERBLINK[%ld] = false;\n", JScnt);
	   
	   if(star[k].FLAG_2MASSCNS3==1)
	     fprintf(fpjs, "catalogFLAG_2MASSCNS3[%ld] = true;\n", JScnt);
	   else
	     fprintf(fpjs, "catalogFLAG_2MASSCNS3[%ld] = false;\n", JScnt);
	   
	   if(star[k].FLAG_CNS3==1)
	     fprintf(fpjs, "catalogFLAG_CNS3[%ld] = true;\n", JScnt);
	   else
	     fprintf(fpjs, "catalogFLAG_CNS3[%ld] = false;\n", JScnt);

	   if(star[k].FLAG_HIPS==1)
	     fprintf(fpjs, "catalogFLAG_HIPS[%ld] = true;\n", JScnt);
	   else
	     fprintf(fpjs, "catalogFLAG_HIPS[%ld] = false;\n", JScnt);

	   if(star[k].FLAG_HIPD==1)
	     fprintf(fpjs, "catalogFLAG_HIPD[%ld] = true;\n", JScnt);
	   else
	     fprintf(fpjs, "catalogFLAG_HIPD[%ld] = false;\n", JScnt);

	   if(star[k].FLAG_RECONS==1)
	     fprintf(fpjs, "catalogFLAG_RECONS[%ld] = true;\n", JScnt);
	   else
	     fprintf(fpjs, "catalogFLAG_RECONS[%ld] = false;\n", JScnt);

	   fprintf(fpjs, "catalogPlanetPhotonRate[%ld] = %f;\n", JScnt, Pphot);
	   fprintf(fpjs, "catalogStarLeakPhotonRate[%ld] = %f;\n", JScnt, Sphot);
	   fprintf(fpjs, "catalogZodiPhotonRate[%ld] = %f;\n", JScnt, Zphot);
	   fprintf(fpjs, "catalogExoZodiPhotonRate[%ld] = %f;\n", JScnt, EZphot);
	   fprintf(fpjs, "catalogSNR1hr[%ld] = %f;\n", JScnt, SNR*60.0);
	   fprintf(fpjs, "catalogName[%ld] = \"%s\";\n", JScnt, star[k].name);
	   JScnt++;
	 }
    }
  fclose(fp);
  if(WRITEJS==1)
    fclose(fpjs);

    }

  printf("------------------------------------------------------------------------------\n");
  
  k = 64163;
  printf("%s %s %s\n", star[k].name1, star[k].name2, star[k].name3);
  printf("Coord         : %f %f\n", star[k].RAdeg, star[k].DECdeg);
  printf("dist          : %f   err: %f\n", star[k].dist, star[k].disterr);

  printf("\n");
  printf("FLAG_SUPERBLINK  : %d\n", star[k].FLAG_SUPERBLINK);
  printf("FLAG_2MASSCNS3   : %d\n", star[k].FLAG_2MASSCNS3);
  printf("FLAG_CNS3        : %d\n", star[k].FLAG_CNS3);
  printf("FLAG_HIPS        : %d\n", star[k].FLAG_HIPS);  
  printf("FLAG_HIPD        : %d\n", star[k].FLAG_HIPD);
  printf("\n");

  printf("VI            : %f\n", star[k].VI);
  printf("VJ            : %f\n", star[k].VJ);
  printf("\n");

  printf("BTmag         : %f\n", star[k].BTmag);
  printf("VTmag         : %f\n", star[k].VTmag);
  printf("BJmag         : %f\n", star[k].BJmag);
  printf("RFmag         : %f\n", star[k].RFmag);
  printf("INmag         : %f\n", star[k].INmag);
  printf("Vemag         : %f\n", star[k].Vemag);
  printf("Jmag          : %f\n", star[k].Jmag);
  printf("Hmag          : %f\n", star[k].Hmag);
  printf("Kmag          : %f\n", star[k].Kmag);
  printf("SpectralIndex : %f\n", star[k].SpectralIndex);

  printf("---------- CNS3 ---------\n");
  printf("Name          : %s %s %s\n", star[k].CNS3name1, star[k].CNS3name2, star[k].CNS3name3);
  printf("SpType        : %f\n", star[k].SpType);
  printf("mV            : %f\n", star[k].mV);
  printf("BV            : %f\n", star[k].BV);
  printf("UB            : %f\n", star[k].UB);
  printf("RI            : %f\n", star[k].RI);
  printf("CNS3dist      : %f\n", star[k].CNS3dist);
  printf("MV            : %f\n", star[k].MV);

  printf("-------- HIPPARCOS ------\n");
  printf("HIPindex      : %ld\n", star[k].HIPindex);
  printf("HIPVmag       : %f\n", star[k].HIPVmag);
  printf("HIPSpType     : %f\n", star[k].HIPSpType);
  printf("HIPBT         : %f\n", star[k].HIPBT);
  printf("HIPVT         : %f\n", star[k].HIPVT);
  printf("HIPBV         : %f\n", star[k].HIPBV);
  printf("HIPVI         : %f\n", star[k].HIPVI);

  return(0);
}






