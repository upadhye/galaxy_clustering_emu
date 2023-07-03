#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "logk.h"
#include "params.h"
#include "cpgplot.h"

int main(int argc, char**argv)
{

  double inputp[nparams]; // upper limit for derivative
  double inputm[nparams]; // lower limit for derivative
  double mid_design[nparams]; 
  char inputs[256]; 
  char paramnames[5][20]; 
  double outputredshift = 0.5; 

  sprintf(paramnames[0], "M\\dcut\\u"); 
  sprintf(paramnames[1], "M\\d1\\u"); 
  sprintf(paramnames[2], "\\gs"); 
  sprintf(paramnames[3], "\\gk"); 
  sprintf(paramnames[4], "\\ga"); 


  /* FILE *fpderiv = fopen("pk_derivs_z0.dat","w"); */
  double *deriv = malloc(sizeof(double)*nk*5);
  double *plus = malloc(nk*sizeof(double)); 
  double *minus = malloc(nk*sizeof(double)); 

  float *k = malloc(nk*sizeof(float)); 
  float *dPdtheta = malloc(nk*sizeof(float)); 

  int m, n; 

  /* float min_design[5] = {12.85, 13.3, 0.5, 0.5, 0.5}; */
  /* float max_design[5] = {13.85, 14.3, 1.2, 1.5, 1.5}; */

  for (n = 0; n < 5; n++)
    {
      mid_design[n] = 0.5;      
    }

  for (n = 0; n < nk; n++)
    k[n] = (float)logk[n];
  
  cpgopen("?");
  cpgsch(1.4);
  cpgslw(3);
  cpgenv(-2, 0., -2, 2, 0, 10); 
  cpglab("k [Mpc\\u-1\\d]", "\\(0683)logP/\\(0683)\\gh\\di\\u",""); 

  int i = 0; 
  for (i = 0; i < 2; i++)
    {
      outputredshift = 0.5*i+0.5; 
      for (m = 0; m < 5 ; m++)
	{
	  for (n = 0; n < 5; n++)
	    {
	      inputp[n] = mid_design[n];
	      inputm[n] = mid_design[n];
	    }


	  inputp[m] = mid_design[m]*1.01; //1% change for the derivative
	  inputm[m] = mid_design[m]*0.99;
	  /* if (m==1) */
	  /* 	{ */
	  /* 	  inputp[m] = midpoint[m]*1.02; //2% change for the derivative */
	  /* 	  inputm[m] = midpoint[m]*0.98; */
	  /* 	} */
	  //      fprintf(stderr, "inputp = %f, inputm=%f\n", inputp[m], inputm[m]);
	  emu(inputp, outputredshift, plus);
	  emu(inputm, outputredshift, minus);

      
	  for (n = 0; n < nk; n++)
	    {
	      float k_now = pow(10.,logk[n]); 
	      //	      plus[n] = pow(10.,plus[n])/pow(k_now,1.5)*4.*M_PI*M_PI;
	      plus[n] = log10(plus[n]); 
	      //	      minus[n] = pow(10.,minus[n])/pow(k_now,1.5)*4.*M_PI*M_PI;
	      minus[n] = log10(minus[n]); 

	      deriv[m*nk+n] = (plus[n]-minus[n])/(inputp[m]-inputm[m]);
	      dPdtheta[n] = (float)deriv[m*nk+n]; 
	    }
      
	  /* cpgsci(2+m); */
	  if (i==0)
	    cpgsci(15);
	  else 
	    cpgsci(1);


	  cpgsls(m+1);
	  cpgslw(5);
	  cpgline(nk, k, dPdtheta); 

	  cpgsci(1);
	  if (i==0)
	    {
	      cpgmove(log10(0.013), 1.75-0.2*m);
	      cpgdraw(log10(0.018), 1.75-0.2*m);

	      cpgsci(1);
	      cpgslw(1);
	      cpgsls(1);

	      cpgtext(log10(0.02), 1.75-0.2*m, paramnames[m]);
	    }

	}

    }  
  cpgsci(15);
  cpgtext(log10(0.04), 1.5, "z = 0");
  cpgsci(1);
  cpgtext(log10(0.04), 1.25, "z = 0.5"); 
  /* cpgsci(1); */
  /* cpgtext(log10(0.04), 1.0, "z = 1"); */

  for (n = 0; n < nk; n++)
    {
      fprintf(stderr, "%f %f %f %f %f %f\n", logk[n], deriv[n], deriv[n+1*nk], deriv[n+2*nk], deriv[n+3*nk], deriv[n+4*nk]);
    }

  cpgend(); 
  
  /* fclose(fpderiv); */



  return(0);
}
