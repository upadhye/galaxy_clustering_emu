#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logk.h"
#include "params.h"
#include "cpgplot.h"
#include "gsl/gsl_spline.h"

int main(int argc, char**argv)
{
  nparams = 5;
  double inputparams[nparams]; 
  double inputp[nparams]; // upper limit for derivative
  double inputm[nparams]; // lower limit for derivative
  double mid_design[nparams]; 
  char inputs[256]; 
  char paramnames[5][20]; 
  double outputredshift = 0; 

  int nread = 2; 
  int nkL = 0; 
  double *kL = malloc(3000*sizeof(double)); 
  double *PkL = malloc(3000*sizeof(double)); 
  FILE *fp_lin = fopen("/users/astro/jkwan/coyote/GPemulate/coyoteplus/cosmomc/camb/newM000_matterpower_z1.dat", "r");
  /* FILE *fp_lin = fopen("/users/astro/jkwan/coyote/GPemulate/data/Mira/M000/pk_coyote_M000_1.0000.dat", "r"); */
  while(nread==2 && nkL < 3000)
    {
      nread=fscanf(fp_lin, "%lf %lf\n", &(kL[nkL]), &(PkL[nkL])); 
      kL[nkL] *= 0.71;
      PkL[nkL] /= pow(0.71,3);
      nkL++; 
    }
  nkL--; 
  fclose(fp_lin); 

  gsl_interp_accel *accL
    = gsl_interp_accel_alloc ();
  gsl_spline *splineL 
    = gsl_spline_alloc (gsl_interp_cspline, nkL);

  gsl_spline_init (splineL, kL, PkL, nkL);


  nread = 2;
  int nkNL = 0;
  double *kNL = malloc(3000*sizeof(double));
  double *PkNL = malloc(3000*sizeof(double));
  FILE *fp_nonlin = fopen("/users/astro/jkwan/coyote/GPemulate/data/Mira/M000/pk_coyote_M000_0.4985.dat", "r");
  while(nread==2 && nkNL < 3000)
    {
      nread=fscanf(fp_nonlin, "%lf %lf\n", &(kNL[nkNL]), &(PkNL[nkNL]));
      nkNL++;
    }
  nkNL--;

  fclose(fp_nonlin);

  gsl_interp_accel *accNL
    = gsl_interp_accel_alloc ();
  gsl_spline *splineNL
    = gsl_spline_alloc (gsl_interp_cspline, nkNL);

  gsl_spline_init (splineNL, kNL, PkNL, nkNL);


  sprintf(paramnames[0], "M\\dcut\\u"); 
  sprintf(paramnames[1], "M\\d1\\u"); 
  sprintf(paramnames[2], "\\gs"); 
  sprintf(paramnames[3], "\\gk"); 
  sprintf(paramnames[4], "\\ga"); 

  int l, m, n; 
  for (n = 0; n < 5; n++)
    {
      inputparams[n] = 0.5;      
      /* fprintf(stderr, "%s = %f %f %f\n", paramnames[n], mid_design[n], max_design[n], min_design[n]); */
    }

  cpgopen("?"); 
  cpgpap(8., 0.8);
  cpgsch(1.5);
  cpgslw(3);
  cpgenv(-2, 0., 0., 10., 0, 10);
  cpglab("k [Mpc\\u-1\\d]", "bias",""); 
  cpgslw(5);
  int param_no = 0;

  for (m = 0; m < 11; m++)
    {
      inputparams[param_no] = 0.1*(m-1);

      double *output_pk = malloc(nk*sizeof(double)); 

      emu(inputparams, outputredshift, output_pk); 

      float  *k = malloc(nk*sizeof(float)); 
      float *Pk = malloc(nk*sizeof(float)); 
      float *biasL = malloc(nk*sizeof(float)); 
      float *biasNL = malloc(nk*sizeof(float)); 


      for (n = 0; n < nk; n++)
	{
	  k[n] = (float)logk[n]; //(float)pow(10.,logk[n]); 
	  Pk[n] = (float)output_pk[n];
	  double pk_interpL   = gsl_spline_eval (splineL, pow(10.,logk[n]), accL);
	  double pk_interpNL   = gsl_spline_eval (splineNL, pow(10.,logk[n]), accNL);
	  biasL[n] = sqrt(Pk[n]/pk_interpL);
	  biasNL[n] = sqrt(Pk[n]/pk_interpNL);
	  //	  fprintf(stderr, "%f %f\n", k[n], Pk[n]);
	}

      cpgline(nk, k, biasL); 
      cpgsls(2); 
      cpgline(nk, k, biasNL);
      cpgsls(1);

      free(output_pk); 
      free(k); free(Pk); free(biasL); free(biasNL); 
    }
  cpgsls(1); 
  cpgmove(log10(0.015), 9.); 
  cpgdraw(log10(0.02), 9.); 
  cpgsls(1);
  cpgslw(1); 
  cpgsch(1.3); 
  cpgtext(log10(0.025), 9., "linear P\\dm\\u(k)");  
  cpgsch(1.5); 

  cpgsls(2); 
  cpgslw(5); 
  cpgmove(log10(0.015), 8.5); 
  cpgdraw(log10(0.02), 8.5); 
  cpgsls(1); 
  cpgslw(1);
  cpgsch(1.3); 
  cpgtext(log10(0.025), 8.45, "non-linear P\\dm\\u(k)");  
  cpgsch(1.5);

  cpgslw(3);
  cpgsci(2); 
  cpgsch(1.5);
  cpgtext(log10(0.015), 7.75, "z=1");  
  cpgslw(1); 
  cpgsci(1); 
  cpgsch(1); 

  gsl_spline_free (splineNL);
  gsl_interp_accel_free (accNL);

  gsl_spline_free (splineL);
  gsl_interp_accel_free (accL);


 cpgend(); 
}


