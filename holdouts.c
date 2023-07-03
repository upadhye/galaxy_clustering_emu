#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "logk.h"
#include "params.h"
#include "cpgplot.h"
#include "gsl/gsl_spline.h"

int main(int argc, char**argv)
{
  double params[5][5] = {
    0.6363, 0.5333, 0.8571, 0.70, 0.70,
    0.2500, 0.2500, 0.2500, 0.25, 0.25,
    0.5000, 0.5000, 0.5000, 0.50, 0.50,
    0.7500, 0.7500, 0.7500, 0.75, 0.75,
    0.4545, 0.4667, 0.4286, 0.40, 0.40}; 

  double scalefactor[6] = {0.4985, 0.6086, 0.6974, 0.8051, 0.9083, 1.0000}; 

  double *pk_emu = malloc(nk*sizeof(double)); 

  float *k_nbody = malloc(nk*sizeof(float)); 
  float *pk_nbody = malloc(nk*sizeof(float)); 

  float *k_plot = malloc(nk*sizeof(float)); 
  float *pk_plot = malloc(nk*sizeof(float)); 
  //  float *ratio = malloc(nk*sizeof(float)); 

  cpgopen("?"); 
  cpgsch(1.2); 
  cpgslw(4); 
  cpgenv(-2., 0., 0.95, 1.05, 0, 10); 
  cpglab("k [Mpc\\u-1\\d]", "Emu/N-body", ""); 

  int i=0, j=0, k=0;   

  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nk); 
  gsl_interp_accel *acc  = gsl_interp_accel_alloc ();


  for (i = 0; i < 6; i++)
    {
      double outputredshift = 1/scalefactor[i] - 1;
      /* double outputredshift = 0; */



      for (j = 0; j < 5; j++)
	{
	  char filename[256]; 
	  sprintf(filename, "holdouts/pk_M000_%03d_bspline_10_%.4f_ratio_test.dat", j, scalefactor[i]); 
	  fprintf(stderr, "%s\n", filename); 
	  FILE *fp = fopen(filename, "r"); 


	  float junk; 
	  for (k = 0; k < 330; k++)
	    {
	      fscanf(fp, "%f %f", &(k_nbody[k]), &(pk_nbody[k])); 
	      //	      pk_nbody[k] = (pk_nbody[k]); 
	      //	      fprintf(stderr, "%f %f\n", k_nbody[k], pk_nbody[k]); 
	    }

	  fclose(fp); 

	  // remember to use version of emu that only outputs the bias
	  //	  double inputparams[5] = {0.5, 0.5, 0.5, 0.5, 0.5};
	  /* for (k = 0; k < 5; k++) */
	  /*   inputparams[k] = 0.5; */

	  //	  fprintf(stderr, "%lf %lf %lf %lf %lf\n", params[0][0], params[0][1], params[0][2], params[0][3], params[0][4]); 

	  emu(params[j], outputredshift, pk_emu); 

	  gsl_spline_init(spline, logk, pk_emu, nk); 


	  for(k = 0; k < 329; k++)
	    {
	      k_plot[k] = k_nbody[k];
	      double pk_emu_interp = gsl_spline_eval(spline, k_nbody[k], acc);
	      pk_plot[k] = (pk_emu_interp)/(pk_nbody[k]); 
	      //	      fprintf(stderr, "%f %lf %f %f\n", k_nbody[k], pk_emu_interp, pk_nbody[k], pk_plot[k]); 
	      /* if (k >= 331) */
	      /* 	pk_plot[k] = (float)log10(pk_emu[k])/pk_nbody[k-331]; */
	      /* else */
	      /* 	pk_plot[k] = 1; */
	      
	    }
	  cpgsci(1+j);
	  cpgline(329, k_nbody, pk_plot);
	  /* cpgsci(3); */
	  /* cpgline(330, k_nbody, pk_plot); */

	  gsl_interp_accel_reset(acc);
	}

    }

  cpgend(); 

  free(pk_emu); 
  free(k_plot); free(pk_plot); 
  free(k_nbody); free(pk_nbody); 
  gsl_spline_free(spline); gsl_interp_accel_free(acc); 
  return(0); 
}
