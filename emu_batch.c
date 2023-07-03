#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logk.h"
#include "params.h"
#include "gsl/gsl_spline.h"
#include "cpgplot.h"

int main(int argc, char**argv)
{


  nparams = 5;

  int  npred = 5;
  char inputs[256];
  char paramnames[5][20];
  double outputredshift = 0;


  sprintf(paramnames[0], "M\\dcut\\u");
  sprintf(paramnames[1], "M\\d1\\u");
  sprintf(paramnames[2], "\\gs");
  sprintf(paramnames[3], "\\gk");
  sprintf(paramnames[4], "\\ga");

  double inputparams[5][5] =
    {
      {0.1, 0.1, 0.1, 0.1, 0.1},
      {0.25, 0.25, 0.25, 0.25, 0.25},
      {0.5, 0.5, 0.5, 0.5, 0.5},
      {0.75, 0.75, 0.75, 0.75, 0.75},
      {0.8, 0.8, 0.8, 0.8, 0.8},
    };


  int param_no = 0;

  int m, n;
  float  *k = malloc(nk*sizeof(float));
  double *output_pk = malloc(nk*sizeof(double));
  float *Pk = malloc(npred*nk*sizeof(float));

  for (m = 0; m < npred; m++)
    {

      emu(&(inputparams[m]), outputredshift, output_pk);

      char filename[256];
      if (m == 4)
	sprintf(filename, "holdouts/pk_M000_005_bspline_10_1.0000_ratio_test.dat");
      else
	sprintf(filename,"holdouts/pk_M000_%03d_bspline_10_1.0000_ratio_test.dat", m);

      FILE *fp = fopen(filename, "r");
      

      int nread = 2;
      float junk, input_pk;
      for (n = 0; n < nk; n++)
	{
	  //	  k[n] = (float)logk[n]; //(float)pow(10.,logk[n]);
	  fscanf(fp, "%f %f", &junk, &input_pk);
	  input_pk = log10(input_pk/pow(10.,junk));
	  Pk[npred*n+m] = output_pk[n]/input_pk;//(float)((pow(10.,output_pk[n])/pow(pow(10.,k[n]), 1.5)*4.*M_PI*M_PI));
	  /* fprintf(stderr, "%f %f\n", k[n], Pk[n]);  */
	}

      fclose(fp);

    }




/*   FILE *fpout = fopen("emu_test_z1.dat", "w"); */

/*   for (n = 0; n < nk; n++) */
/*      fprintf(fpout, "%f %f %f %f %f %f\n", (float)logk[n], Pk[npred*n+0],Pk[npred*n+1],Pk[npred*n+2],Pk[npred*n+3],Pk[npred*n+4]); */

  cpgopen("?");
  cpgslw(5);
  cpgsch(1.2);
  cpgenv(-2., 0., 0.978, 1.02, 0, 10);
  cpglab("k [Mpc\\u-1\\d]", "Emu/N-body","");

  float *kplot = malloc(nk*sizeof(float)); 
  float *pkplot = malloc(nk*sizeof(float)); 

  for (m = 0; m < 5; m++)
    {
      for (n = 0; n < nk; n++)
	{
	  kplot[n] = (float)(logk[n]); 
	  pkplot[n] = (float)Pk[npred*n+m]; 
	}
      cpgsci(m+2);
      cpgline(nk,kplot, pkplot);
    }

  cpgsci(1);
  cpgmove(-2.,1);
  cpgdraw(0.,1);

  cpgsls(2);
  cpgmove(-2.,1.01);
  cpgdraw(0.,1.01);

  cpgmove(-2.,0.99);
  cpgdraw(0.,0.99);


  cpgend();

/*   fclose(fpout);  */
  free(output_pk);
  free(k); free(Pk);


}


