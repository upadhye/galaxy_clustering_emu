#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logk.h"
#include "params.h"
#include "cpgplot.h"

int read_colour_table(int option); 

int main(int argc, char**argv)
{
  nparams = 5;
  double inputparams[nparams]; 
  double inputp[nparams]; // upper limit for derivative
  double inputm[nparams]; // lower limit for derivative
  double mid_design[nparams]; 
  char inputs[256]; 
  char paramnames[5][20]; 
  double outputredshift; 


  sprintf(paramnames[0], "M\\dcut\\u"); 
  sprintf(paramnames[1], "M\\d1\\u"); 
  sprintf(paramnames[2], "\\gs"); 
  sprintf(paramnames[3], "\\gk"); 
  sprintf(paramnames[4], "\\ga"); 


  // this takes care of the midpoint which is used as a baseline for the comparison
  int l, m, n; 
  for (n = 0; n < 5; n++)
      mid_design[n] = 0.5;


  /* midpoint[0] = 0.1296;  */
  /* midpoint[1] = 0.0225;  */
  /* midpoint[2] = 0.95;  */
  /* midpoint[3] = 1.0; */
  /* midpoint[4] = 0.8; */

  outputredshift = 0.; 

  double *midpoint = malloc(nk*sizeof(double)); 
  double *var = malloc(nk*6*sizeof(double)); 
  double *plus = malloc(nk*sizeof(double)); 
  double *minus = malloc(nk*sizeof(double)); 

  emu(mid_design, outputredshift, midpoint); 

  // now do the other models to compare with it

  for (n = 0; n < 5; n++)
    inputparams[n] = 0.5;

  char filename[256]; 

  cpgopen("?"); 
  cpgpap(12., 0.25);
  cpgsch(2);
  cpgslw(2);
  // read_colour_table(2); 
  cpgscr(16, 0.878,0.878,0.878);
  cpgscr(17, 0.627,0.627,0.627);
  cpgscr(18, 0.5,0.5,0.5);
  cpgscr(19, 0.376, 0.376, 0.376); 
  cpgscr(20, 0.25, 0.25, 0.25); 

  cpgsci(1); 
  for (n = 0; n < 5; n++)
    {
      float xmin, xmax; 
      float ymin, ymax; 
      xmin = 0.07+n*0.182; 
      xmax = 0.07+(n+1)*0.182; 
      ymin = 0.22; 
      ymax = 0.95; 
      cpgsvp(xmin,xmax,ymin,ymax);      
      if (n==0)
	{
	  cpgswin(-2.,0,2.7,6.); 
	  //	  cpgswin(12.001,15.,-1.,1.); 
	  cpgbox("BCSTNL", 0.0, 0., "BCSTNL", 0.2, 0);
	  cpglab("", "P(k)",""); 
	  

	}
      else
	{
	  cpgswin(-2.,0.,3.,6.); 
	  cpgbox("BCSTNL", 0.0, 0., "BCSTL", 0.2, 0);
	}

      if (n==2)
	{
	  cpgtext(-1.,2.5,"k"); 
	}

      //  n = 0; 
      for (m = 0; m < 10 ; m++)
	{
	  inputparams[n] = 0.1*(m+0.5); //there are five intervals between max and min
	  //	  emu(inputparams, outputredshift, &(cM_var[m*nM]), &h);
	  emu(inputparams, outputredshift, var);
	  //	  fprintf(stderr, "%f %f %f %f %f\n", inputparams[0], inputparams[1], inputparams[2], inputparams[3], inputparams[4]); 
	  inputparams[n] = 0.1*(m+1); 
	  emu(inputparams, outputredshift, plus);
	  //	  fprintf(stderr, "%f %f %f %f %f\n", inputparams[0], inputparams[1], inputparams[2], inputparams[3], inputparams[4]); 
	  inputparams[n] = 0.1*(m); 
	  emu(inputparams, outputredshift, minus);
	  //	  fprintf(stderr, "%f %f %f %f %f\n", inputparams[0], inputparams[1], inputparams[2], inputparams[3], inputparams[4]); 
	  float logkforplot[2*nk], pkforplot[2*nk], pkforplotp[nk], pkforplotm[nk]; 
	  int nmodes = 0, nmodes2=0; 

	  for (l = 0; l < nk; l++)
	    {
	      /* if (logk[l] > -2 && logk[l] < 0) */
		{
		  logkforplot[nmodes] = (float)logk[l];
		  pkforplot[nmodes] = (float)log10(var[l]);//(log10(pow(10.,var[l])/pow(pow(10.,logk[l]), 1.5)*4.*M_PI*M_PI*M_PI)); 
		  nmodes++;
		}
	    }

	  /* fprintf(stderr, "%d\n", nmodes); */
	  /* for (l = 0; l < nk; l++) */
	  /*   if (logk[nk-l-1] > -2 && logk[nk-l-1] < 0) */
	  /*     { */
	  /* 	logkforplot[nmodes+nmodes2] = (float)logk[nk-l-1]; */
	  /* 	pkforplot[nmodes+nmodes2] = (float)(minus[nk-l-1]/midpoint[nk-l-1]); */
	  /* 	nmodes2++; */
	  /* 	//		fprintf(stderr, "%f %f %f %f\n", logk[nk-l-1], minus[nk-l-1]/midpoint[nk-l-1], logkforplot[nmodes+l], pkforplot[nmodes+l]);  */
	  /*     } */

	  /* FILE *fp_test = fopen("testpoly.dat","w"); */
	  /* for (l = 0; l < nmodes+nmodes2; l++) */
	  /*   fprintf(fp_test, "%f %f\n", logkforplot[l], pkforplot[l]);  */
	  /* fclose(fp_test); */

	  /* for (l = 0; l < nM; l++) */
	  /*   { */
	  /*     logmforplot[l] = (float)logm[l];  */
	  /*     // 	      cMforplot[l] = (float)(cM_var[l+nM*m]-cM_midpoint[l]);  */
 	  /*     cMforplot[l] = (float)(cM_var[l]-cM_midpoint[l]);  */
 	  /*     cMforplotp[l] = (float)(cM_plus[l]-cM_midpoint[l]);  */
 	  /*     cMforplotm[l] = (float)(cM_minus[l]-cM_midpoint[l]);  */
	  /*     if (l%4==0) */
	  /* 	{ */
	  /* 	  cpgerr1(2,logmforplot[l], cMforplot[l], (float)(cM_plus[l]-cM_var[l]), 1.0); */
	  /* 	  cpgerr1(4,logmforplot[l], cMforplot[l], (float)(cM_var[l]-cM_minus[l]), 1.0); */
	  /* 	} */
	  /*     /\* fprintf(stderr, "%f %f\n", logm[l], cMforplot[l]);  *\/ */
	  /*   } */

	  /* for(l=0;l < nk; l++) */
	  /*   { */
	  /*     logkforplot[l] = (float)logk[l]; */
 	  /*     pkforplot[l] = (float)(plus[l]/midpoint[l]); */
	  /*     logkforplot[l+nk] = (float)logk[nk-l-1]; */
 	  /*     pkforplot[l+nk] = (float)(minus[nk-l-1]/midpoint[nk-l-1]); */
	  /*     pkforplotp[l] = (float)(var[l]/midpoint[l]);  */
	  /*   } */


	  cpgslw(4);
	  /* cpgline(nM, logmforplot, cMforplot);     */

	  /* cpgline(nM, logmforplot, cMforplotm);     */
	  cpgsfs(1);
	  cpgsci(16+m);
	  /* if (m%2) */
	  /*   cpgshs(45.0, 2.0, 0.); */
	  /* else */
	  /*   cpgshs(-45.0, 2.0, 0.); */
	  /* cpgpoly(nmodes+nmodes2, logkforplot, pkforplot); */
	  cpgsci(1);
	  cpgslw(3);
	  cpgline(nmodes, logkforplot, pkforplot);
	  cpgsci(1);
	  cpgslw(2);
	  cpgtext(-0.5,5.7,paramnames[n]); 

	}

      inputparams[n] = mid_design[n];  // return inputparams back to original midpoint params
      /* sprintf(filename, "sensitivities_%d.dat", n); */
      /* FILE *fp = fopen(filename, "w"); */

      /* for (l = 0; l < nM; l++) */
      /* 	fprintf(fp,"%f %f %f %f %f %f %f %f\n", logm[l], cM_midpoint[l], cM_var[l], cM_var[l+nM], cM_var[l+nM*2], cM_var[l+nM*3], cM_var[l+nM*4], cM_var[l+nM*5]); */

      /* fclose(fp); */




    }

  for(m=0;m<5;m++)
    {
      inputp[m] = mid_design[m];
      inputm[m] = mid_design[m];
    }

  FILE *fpderiv = fopen("pk_derivs_z0.dat","w");
  double *deriv = malloc(sizeof(double)*nk*5);
  for (m = 0; m < 5 ; m++)
    {
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
  	  deriv[m*nk+n] = (plus[n]-minus[n])/(inputp[m]-inputm[m]);
  	}
    }
  
  /*     // restore original settings */
  /*     inputp[m] = midpoint[m];  */
  /*     inputm[m] = midpoint[m];  */

  /*     /\* float logmforplot[nM], cMforplot[nM];  *\/ */
  /*     /\* for (l = 0; l < nM; l++) *\/ */
  /*     /\* 	{ *\/ */
  /*     /\* 	  logmforplot[l] = (float)logm[l];  *\/ */
  /*     /\* 	  cMforplot[l] = (float)(cM_var[l+nM*m]/cM_midpoint[l]);  *\/ */
  /*     /\* 	} *\/ */
  /*     /\* cpgsci(7+m); *\/ */
  /*     /\* cpgslw(4); *\/ */
  /*     /\* cpgline(nM, logmforplot, cMforplot);     *\/ */
  /*     /\* cpgsci(1); *\/ */
  /*     /\* cpgtext(14.25,0.75,paramnames[n]);  *\/ */
  /*     /\* cpgslw(2); *\/ */
  /*   } */

  for (n = 0; n < nk; n++)
    {
      fprintf(fpderiv, "%f %f %f %f %f %f\n", logk[n], deriv[n], deriv[n+1*nk], deriv[n+2*nk], deriv[n+3*nk], deriv[n+4*nk]);
    }
  
  
  fclose(fpderiv);
  cpgend();
  free(midpoint); free(var);  
  return(0);
}



int read_colour_table(int option)
{
  FILE *fp;
  int Nindices; 
  char filename[256];
  
  if (option == 1)
    {
      sprintf(filename,"/home/jkwan/Analysis/HEAOB.lut"); 
      Nindices=200;
    }
  else if (option == 2)
    {
      sprintf(filename,"/home/jkwan/Analysis/SLS.lut");
      Nindices=200;
    }
  else if (option == 3)
    {
      sprintf(filename, "/home/jkwan/Analysis/idl5.lut");
      Nindices=256;
    }
  else
    {
      fprintf(stderr, "Invalid table selected; must be between 1-3");
      return(1);
    }
  fp = fopen(filename,"r");
  int i,count; 
  float rgb[3];

  for (i = 0,count=0; i < Nindices; i++)
    {
      fscanf(fp, "%f %f %f\n", &(rgb[0]), &(rgb[1]), &(rgb[2])); 
      if (i%2 == 0)
	{
	  cpgscr(count, rgb[0], rgb[1], rgb[2]);
	  cpgsci(count++); 
	}
    }

  cpgscr(0, 1,1,1);
  cpgscr(1, 0,0,0);


  fclose(fp); 
  return(0);
}
