#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_spline.h>
#include"junk.h"

int main(int argc, char**argv)
{
  int nk = 661; 
  int i, j, k;
  double h = 0.71; 

  FILE *fp_k, *fp_pk;
  char intermediate_k_file[256]; 
  char intermediate_pk_file[256], pk_z0_file[256], pk_z1_file[256]; 

  sprintf(intermediate_k_file, "logk.txt"); 
  sprintf(intermediate_pk_file, "ymean.txt"); 
  sprintf(pk_z0_file, "pk_smooth_z0_new.dat"); 
  sprintf(pk_z1_file, "pk_smooth_z1_new.dat"); 


  // Read in new files from Earl

  fp_k = fopen(intermediate_k_file, "r"); 
  fp_pk = fopen(intermediate_pk_file, "r"); 

  double *k_int = malloc(1228*sizeof(double)); 
  double **pk_int = (double**)malloc(4*sizeof(double*));

  for (i = 0; i < 4; i++)
    pk_int[i] = malloc(1228*sizeof(double));
  

  int nk_int = 0; 
  int nread = 1; 
  while (nread == 1)
    {
      nread = fscanf(fp_k, "%lf", &(k_int[nk_int]));
      fscanf(fp_pk, "%lf %lf %lf %lf", &(pk_int[0][nk_int]),&(pk_int[1][nk_int]),&(pk_int[2][nk_int]),&(pk_int[3][nk_int]));
      double know = pow(10.,k_int[nk_int]);
      k_int[nk_int] = know*h;
      k_int[nk_int] = log10(k_int[nk_int]); 
      for (k = 0; k < 4; k++)
	{
	  pk_int[k][nk_int] = pow(10.,pk_int[k][nk_int]);

	  pk_int[k][nk_int] *= (2*M_PI*M_PI*pow(know,-1.5));
	  pk_int[k][nk_int] /= h*h*h;
	}

      //      fprintf(stderr, "%lf %lf %lf %lf %lf\n", k_int[nk_int], pk_int[0][nk_int], pk_int[1][nk_int] , pk_int[2][nk_int], pk_int[3][nk_int]);
	//log10(k^1.5*P(k)/2pi^2) in [Mpc/h]^1.5
      nk_int++;
    }
									
  fclose(fp_k); fclose(fp_pk);

  nk_int--;


  // Read in old files from Earl
  FILE *fp_z0 = fopen("pk_smooth_z0_new.dat", "r"); 
  FILE *fp_z1 = fopen("pk_smooth_z1_new.dat", "r"); 

  int nk_old = 582; 
  double  *k_z0 = (double*)malloc(nk_old*sizeof(double)); 
  double  *k_z1 = (double*)malloc(nk_old*sizeof(double)); 
  double *pk_z0 = (double*)malloc(nk_old*sizeof(double));
  double *pk_z1 = (double*)malloc(nk_old*sizeof(double));


  for (i = 0; i < nk_old; i++)
    {
      fscanf(fp_z0, "%lf %lf", &(k_z0[i]), &(pk_z0[i]));
      fscanf(fp_z1, "%lf %lf", &(k_z1[i]), &(pk_z1[i]));
      k_z0[i] = log10(k_z0[i]); 
      k_z1[i] = log10(k_z1[i]); 
    }
									
  fclose(fp_z0); fclose(fp_z1);


  double **pk = (double**)malloc(6*sizeof(double*));
  for (i = 0; i < 6; i++)
    pk[i] = malloc(nk*sizeof(double)); 

  // This interpolates the four new redshifts
  for (i = 0; i < 4; i++)
    {
      gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nk_int);
      gsl_interp_accel *acc  = gsl_interp_accel_alloc ();

      gsl_spline_init (spline, k_int, pk_int[i], nk_int);

      for (j = 0; j < nk; j++)
          pk[i+1][j] = gsl_spline_eval(spline, logk[j], acc);

      gsl_interp_accel_free(acc);
      gsl_spline_free(spline);
    }

  //This interpolates the two old redshifts
  gsl_spline *spline_z0 = gsl_spline_alloc(gsl_interp_cspline, nk_old);
  gsl_interp_accel *acc_z0  = gsl_interp_accel_alloc ();

  gsl_spline *spline_z1 = gsl_spline_alloc(gsl_interp_cspline, nk_old);
  gsl_interp_accel *acc_z1  = gsl_interp_accel_alloc ();

  gsl_spline_init (spline_z0, k_z0, pk_z0, nk_old);
  gsl_spline_init (spline_z1, k_z1, pk_z1, nk_old);

  for (j = 0; j < nk; j++)
    {
      pk[0][j] = gsl_spline_eval(spline_z1, logk[j], acc_z1);
      pk[5][j] = gsl_spline_eval(spline_z0, logk[j], acc_z0);
    }


  gsl_interp_accel_free(acc_z0);       gsl_interp_accel_free(acc_z1);
  gsl_spline_free(spline_z0);       gsl_spline_free(spline_z1);

  FILE *fp_out = fopen("morejunk.h","w");

  for (j = 0; j < nk; j++)
    fprintf(fp_out, "%lf, %lf, %lf, %lf, %lf, %lf,\n", pk[0][j], pk[1][j], pk[2][j] , pk[3][j], pk[4][j], pk[5][j]);

  /* for (j = 0; j < nk; j++) */
  /*   fprintf(fp_out, "%lf %lf %lf %lf %lf %lf %lf\n", logk[j], pk[0][j], pk[1][j], pk[2][j] , pk[3][j], pk[4][j], pk[5][j]); */

  fclose(fp_out);

  free(k_int); free(pk_int); free(pk); 
  free(k_z0); free(pk_z0); free(k_z1); free(pk_z1);

  return(0); 
}

