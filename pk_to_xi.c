#include <math.h>
#include "gsl/gsl_sf_bessel.h"
#include"gsl/gsl_spline.h"

int pk_to_xi(int nk, double *k, double *Pk, double *r_out, double *xi_out)
{
  int n; 
  int max_sum = 500; 

  double *phi     = malloc(max_sum*sizeof(double)); 
  double *phidash = malloc(max_sum*sizeof(double)); 

  double sigma_smooth = 1.;
  int i;
  for (i =0; i < nk; i++)
    Pk[i] *= exp(-k[i]*k[i]*sigma_smooth*sigma_smooth);

  gsl_interp_accel *acc 
    = gsl_interp_accel_alloc ();
  gsl_spline *spline 
    = gsl_spline_alloc (gsl_interp_cspline, nk);

  gsl_spline_init (spline, k, Pk, nk);

  int m; 

  for (m = 0; m < 200; m++)
    {
      double r = 10+m; // r = 80 is the limit for k = 1; 
      double sum = 0; 
      /* for (n=0; n < nk; n++) */
      /* 	Pk[n] /= r;  */

      for (n = 1, sum=0; n < max_sum; n++)
	{
	  double nu = 0.5; 
	  int nzero = n; 
	  double zero = (double)n; // zero of Bessel/PI
 
	  double h = 1./150; //120.
	  double hr = h*zero;
	  phi[n] = hr *tanh(0.5*M_PI * sinh(hr));
	  phidash[n] = 1-(tanh(0.5*M_PI*sinh(hr))*tanh(0.5*M_PI*sinh(hr)));
	  phidash[n] *= 0.5*(sinh(M_PI*sinh(hr))+M_PI*hr*cosh(hr));
	  /* phidash[n] = M_PI*0.5*hr*cosh(hr)/cosh(0.5*M_PI*sinh(hr))/cosh(0.5*M_PI*sinh(hr))+tanh(0.5*M_PI*hr);  */
	  /* double weights = gsl_sf_bessel_Ynu(0.5,M_PI*zero)/gsl_sf_bessel_Jnu(1.5,M_PI*zero); */
	  double weights = -1.*sqrt(2./(M_PI*M_PI*zero))*cos(M_PI*zero)/(gsl_sf_bessel_Jnu(1.5,M_PI*zero));
	  /* fprintf(stderr, "%lf %lf\n", gsl_sf_bessel_Ynu(1.5,M_PI*zero),gsl_sf_bessel_Jnu(1.5,M_PI*zero)); */
	  double x = M_PI*phi[n]/h; 
	  double pkterm   = pow(x,1.5)*gsl_spline_eval (spline, x/r, acc);
	  double besselterm = gsl_sf_bessel_Jnu(0.5,x);
	  double hr_p = h*(zero+1.); 
	  double hr_m = h*(zero-1.); 
	  double hr_pp = h*(zero+2.); 
	  double hr_mm = h*(zero-2.); 

	  double phi_p = hr_p *tanh(0.5*M_PI * sinh(hr_p));
	  double phi_m = hr_m *tanh(0.5*M_PI * sinh(hr_m));
	  double phi_pp = hr_pp *tanh(0.5*M_PI * sinh(hr_pp));
	  double phi_mm = hr_mm *tanh(0.5*M_PI * sinh(hr_mm));
	  //	  double derivs = (phi[n]-phi[n-1])/(h); 
	  double derivs = (-phi_pp+8*phi_p-8*phi_m+phi_mm)/(12.*h); 
	  sum += M_PI*weights*pkterm*besselterm*derivs;
	  /* fprintf(stderr, "%d %f %f %f %f %f %f \n", n, sum, pkterm, besselterm, weights,  phi[n], derivs); */
	  /* fprintf(stderr, "%d %f %f %f %f \n", n, x, x/r, phi[n],  phidash[n] ); */
	}
      double constant = sqrt(0.5*M_PI)/(2.*M_PI*M_PI)/(r*r*r);
      
      r_out[m] = r; 
      xi_out[m] = sum*constant; 

    }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);


  free(phi); free(phidash);

  return(0); 
}
