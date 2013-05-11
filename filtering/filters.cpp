#include <math.h>
#include "cholesky.h"

/**
   This file contains the function parameters for model parameter estimators 
   as defined in Alireza Javaheri's book, "Inside Volatility Arbitrage"
**/

/*
  log_stock_prices are the log of stock prices
  muS is the real-world stock drift
  n_stock_prices is the number of the above stock prices
  (omega, theta, xi, rho) are the Heston parameters
  u[] is the set of means of observation errors
  v[] is the set of variances of observation errors
  estimates[] are the estimated observations from the filter
*/
void estimate_extended_kalman_parameters_1_dim (
                                                double *log_stock_prices,
                                                double muS,
                                                int n_stock_prices,
                                                double omega,
                                                double theta,
                                                double xi,
                                                double rho,
                                                double *u,
                                                double *v,
                                                double *estimates)
{
  int i1;
  double x, x1, W, H, A;
  double P, P1, z, U, K;
  const double delt=1.0/252.0;
  const double eps=0.00001;
  x = 0.04;
  P=0.01;
  u[0]=u[n_stock_prices-1]=0.0;
  v[0]=v[n_stock_prices-1]=1.0;
  estimates[0]=estimates[1]=log_stock_prices[0]+eps;

  for (i1=1;i1<n_stock_prices-1;i1++)
	{
      if (x<0) 
        x=0.00001;
		
      x1 = x + (omega-rho*xi*muS - (theta-0.5*rho*xi) * x) * delt + rho*xi* (log_stock_prices[i1]-log_stock_prices[i1-1]);
      A = 1.0-(theta-0.5*rho*xi)*delt;
      W = xi*sqrt((1-rho*rho) * x * delt);
      P1 = W*W + A*P*A;
      if (x1<0) 
        x1=0.00001;
      H = -0.5*delt;
      U = sqrt(x1*delt);
      K = P1*H/( H*P1*H + U*U);
      z = log_stock_prices[i1+1];
      x = x1 + K * (z - (log_stock_prices[i1] + (muS-0.5*x1)*delt));
      u[i1] = z - (log_stock_prices[i1] + (muS-0.5*x1)*delt);
      v[i1] = H*P1*H + U*U;
      estimates[i1+1] = log_stock_prices[i1] + (muS-0.5*x1)*delt;
      P=(1.0-K*H)*P1;
	}
}
// Having u[] and v[] we can evaluate the (minus log) Likelihood as
// the sum of log(v[i1])+u[i1]*u[i1]/v[i1]
// and minimize the sum in order to obtain the optimal parameters
// the minimization could be done for instance via the direction set
// routine available in the Numerical Recipes in C

void estimate_unscented_kalman_parameters_1_dim(double *log_stock_prices,
                                                double muS,
                                                int n_stock_prices,
                                                double omega,
                                                double theta, double xi,
                                                double rho,
                                                double *u,
                                                double *v,
                                                double *estimates)
{
  int     i1,i2, i3, t1;
  int     ret;
  int     na=3;
  double  x, xa[3];
  double  X[7], Xa[3][7];
  double  Wm[7], Wc[7], Z[7];
  double  x1;
  double  prod, prod1;
  double  P, P1;
  double **Pa, **proda;
  double  z, U, Pzz, K;
  double  delt=1.0/252.0;
  double  a=0.001 , b=0.0, k=0.0, lambda;
  double eps=0.00001;
  lambda = a*a*(na +k)-na;
  proda= new double * [na];
  Pa =   new double * [na];
  for (i1=0;i1<na;i1++)
    {
      Pa[i1]= new double [na];
      proda[i1]= new double [na];
    }
  xa[1]=xa[2]=0.0;
  x= 0.04;
  u[0]=u[n_stock_prices-1]=0.0;
  v[0]=v[n_stock_prices-1]=1.0;
  estimates[0]=estimates[1]=log_stock_prices[0]+eps;
  xa[0]=x;
  Pa[0][0]= Pa[1][1]= Pa[2][2] = 1.0;
  Pa[1][0]= Pa[0][1]= Pa[1][2]=Pa[2][1]= Pa[0][2]=Pa[2][0]=0;
  for (i1=0;i1<na;i1++)
    {
      for (i2=0;i2<na;i2++)
        {
          proda[i1][i2]=0.0;
        }
    }
  Wm[0]=lambda/(na+lambda);
  Wc[0]=lambda/(na+lambda) + (1-a*a+b);
  for (i3=1;i3<(2*na+1);i3++)
    {
      Wm[i3]=Wc[i3]=1/(2*(na+lambda));
    }
  for (t1=1;t1<n_stock_prices-1;t1++)
    {
      for (i1=0;i1<na;i1++)
        {
          Xa[i1][0]= xa[i1];
        }
      for (i1=0;i1<na;i1++)
        {
          for (i2=0;i2<na;i2++)
            {
              if (i1==i2)
                {
                  if (Pa[i1][i2] < 1.0e-10)
                    Pa[i1][i2]= 1.0e-10;
                } else {
                if (Pa[i1][i2] < 1.0e-10)
                  Pa[i1][i2]= 0.0;
              }
            }
        }
      ret = sqrt_matrix(Pa, proda, na);
      for (i3=1;i3<(1+na);i3++)
        {
          for (i1=0;i1<na;i1++)
            {
              Xa[i1][i3]= xa[i1] + sqrt(na+lambda) * proda[i1][i3-1];
            }
        }
      for (i3=(1+na);i3<(2*na+1);i3++)
        {
          for (i1=0;i1<na;i1++)
            {
              Xa[i1][i3]= xa[i1] - sqrt(na+lambda) * proda[i1][i3-na-1];
            } }
      for (i3=0;i3<(2*na+1);i3++)
        {
          if (Xa[0][i3]<0) Xa[0][i3]=0.0001;
          X[i3]= Xa[0][i3] + (omega-muS*rho*xi
                              -
                              (theta-0.5*rho*xi) *Xa[0][i3])*delt +
            rho*xi* (log_stock_prices[t1]-log_stock_prices[t1-1]) +
            xi*sqrt((1-rho*rho)*delt*Xa[0][i3])*Xa[1][i3];
        }
      x1 = 0;
      for (i3=0;i3<(2*na+1);i3++)
        {
          x1 += Wm[i3]*X[i3];
        }
      P1=0.0;
      for (i3=0;i3<(2*na+1);i3++)
        {
          P1 += Wc[i3]*(X[i3]-x1)*(X[i3]-x1);
        }
      z=0;
      for (i3=0;i3<(2*na+1);i3++)
        {
          if (X[i3]<0) X[i3]=0.00001;
          Z[i3] = log_stock_prices[t1] + (muS-0.5*X[i3])*delt +
            sqrt(X[i3]*delt)*Xa[2][i3];
          z += Wm[i3]*Z[i3];
        }
      Pzz=0;
      for (i3=0;i3<(2*na+1);i3++)
        {
          Pzz +=  Wc[i3]*(Z[i3]-z)*(Z[i3]-z);
        }
      prod=0.0;
      for (i3=0;i3<(2*na+1);i3++)
        {
          prod += Wc[i3]*(X[i3]-x1)* (Z[i3]-z);
        }
      K= prod/Pzz;
      u[t1] = log_stock_prices[t1+1] - z;
      v[t1] = Pzz;
      estimates[t1+1] = z;
      x = x1 + K*(log_stock_prices[t1+1] - z);
      P = P1 - K*K * Pzz;
      xa[0]=x;
      Pa[0][0] = P;
      if (x<0) x=0.0001;
      Pa[1][0]= Pa[0][1]= Pa[1][2]=Pa[2][1]= Pa[0][2]=Pa[2][0]=0;
    }
  for (i1=0;i1<na;i1++)
    {
      delete [] Pa[i1];
      delete [] proda[i1];
    }
  delete [] Pa;
  delete [] proda;
}

/*
// the routine sqrt_matrix() can be constructed via the Cholesly decomposition
// also available as choldc() in the Numerical Recipes in C
// log_stock_prices  are the log of stock prices
// muS is the real-world stock drift
// n_stock_prices is the number of the above stock prices
// (omega, theta, xi, rho) are the Heston parameters
// ll is the value of (negative log) Likelihood function
// estimates[] are the estimated observations from the filter
// The function ran2() is from Numerical Recipes in C
// and generates uniform random variables
// The function Normal_inverse() can be found from many sources
// and is the inverse of the Normal CDF
// Normal_inverse(ran2(.)) generates a set of Normal random variables
void estimate_particle_extended_kalman_parameters_1_dim(
                                                        double *log_stock_prices,
                                                        double muS,
                                                        int n_stock_prices,
                                                        double omega,
                                                        double theta,
                                                        double xi,
                                                        double rho,
                                                        double *ll,
                                                        double *estimates)
{
  int     i1, i2, i3;
  double  H, A, x0, P0, z;
  int     M=1000;
  double  x[1000], xx[1000], x1[1000], x2[1000];
  double  P[1000], P1[1000], U[1000], K[1000], W[1000];
  double  w[1000],  u[1000], c[1000];
  double  q, pz, px, s, m, l;
  double  delt=1.0/252.0, x1_sum;
  long    idum=-1;
  A = 1.0-(theta-0.5*rho*xi)*delt;
  H = -0.5*delt;
  x0 = 0.04;
  P0 = 0.000001;
  for (i2=0; i2<M; i2++)
    {
      x[i2] = x0 + sqrt(P0)* Normal_inverse(ran2(&idum));
      P[i2] = P0;
    }
  *ll=0.0;
  for (i1=1;i1<n_stock_prices-1;i1++)
    {
      l = 0.0;
      x1_sum=0.0;
      for (i2=0; i2<M; i2++)
        {
          // EKF for the proposal distribution 
          if (x[i2]<0) x[i2]=0.00001;
          x1[i2] = x[i2] + ( omega-rho*xi*muS - (theta-
                                                 0.5*rho*xi) * x[i2]) * delt + rho*xi*
            (log_stock_prices[i1]-log_stock_prices[i1-1]);
          W[i2]  = xi*sqrt((1-rho*rho) * x[i2] * delt);
          P1[i2] = W[i2]*W[i2] + A*P[i2]*A;
          if (x1[i2]<0) x1[i2]=0.00001;
          U[i2] = sqrt(x1[i2]*delt);
          K[i2] = P1[i2]*H/( H*P1[i2]*H + U[i2]*U[i2]);
          z = log_stock_prices[i1+1];
          x2[i2] = x1[i2] + K[i2] * (z - (log_stock_prices[i1]
                                          + (muS-0.5*x1[i2])*delt));
          x1_sum+= x1[i2];
          P[i2]=(1.0-K[i2]*H)*P1[i2];
          // sample 
          xx[i2] = x2[i2]+sqrt(P[i2])*Normal_inverse(ran2(&idum));
          if (xx[i2]<0) xx[i2]=0.00001;
          // calculate weights 
          m = x2[i2];
          s = sqrt(P[i2]);
          q = 0.39894228/s * exp( - 0.5* (xx[i2] - m)*
                                  (xx[i2] - m)/(s*s) );
          m = log_stock_prices[i1] + (muS-0.5*xx[i2])*delt;
          s = sqrt(xx[i2]*delt);
          pz = 0.39894228/s * exp( - 0.5* (z - m)*(z - m)/(s*s) );
          m = x[i2] + ( omega-rho*xi*muS - (theta-0.5*
                                            rho*xi) * x[i2]) * delt + rho*xi*
            (log_stock_prices[i1]-log_stock_prices[i1-1]);
          s = xi*sqrt((1-rho*rho) * x[i2] * delt);
          px = 0.39894228/s * exp( - 0.5* (xx[i2] - m)*
                                   (xx[i2] - m)/(s*s) );
          w[i2] = pz * px / MAX(q, 1.0e-10);
          l += w[i2]; }
      *ll += log(l);
      estimates[i1+1]= log_stock_prices[i1] +
        (muS-0.5*x1_sum/M)*delt;
      // normalize weights 
      for (i2=0; i2<M; i2++)
        w[i2] /= l;
      // resample and reset weights 
      c[0]=0;
      for (i2=1; i2<M; i2++)
        c[i2] = c[i2-1] + w[i2];
      i2=0;
      u[0] = 1.0/M * ran2(&idum);
      for (i3=0; i3<M; i3++)
        {
          u[i3] = u[0] + 1.0/M *i3;
          while (u[i3] > c[i2])
            i2++;
          x[i3] = xx[i2];
          w[i3] = 1.0/M;
        }
    }
  *ll *= -1.0; }
// *ll is the value of (negative log) Likelihood function
// we can minimize it to obtain the optimal parameter-set

void estimate_particle_unscented_kalman_parameters_1_dim(
                                                         double *log_stock_prices,
                                                         double muS,
                                                         int n_stock_prices,
                                                         double omega,
                                                         double theta,
                                                         double xi,
                                                         double rho,
                                                         double *ll,
                                                         double *estimates)
{
  int      i1, i2, i3, i4;
  int      na=3;
  double   x0, P0;
  double   Wm[7], Wc[7];
  int      M=1000;
  double   x[1000], xx[1000], x1[1000], x2[1000],
    zz[1000], Z[1000][7];
  double   X[1000][7], Xa[1000][3][7];
  double   xa[1000][3], prod[1000];
  double   P[1000], P1[1000], U[1000], K[1000],
    W[1000], Pzz[1000];
  double w[1000], u[1000], c[1000]; double ***Pa, ***proda;
  double q,pz,px,s,m, l, z;
  double delt=1.0/252.0;
  long idum=-1;
  int ret;
  double a=0.001 , b=0.0, k=0.0, lambda;
  proda= new double ** [M];
  Pa =   new double ** [M];
  for (i2=0;i2<M;i2++)
    {
      Pa[i2]= new double * [na];
      proda[i2]= new double * [na];
      for (i1=0;i1<na;i1++)
        {
          Pa[i2][i1]= new double [na];
          proda[i2][i1]= new double [na];
        }
    }
  for (i2=0;i2<M;i2++)
    {
      for (i1=0;i1<na;i1++)
        {
          for (i3=0;i3<na;i3++)
            {
              proda[i2][i1][i3]=0.0;
            } }
    }
  lambda = a*a*(na +k)-na;
  Wm[0]=lambda/(na+lambda);
  Wc[0]=lambda/(na+lambda) + (1-a*a+b);
  for (i3=1;i3<(2*na+1);i3++)
    {
      Wm[i3]=Wc[i3]=1/(2*(na+lambda));
    }
  x0 = 0.04;
  P0 = 0.000001;
  for (i2=0; i2<M; i2++)
    {
      x[i2] = x0 + sqrt(P0)* Normal_inverse(ran2(&idum));
      P[i2] = P0;
      xa[i2][0]=x[i2];
      xa[i2][1]=xa[i2][2]=0.0;
      Pa[i2][0][0]= P[i2];
      Pa[i2][1][1]= Pa[i2][2][2] = 1.0;
      Pa[i2][1][0]= Pa[i2][0][1]= Pa[i2][1][2] =
        Pa[i2][2][1] =
        Pa[i2][0][2] = Pa[i2][2][0] = 0.0;
    }
  *ll=0.0;
  for (i1=1;i1<n_stock_prices-1;i1++)
    {
      l = 0.0;
      estimates[i1+1]=0.0;
      for (i2=0; i2<M; i2++)
        {
          // UKF for the proposal distribution 
          for (i3=0;i3<na;i3++)
            {
              Xa[i2][i3][0]= xa[i2][i3];
            }
          for (i3=0;i3<na;i3++)
            {
              for (i4=0;i4<na;i4++)
                {
                  if (i3==i4)
                    {
                      if (Pa[i2][i3][i4] < 1.0e-10)
                        Pa[i2][i3][i4]= 1.0e-10;
                    } else {
                    if (Pa[i2][i3][i4] < 1.0e-10)
                      Pa[i2][i3][i4] = 0.0;
                  } }
            }
          ret = sqrt_matrix(Pa[i2],proda[i2],na);
          for (i3=1;i3<(1+na);i3++)
            {
              for (i4=0;i4<na;i4++)
                {
                  Xa[i2][i4][i3]= xa[i2][i4] + sqrt(na+lambda) *
                    proda[i2][i4][i3-1];
                } }
          for (i3=(1+na);i3<(2*na+1);i3++)
            {
              for (i4=0;i4<na;i4++)
                {
                  Xa[i2][i4][i3]= xa[i2][i4] - sqrt(na+lambda) *
                    proda[i2][i4][i3-na-1];
                } }
          for (i3=0;i3<(2*na+1);i3++)
            {
              if (Xa[i2][0][i3]<0) Xa[i2][0][i3]=0.0001;
              X[i2][i3]= Xa[i2][0][i3] + (omega-muS*rho*xi   -
                                          (theta-0.5*rho*xi) *Xa[i2][0][i3])*delt +
                rho*xi* (log_stock_prices[i1]-
                         log_stock_prices[i1-1]) +
                xi*sqrt((1-rho*rho)*delt*Xa[i2][0][i3])*
                Xa[i2][1][i3];
            }
          x1[i2] = 0;
          for (i3=0;i3<(2*na+1);i3++)
            {
              x1[i2] += Wm[i3]*X[i2][i3];
            }
          P1[i2]=0.0;
          for (i3=0;i3<(2*na+1);i3++)
            {
              P1[i2] += Wc[i3]*(X[i2][i3]-x1[i2])*(X[i2][i3]-
                                                   x1[i2]);
            }
          zz[i2]=0;
          for (i3=0;i3<(2*na+1);i3++)
            {
              if (X[i2][i3]<0) X[i2][i3]=0.00001;
              Z[i2][i3] = log_stock_prices[i1] +
                (muS-0.5*X[i2][i3])*delt + sqrt(X[i2][i3]*delt)*Xa[i2][2][i3];
              zz[i2] += Wm[i3]*Z[i2][i3];
            }
          Pzz[i2]=0;
          for (i3=0;i3<(2*na+1);i3++)
            {
              Pzz[i2] +=  Wc[i3]*(Z[i2][i3]-zz[i2])*(Z[i2][i3]-
                                                     zz[i2]);
            }
          prod[i2]=0.0;
          for (i3=0;i3<(2*na+1);i3++)
            {
              prod[i2] += Wc[i3]*(X[i2][i3]-x1[i2])* (Z[i2][i3]-
                                                      zz[i2]);
            }
          K[i2]= prod[i2]/Pzz[i2];
          z = log_stock_prices[i1+1];
          estimates[i1+1] += zz[i2]/M;
          x2[i2] = x1[i2] + K[i2]*(z - zz[i2]);
          P[i2] = P1[i2] - K[i2]*K[i2] * Pzz[i2];
          xa[i2][0]=x2[i2];
          Pa[i2][0][0] = P[i2];
          if (x2[i2]<0) x2[i2]=0.0001;
          Pa[i2][1][0]= Pa[i2][0][1]= Pa[i2][1][2]
            =Pa[i2][2][1]= Pa[i2][0][2]=Pa[i2][2][0]=[0];
          // sample 
          xx[i2] = x2[i2] + sqrt(P[i2])*
            Normal_inverse(ran2(&idum));
          if (xx[i2]<0) xx[i2]=0.00001;
          // calculate weights 
          m = x2[i2];
          s = sqrt(P[i2]);
          q = 0.39894228/s * exp( - 0.5* (xx[i2] - m)*
                                  (xx[i2] - m)/(s*s) );
          m= log_stock_prices[i1] + (muS-0.5*xx[i2])*delt;
          s= sqrt(xx[i2]*delt);
          pz= 0.39894228/s * exp( - 0.5* (z - m)*
                                  (z - m)/(s*s) );
          m= x[i2] + ( omega-rho*xi*muS -
                       (theta-0.5*rho*xi) * x[i2]) * delt +
            rho*xi* (log_stock_prices[i1]-
                     log_stock_prices[i1-1]);
          s= xi*sqrt((1-rho*rho) * x[i2] * delt);
          px= 0.39894228/s * exp( - 0.5* (xx[i2] - m)*
                                  (xx[i2] - m)/(s*s) );
          w[i2]= MAX(pz, 1.0e-10) *
            MAX(px, 1.0e-10) / MAX(q, 1.0e-10);
          l += w[i2];
        }
      *ll += log(l);
      // normalize weights 
      for (i2=0; i2<M; i2++)
        w[i2] /= l;
      // resample and reset weights 
      c[0]=0;
      for (i2=1; i2<M; i2++)
        c[i2] = c[i2-1] + w[i2];
      i2=0;
      u[0] = 1.0/M * ran2(&idum);
      for (i3=0; i3<M; i3++)
        {
          u[i3] = u[0] + 1.0/M *i3;
          while (u[i3] > c[i2])
            i2++;
          x[i3]= xx[i2];
          w[i3]=1.0/M;
        }
    }
  *ll *= -1.0;
  for (i2=0;i2<M;i2++)
    {
      for (i1=0;i1<na;i1++)
        {
          delete [] Pa[i2][i1];
          delete [] proda[i2][i1];
        }
    }
  for (i2=0;i2<M;i2++)
    {
      delete [] Pa[i2];
      delete [] proda[i2];
    }
  delete [] Pa;
  delete [] proda;
}
*/
