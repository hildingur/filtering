#include <math.h>

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