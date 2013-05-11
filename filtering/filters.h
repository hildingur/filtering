#ifndef filters_h
#define filters_h

/**
This file contains the function headers for model parameter estimators 
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
void estimate_extended_kalman_parameters_1_dim(
	double *log_stock_prices,
	double muS,
	int n_stock_prices,
	double omega,
	double theta,
	double xi,
	double rho,
	double *u,
	double *v,
	double *estimates);

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
void estimate_unscented_kalman_parameters_1_dim(
                                                double *log_stock_prices,
                                                double muS,
                                                int n_stock_prices,
                                                double omega,
                                                double theta, double xi,
                                                double rho,
                                                double *u,
                                                double *v,
                                                double *estimates);

/*
void estimate_particle_extended_kalman_parameters_1_dim(
                                                        double *log_stock_prices,
                                                        double muS,
                                                        int n_stock_prices,
                                                        double omega,
                                                        double theta,
                                                        double xi,
                                                        double rho,
                                                        double *ll,
                                                        double *estimates);

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
                                                         double *estimates);
														 */

#endif
