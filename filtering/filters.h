#ifndef filters_h
#define filters_h
// log_stock_prices are the log of stock prices
// muS is the real-world stock drift
// n_stock_prices is the number of the above stock prices
// (omega, theta, xi, rho) are the Heston parameters
// u[] is the set of means of observation errors
// v[] is the set of variances of observation errors
// estimates[] are the estimated observations from the filter
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

#endif