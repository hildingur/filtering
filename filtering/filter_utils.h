#ifndef filter_utils_h
#define filter_utils_h

#include <string>
#include <vector>
#include "recipes/nr.h"

using namespace std;

class vol_params 
{
private:
	double omega, theta, xi, roe, p, mle, chi2;
	double *u, *v, *estimates;
public:
	double get_omega();
	double get_theta();
	double get_xi();
	double get_roe();
	double get_p();
	double get_mle();
	double get_chi2();
	void set_omega(double omega_);
	void set_theta(double theta_);
	void set_xi(double xi_);
	void set_roe(double roe_);
	void set_p(double p_);
	void set_mle(double mle_);
	void set_chi2(double chi2_);
	void copy_params(vol_params& params);
	
	double* get_u();
	double* get_v();
	double* get_estimates();

	vol_params();
	vol_params(double omega_, double theta_, double xi_, double roe_, double p_, double mle_, double chi2_);
	vol_params(const vol_params& copy); //copy constructor
	/*
	Function that populates a vector with the initial parameters of the model
	*/
	void populate_starting_vector(Vec_IO_DP& starting_vector);
	void extract_params_from_vector(Vec_IO_DP& starting_vector);

	void set_best_estimate(const vol_params& params, 
		const double mle, 
		const double chi2, 
		const double* u,
		const double* v,
		const double* estimates,
		const int size);
};

void read_lines(string& fname, vector<double>& out);

int sqrt_matrix(double** pa, double** proda, int& N);

/**
Checks to see if the data is normal using the chi squared test with 20 degrees of freedom and alpha = 0.05
*/
bool is_normal(double* data, int data_count, double& chi2_);

double cumulative_normal(const double x);

double gaussrand();

void normalize_data(const double* data_in, double* data_out, const int length);

#endif