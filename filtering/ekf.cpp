#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "./recipes/nr.h"
#include "filters.h"

using namespace std;

void read_lines(char * fname, double * out);
DP minimize_target(Vec_I_DP & input);
DP minimize_target_w_x_y_z_square(Vec_I_DP & input);

DP minimize_target_extended_kalman_parameters_1_dim(Vec_I_DP & input);

const int n_stock_prices = 4617;
double log_stock_prices[n_stock_prices], u[n_stock_prices], v[n_stock_prices], estimates[n_stock_prices + 1];
//double muS = 0.687;
double muS = 0.0;
int call_counter = 0;

int main() {
	//TODO: Make the filename a command line arg, make the output a 
	//Command Line Argument and read the input into a vector
	//changed
	double z[n_stock_prices];
	read_lines("C:\\Users\\hildingur\\Dropbox\\vol_filter (2)\\ibm_trimmed.csv", z);
	for(int i = 0; i < n_stock_prices; i++) {
		log_stock_prices[i] = log(z[i]);
	}

	//Initializing the starting point
	double a[4] = {0.2, 1.00, 0.5, -0.2};
	Vec_IO_DP starting_point(a, 4);
	
	//Initializing the identity matrix, don't know a more elegant
	//way to do it using the Numerical recipies api.
	Mat_IO_DP identity_matrix(4, 4);
	for(int i = 0; i < 4; i++) {
		for(int j = 0; j < 4; j++) {
			if(i==j) identity_matrix[i][j] = 1.00;
			else identity_matrix[i][j] = 0.00;
		}
	}
	
	DP ftol = 1.00e-6;
	int iter;
	DP fret;

	NR::powell(starting_point, identity_matrix, ftol, iter, fret, minimize_target_extended_kalman_parameters_1_dim);
	
	cout<<"Ran succesfully in "<<call_counter<<" iterations with return value "<<fret<<endl;

	double omega = starting_point[0];
	double theta = starting_point[1];
	double xi = starting_point[2];
	double rho = starting_point[3];

	cout<<"Parameters are "<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho<<endl<<endl;

	cout<<"Press any key to continue"<<endl;
	cin.get();
	return 0;
}

DP minimize_target_extended_kalman_parameters_1_dim(Vec_I_DP & input) {
	double omega = input[0];
	double theta = input[1];
	double xi = input[2];
	double rho = input[3];

	estimate_extended_kalman_parameters_1_dim(log_stock_prices, 
		muS, 
		n_stock_prices, 
		omega, 
		theta, 
		xi, 
		rho, 
		u, 
		v, 
		estimates);

	double sum = 0;
	for(int i1 = 0; i1 < n_stock_prices; i1++)
		sum+=(log(v[i1])+u[i1]*u[i1]/v[i1]);

	//Printing out a status message to see whats going on
	if(call_counter++ % 100 == 0)
		cout<<" call_counter = "<<call_counter 
			<<" omega = "<<omega
			<<" theta = "<<theta
			<<" xi = "<<xi
			<<" rho = "<<rho
			<<" sum = "<<sum
			<<endl;

	return sum;
}


DP minimize_target_x_squre(Vec_I_DP & input) {
	DP in = input[0];
	DP out = in * in;
	cout<<"Input = "<<in<<" Returning "<<out<<endl;
	return out;
}

DP minimize_target_w_x_y_z_square(Vec_I_DP & input) {

	DP x = input[0];
	DP y = input[1];
	DP z = input[2];
	DP w = input[4];
	DP out = x*x+y*y+z*z;
	cout<<"Input = ["<<w<<","<<x<<","<<y<<","<<z<<"] Returning "<<out<<endl;
	return out;
}

void read_lines(char * fname, double * out) {
	std::ifstream fhandle(fname);
	char line[1000];
	bool first = true;
	int index = 0;
	do {
		fhandle.getline(line, 1000);
		if(first) {
			cout<<"Skipping the first line"<<endl;
			first = false;
			continue;
		}
		if(!fhandle.eof())
			out[index++] = atof(line);
	} while(!fhandle.eof());

}