#ifndef filter_utils_h
#define filter_utils_h

#include <string>
#include <vector>

using namespace std;

void read_lines(string& fname, vector<double>& out);

int sqrt_matrix(double** pa, double** proda, int& N);

/**
Checks to see if the data is normal using the chi squared test with 20 degrees of freedom and alpha = 0.05
*/
bool is_normal(double* data, int data_count);

double cumulative_normal(const double x);

double gaussrand();

void normalize_data(const double* data_in, double* data_out, const int length);

#endif