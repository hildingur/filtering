#include "filter_utils.h"
#include "recipes/nr.h"

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

void read_lines(string& fname, vector<double>& out) {
	std::ifstream fhandle(fname.c_str());
	char line[1000];
	int index = 0;
	do {
		fhandle.getline(line, 1000);
		if(!fhandle.eof())
			out.push_back(atof(line));

	} while(!fhandle.eof());

}

// Driver for numerical recipies choldc routine
//Cholesky to conform to the book.
int sqrt_matrix(double** pa, double** proda, int& N)
{
       
        int i,j;
        Vec_DP p(N),x(N);
        //Mat_DP a(*pa,N,N), atest(N,N), chol(N,N);
		Mat_DP a(N,N), atest(N,N), chol(N,N);
		
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				a[i][j] = pa[i][j];
			}
		}

		NR::choldc(a,p);
        for (i=0;i<N;i++) {
          for (j=0;j<N;j++) {
            if (i == j) chol[i][i]=p[i];
            else chol[i][j]=(i > j ? a[i][j] : 0.0);
			proda[i][j] = chol[i][j];
          }
        }
		return 0;
}

/**
Checks to see if the data is normal using the k-s test
*/
bool is_normal(double* data, int data_count)
{
	const int bin_count = 22;
	const int edge_count = bin_count - 1;
	double bin_edges[edge_count];
	const double starting_bin = -2.5;
	const double increment = 0.25;
	
	cout.setf(ios::fixed);
	cout<<setprecision(10);

	//ebins - expected bins, bins - actual bins
	Vec_IO_DP ebins(bin_count), bins(bin_count);
	for(int i = 0; i < edge_count; i++) 
	{
		bin_edges[i] = starting_bin + increment * ((double) i);
		if(i==0)  //The first bin
		{			
			ebins[i] = ((double) data_count) * cumulative_normal(bin_edges[i]);
		}
		else
		{
			ebins[i] = ((double) data_count) * (cumulative_normal(bin_edges[i]) - cumulative_normal(bin_edges[i-1]));
		} 

		bins[i] = 0.00; //initialing to zero.
	}

	//the last ebin: ebin[bin_count - 1] =  is (1 - cumulative_normal(bin_edges[edge_count - 1])) * data_count
	ebins[bin_count - 1] = (1 - cumulative_normal(bin_edges[edge_count - 1])) * data_count;
	bins[bin_count - 1] = 0;

	//binning the data
	for(int i = 0; i < data_count; i++)
	{
		bool binned = false;
		double num = data[i];
		//iterate thru the edges looking for the appropriate bin to bin it into
		for(int j = 0; j < edge_count; j++)
		{
			if(num <= bin_edges[j]) {
				binned = true;
				bins[j]++;
				break;
			} 
		} 
		if(!binned)
		{
			bins[edge_count]++;
		}
	}

	double* normalized_data = new double[data_count];
	normalize_data(data, normalized_data, data_count);
	

	double df, chsq, prob;

	NR::chsone(bins, ebins, 2, df, chsq, prob);
	cout<<"df = "<<df<<endl
		<<"chsq = "<<chsq<<endl
		<<"prob = "<<prob<<endl;
	
	delete[] normalized_data;

	return chsq < 31.50;
}

double cumulative_normal(const double x_)
{
	
    //return 1.0 - NR::erfcc(x_/sqrt(2.0));

	double x = x_;
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

void normalize_data(const double* data_in, double* data_out, const int length)
{
	//get the moments of the data for normalization
	Vec_I_DP data_vecor(data_in, length);
	double ave, adev, sdev, var, skew, curt;
	NR::moment(data_vecor, ave, adev, sdev, var, skew, curt);
	for(int i = 0; i < length; i++) {
		data_out[i] = (data_in[i] - ave) / sdev;
	}
}

double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}
