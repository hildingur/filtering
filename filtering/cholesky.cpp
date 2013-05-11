#include <iostream>
#include <iomanip>
#include <random>

#include "recipes/nr.h"
using namespace std;

// Driver for routine cholsl
//Cholesky to conform to the book.
int sqrt_matrix(double** pa, double** proda, int& N)
{
       
        int i,j,k;
        DP sum;
        Vec_DP p(N),x(N);
        //Mat_DP a(*pa,N,N), atest(N,N), chol(N,N);
		Mat_DP a(N,N), atest(N,N), chol(N,N);
		
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				a[i][j] = pa[i][j];
			}
		}

		NR::choldc(a,p);
        cout << "Cholesky:" << endl;
        cout << scientific << setprecision(6);
        for (i=0;i<N;i++) {
          for (j=0;j<N;j++) {
            if (i == j) chol[i][i]=p[i];
            else chol[i][j]=(i > j ? a[i][j] : 0.0);
			proda[i][j] = chol[i][j];
            //cout << setw(16) << a_d[i*N+j];
			cout<<setw(16) << chol[i][j];
          }
          cout << endl;
        }

		return 0;
}

int main() {
	int N = 3;
	double **pa = new double*[N];
	double **proda = new double*[N];

	for(int i=0; i < N; i++) {
		pa[i] = new double[N];
		proda[i] = new double[N];
	}

	pa[0][0] = 100.00; pa[0][1] = 15.00; pa[0][1] = 00.01; 
	pa[1][0] = 15.00;  pa[1][1] = 2.3;   pa[1][1] = 00.01; 
	pa[2][0] =  0.01;  pa[2][1] = 0.01;  pa[1][1] = 1.0; 
	/*
	double pa[3][3] = {
		{100.00, 15.00, 00.01},
		{15.00, 2.3, 0.01},
		{0.01, 0.01, 1.0}
	};
	double proda[3][3] = {
		{100.00, 15.00, 00.01},
		{15.00, 2.3, 0.01},
		{0.01, 0.01, 1.0}
	};
	*/

	sqrt_matrix(pa, proda, N);

	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			cout<<setprecision(16)<<proda[i][j];
		}
		cout<<endl;
	}
	
	cin.get();

}