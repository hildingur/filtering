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
        for (i=0;i<N;i++) {
          for (j=0;j<N;j++) {
            if (i == j) chol[i][i]=p[i];
            else chol[i][j]=(i > j ? a[i][j] : 0.0);
			proda[i][j] = chol[i][j];
          }
        }
		return 0;
}
