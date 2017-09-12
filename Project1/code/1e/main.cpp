#include <iostream>
#include "armadillo"
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace arma;

int main()
{
    int n = 1000, i;
    double h = 1./(n+1);
    vec x(n), temp(n), f(n), b_tilde(n), u(n), v(n), eps(n);
    mat A = zeros<mat>(n,n);

    for(i=0 ; i < n; i++) {
        x(i) = (i+1)*h;
        f(i) = 100*exp(-10*x(i));
        u(i) = 1 - (1-exp(-10))*x(i)-exp(-10*x(i));
        b_tilde(i)= pow(h,2) * f(i);
    }

    double btemp = 2;
    v(0) = b_tilde(0)/btemp;

    for(i=1 ; i < n ; i++) {
        temp(i) = -1./btemp;
        btemp = 2 + temp(i);
        v(i) = (b_tilde(i) + v(i-1))/btemp;
    }

    v(n-1) = v(n-1)/btemp;

    for(i=n-2 ; i > 0 ; i--){
        v(i) -= temp(i+1)*v(i+1);
    }

    // LU decomp with solve function
    for(i = 1; i< n-1; i++){
        A(i,i-1) = -1.;
        A(i,i) = 2.;
        A(i,i+1) = -1.;
    }
    A(0,0) = A(n-1,n-1) = 2;
    A(0,1) = A(n-1,n-2) = -1;

    clock_t t;
    t = clock();

    vec v_tilde = solve(A, b_tilde);
    mat L, U, P;
    lu(L, U, P, A);

    t = clock()-t;
    double sec = ((double)t)/CLOCKS_PER_SEC;
    cout << setprecision(30) << sec << endl;
    //cout << (A-L*U*P); //checking that A-L*U*P is the zero matrix

    //Print to output file
    const char *path="/uio/hume/student-u69/pederbh/FYS4150/Project1/results1e.txt";
    ofstream myfile(path);
    myfile << n << "\n";
    for(i=0; i<n; i++){
        myfile << x[i] << " " << u[i] << " " << v_tilde[i] << "\n";
    }
    myfile.close();
}
