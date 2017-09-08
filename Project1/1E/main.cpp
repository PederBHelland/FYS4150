#include <iostream>
#include "armadillo"
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"
#include <iomanip>      // std::setprecision
#include <algorithm>

using namespace std;
using namespace arma;

int main()
{
    int n = 100000, i;
    int flops;
    double h = 1./(n+1);
    vec x(n), temp(n), f(n), b_(n), u(n), v(n), eps(n);
    vec a, b, c;
    a = zeros<vec>(n)-1;
    b = zeros<vec>(n)+2;
    c = zeros<vec>(n)-1;


    for(i=0 ; i < n; i++) {
        x(i) = (i+1)*h;
        f(i) = 100*exp(-10*x(i));
        u(i) = 1 - (1-exp(-10))*x(i)-exp(-10*x(i));
        b_(i)= pow(h,2) * f(i);
    }

    double btemp = b(0);
    v(0) = b_(0)/btemp;

    clock_t t;
    t = clock();
    for(i=1 ; i < n ; i++) {
        temp[i] = c[i-1]/btemp;
        btemp = b[i] - a[i]*temp[i];
        v[i] = (b_[i] - a[i]*v[i-1])/btemp;
    }

    v(n-1) = v(n-1)/btemp;
    for(i=n-2 ; i >= 0 ; i--){
        v[i] -= temp[i+1]*v[i+1];
    }





/*
    double eps_max = 0;
    for(i=1; i <= n ; i++) {
        eps[i] =log10(fabs((v[i]-u[i])/u[i]));
        //cout << eps[i] << endl;
        if(fabs(eps[i])>eps_max) {
            eps_max = fabs(eps[i]);
        }
    }
*/

    //double eps_max = max_element(eps, 4);
    //cout << eps_max;

    // LU decomp with solve function
    mat A = zeros<mat>(n,n);
    for(i = 1; i< n-1; i++){
        A(i,i-1) = -1.;
        A(i,i) = 2.;
        A(i,i+1) = -1.;
    }
    A(0,0) = A(n-1,n-1) = 2;
    A(0,1) = A(n-1,n-2) = -1;
    //cout << A;
    vec v_ = solve(A, b_);
    //cout << v_ << '\n '<< v<< '\n ' << u;
    t = clock()-t;
    double sec = ((double)t)/CLOCKS_PER_SEC;
    cout << setprecision(30) << sec << endl;
    mat L, U, P;
    lu(L, U, P, A);

    //cout << L << U;
    //cout << (A-L*U*P);


/*
    //Print to output file
    const char *path="/uio/hume/student-u69/pederbh/FYS4150/Project1/results1d.txt";
    ofstream myfile(path);
    myfile << n << "\n";
    for(i=0; i<n; i++){
        myfile << x[i] << " " << u[i] << " " << v[i] << "\n";
    }
    myfile.close();
    */
}
