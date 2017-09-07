#include <iostream>
#include "armadillo"
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"

using namespace std;
using namespace arma;

int main()
{
    //double start, finish;
    clock_t t;
    t = clock();
    int n = 1000, i;
    int flops;
    double h = 1./(n+1);
    vec x(n), temp(n), f(n), b_(n), u(n), v(n);
    vec a, b, c;
    a = zeros<vec>(n) - 1;
    b = zeros<vec>(n) +2;
    c = zeros<vec>(n) - 1;


    for(i=1 ; i <= n ; i++) {
        x[i] = i*h;
        f[i] = 100*exp(-10*x[i]);
        u[i] = 1 - (1-exp(-10))*x[i]-exp(-10*x[i]);
        b_[i]= pow(h,2) * f[i];
    }

    double btemp = b[1];
    v[1] = b_[1]/btemp;
    int counter = 0;

    for(i=2 ; i <= n-1 ; i++) {
        counter++;
        temp[i] = c[i-1]/btemp;
        btemp = b[i] - a[i]*temp[i];
        v[i] = (b_[i] - a[i]*v[i-1])/btemp;
    }
    //cout << counter;

    int counter2=0;
    for(i=n-1 ; i >= 1 ; i--){
        counter2++;
        v[i] -= temp[i+1]*v[i+1];
    }

    //cout << clock();
    //cout <<counter2;
    //cout << u << v;

    //FLOPS
    flops = 8*n - 14;
    //cout << flops;

    //Print to output file
    const char *path="/uio/hume/student-u69/pederbh/FYS4150/Project1/results.txt";
    ofstream myfile(path);
    myfile << n << "\n";
    for(i=0; i<n; i++){
        myfile << x[i] << " " << u[i] << " " << v[i] << "\n";
    }
    myfile.close();

    t = clock()-t;
    double sec = ((float)t)/CLOCKS_PER_SEC;
    cout << sec;
}
