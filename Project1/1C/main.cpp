#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"
#include <iomanip>      // std::setprecision

using namespace std;
using namespace arma;

int main()
{
    int n = 10000000, i;
    double h = 1./(n+1);
    vec x(n), temp(n), f(n), b_(n), u(n), v(n);

    for(i=1 ; i <= n ; i++) {
        x[i] = i*h;
        f[i] = 100*exp(-10*x[i]);
        u[i] = 1 - (1-exp(-10))*x[i]-exp(-10*x[i]);
        b_[i]= pow(h,2) * f[i];
    }

    double btemp = 2;
    v[1] = b_[1]/btemp;

    clock_t t;
    t = clock();
    for(i=2 ; i <= n-1 ; i++) {
        temp[i] = -1./btemp;
        btemp = 2 + temp[i];
        v[i] = (b_[i] + v[i-1])/btemp;
    }

    for(i=n-1 ; i >= 1 ; i--){
        v[i] -= temp[i+1]*v[i+1];
    }

    t = clock()-t;
    double sec = ((double)t)/CLOCKS_PER_SEC;
    cout << setprecision(30) << sec << endl;

    //Print to output file
    const char *path="/uio/hume/student-u69/pederbh/FYS4150/Project1/results1c.txt";
    ofstream myfile(path);
    myfile << n << "\n";
    for(i=0; i<n; i++){
        myfile << x[i] << " " << u[i] << " " << v[i] << "\n";
    }
    myfile.close();
}
