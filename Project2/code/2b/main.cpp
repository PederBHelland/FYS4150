#include <iostream>
#include "armadillo"
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"
#include <iomanip>
#include <algorithm>

using namespace std;

int main()
{
    int n = 1000, i;
    mat A = zeros<mat>(n,n);
    vec V =  ;

    for(i = 1; i< n-1; i++){
        A(i,i-1) = -1./power(hbar,2);
        A(i,i) = 2./power(hbar,2) + V(i);
        A(i,i+1) = -1./power(hbar,2);
    }
    A(0,0) = A(n-1,n-1) = 2./power(hbar,2) + V(i);
    A(0,1) = A(n-1,n-2) = -1./power(hbar,2);

    cout << A;

}
