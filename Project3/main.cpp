#include <iostream>
#include "armadillo"
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"
#include <iomanip>
#include <algorithm>

using namespace std;
//using namespace arma;


void Euler()
{
    int N = 100;
    arma::vec a, v, x, t;
    double dt = 1/N;
    double M_earth = 6e24;
    double M_sun = 2e30;
    double G = 6.67408e-11;
    double r = 1.5e11;
//    double m = M_earth;
    double F = (G*M_sun*M_earth)/(r*r);
    for (int i = 0; i < N-1; i++){
        a(i) = F/M_earth;
        v(i+1) = v(i) + a(i+1)*dt;
        x(i+1) = x(i)  + v(i+1)*dt;
        t(i+1) = t(i) + dt;
    }
}

int main()
{
//Euler();
}
