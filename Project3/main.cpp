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

void Euler(double x_0, double y_0, double M0, double x_1, double y_1, double M1, int N)
{
    //Euler method to calculate the movements the Earth around the sun
    // x_0, y_0 are the intial conditions for the sun
    // x_1, y_1 are the initial conditions for the Earth

    mat a1(N,2), v1(N,2), x1(N,2), F(N,2);
    vec t(N), ax(N), ay(N);
    double G = 6.67408e-11;
    double dt = 1./N;

    //Initial conditions
    x1(0,0) = x_1;
    x1(0,1) = y_1;
    double r0 = sqrt((x1(0,0)-x_0)*(x1(0,0)-x_0) + (x1(0,1)-y_0)*(x1(0,1)-y_0));
    double v0 = sqrt(G*M0/r0);
    double theta = asin(x1(0,1)/r0);
    v1(0,0) = 0.1;
    v1(0,1) = v0;
    t(0) = 0;

    for (int i = 0; i < N-1; i++){
        //Distance between the earth and the sun
        double rx = x1(i,0)-x_0;
        double ry = x1(i,1)-y_0;
        double r = sqrt(rx*rx + ry*ry);

        //Gravitational force in x and y direction
        double F_ = -(G*M0*M1)/(r*r*r); // F % r
        F(i,0) = rx*F_;
        F(i,1) = ry*F_;

        //Motion of the earth
        a1(i,0) = F(i,0)/M1;
        a1(i,1) = F(i,1)/M1;
        v1(i+1,0) = v1(i,0) + a1(i,0)*dt;
        v1(i+1,1) = v1(i,1) + a1(i,1)*dt;
        x1(i+1,0) = x1(i,0) + v1(i+1,0)*dt;
        x1(i+1,1) = x1(i,1) + v1(i+1,1)*dt;

        t(i+1) = t(i) + dt;
        cout << x1(i+1,0) << endl;
    }
    //cout << x1 << endl;
}

int main()
{
    int N = 2000;
    double M_earth = 6e24; //kg
    double M_sun = 2e30;
    double x_0 = 0, y_0 = 0, x_1 = 1.5e11, y_1 = 0; //m
    Euler(x_0, y_0, M_sun, x_1, y_1, M_earth, N);
}
