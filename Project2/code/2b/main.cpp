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
    int N = 5, i, j, k, l;
    double tau;
    double theta;
    mat A = zeros<mat>(N,N);
    mat B = zeros<mat>(N,N);
    vec rho(N), V(N);

    double t, c, s; //tan, cos, sin

    double eps = pow(10, -8); //we want all the elements to be lower than this
    double rho_min = 0;
    double rho_max = 5; //run program with different values
    double h = (rho_max-rho_min)/N;

    for(i = 1; i<N+1; i++){
        rho(i-1) = rho_min + i*h;
        V(i-1) = pow(rho(i-1),2);
    }

    for(i = 1; i< N-1; i++){
        A(i,i-1) = -1./pow(h,2);
        A(i,i) = 2./pow(h,2) + V(i);
        A(i,i+1) = -1./pow(h,2);
    }

    A(0,0) = 2./pow(h,2) + V(0);
    A(N-1,N-1) = 2./pow(h,2) + V(N-1);
    A(0,1) = A(N-1,N-2) = -1./pow(h,2);


    for(k = 0; k < N; k++){
        for(l = 0; l < N; l++){
            if( l > k){ //check all the elements over the diagonal
                if(fabs(A(k, l)) > eps){
                    //cout << fabs(A(k,l));
                    tau = (A(l,l) - A(k, k))/(2.*A(k, l));
                    //cout << tau << "   ";
                    if(tau > 0){ //choose smallest absolute value of t
                        t = -tau + sqrt(1+pow(tau,2));
                    }
                    else {
                        t = -tau - sqrt(1+pow(tau,2));
                       // cout << t;
                    }
                    c = 1./(sqrt(1+pow(t,2)));
                    s = t*c;

                    //Change the elements
                    B(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
                    B(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + A(k,k)*pow(s,2);
                    B(k,l) = (A(k,k)-A(l,l))*c*s + A(k,l)*(pow(c,2)-pow(s,2));
                    for(i = 0; i < N; i++){
                        if(i != k || i != l){
                            B(i,i) = A(i,i);
                            B(i,k) = A(i,k)*c - A(i,l)*s;
                            B(i,l) = A(i,l)*c + A(i,k)*s;
                        }
                    }



                    //cout << k << l <<" " << A(l,l) << " " << A(k,k) << " " << A(k,l) << " " << tau << " "<< theta << " \n \n";


                }
            }
        }
    }
cout << A << "\n \n \n" << B;

    /*
    double tau = ;
    double t = -tau +- sqrt(1+power(tau,2));
    double c = 1./sqrt(1+power(t, 2));
    double s = t*c;*/

}
