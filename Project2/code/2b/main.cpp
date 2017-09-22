#include <iostream>
#include "armadillo"
#include <math.h>
#include <fstream>
#include <chrono>
#include "time.h"
#include <iomanip>
#include <algorithm>


#define CATCH_CONFIG_RUNNER
#include "../catch.hpp"

using namespace std;
using namespace arma;

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}
TEST_CASE( "Factorials are computed", "[factorial]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}


int main()
{
    int dumb = 0;
    char* dumb1[1];
    int result = Catch::Session().run(dumb, dumb1);

    int N = 4, i, j, k, l, x, counter_zero_elements, u, v;
    int M = 100; //number of times we do the algorithm
    int number_of_elements_over_diagonal;
    double tau;
    double theta;
    mat A = zeros<mat>(N,N);
    mat B = zeros<mat>(N,N);
    mat A_original = zeros<mat>(N,N);
    vec rho(N), V(N);

    double t, c, s; //tan, cos, sin

    double eps = pow(10, -7); //we want all the elements to be lower than this
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
    A_original = A;

    B=A;
    x = 0;

    //for(x = 0; x < M; x++){
    while(1){
        for(k = 0; k < N; k++){
            for(l = k+1; l < N; l++){
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
                    B(l,k) = B(k,l);

                    for(i = 0; i < N; i++){
                        if(i != k && i != l){
                            B(i,i) = A(i,i);
                            B(i,k) = A(i,k)*c - A(i,l)*s;
                            B(i,l) = A(i,l)*c + A(i,k)*s;
                            B(k,i) = A(i,k)*c - A(i,l)*s;
                            B(l,i) = A(i,l)*c + A(i,k)*s;
                        }
                    }

                    A=B;
                    //cout << k << l <<" " << A(l,l) << " " << A(k,k) << " " << A(k,l) << " " << tau << " "<< theta << " \n \n";
                }
            }
        }
        //Check if elements are zero
        counter_zero_elements = 0;
        number_of_elements_over_diagonal = (N*N - N)/2;
        for(u = 0; u < N; u++){
            for(v = u+1; v < N; v++){
                if(fabs(A(u, v)) < eps){
                    //cout << u << " " << v << "\n \n \n";
                    counter_zero_elements++;
                }
                //cout << number_of_elements_over_diagonal << "   " << counter_zero_elements << " \n\n";
                if(counter_zero_elements == number_of_elements_over_diagonal){
                    cout << "hei \n \n";
                    goto theEnd;
                }
            }
        }
    }

    theEnd:
        cout << A_original << " \n" << B;
        cout << "";
}
