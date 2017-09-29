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

void setup(int N, int M, double rho_max, string potential, double w_r=0) {
    int number_of_elements_over_diagonal;
    double rho_min = 0;

    mat A = zeros<mat>(N,N);
    mat B = zeros<mat>(N,N);
    mat A_original = zeros<mat>(N,N);
    mat R = zeros<mat>(N,N);
    vec rho(N), V(N), eig_vec(N);
    //eig(N)


    double h = (rho_max-rho_min)/N;


    for(int i = 1; i<N+1; i++){
        rho(i-1) = rho_min + i*h;
        if(potential == "harmonic"){
            V(i-1) = rho(i-1)*rho(i-1);
        }
        if(potential == "coulomb"){
            V(i-1) = w_r*w_r*rho(i-1)*rho(i-1) + 1./rho(i-1);
        }
    }

    for(int i = 1; i< N-1; i++){
        A(i,i-1) = -1./(h*h);
        A(i,i) = 2./(h*h) + V(i);
        A(i,i+1) = -1./(h*h);
    }

    A(0,0) = 2./(h*h) + V(0);
    A(N-1,N-1) = 2./(h*h) + V(N-1);
    A(0,1) = A(N-1,N-2) = -1./(h*h);
    A_original = A;

    B=A;
    //eigenvalues(A, B, N, k, l, tau, t, c, s, i, counter_zero_elements, number_of_elements_over_diagonal, u, v, A_original, eps, j, R, eig);

    //return eig;
}


void Jacobi_rotation(mat &A, mat &B, int N, double &tau, mat A_original, double eps, mat &R)
{
    R = eye(N,N);
    double t, c, s;
    while(1){
        for(int k = 0; k < N; k++){
            for(int l = k+1; l < N; l++){
                if(fabs(A(k, l)) > eps){
                    //cout << fabs(A(k,l));
                    double tau = (A(l,l) - A(k, k))/(2.*A(k, l));
                    //cout << tau << "   ";
                    if(tau > 0){ //choose smallest absolute value of t
                        t = -tau + sqrt(1+(tau*tau));
                    }
                    else {
                        t = -tau - sqrt(1+(tau*tau));
                        // cout << t;
                    }
                    c = 1./(sqrt(1+(t*t)));
                    s = t*c;

                    //Change the elements
                    B(k,k) = A(k,k)*(c*c) - 2*A(k,l)*c*s + A(l,l)*(s*s);
                    B(l,l) = A(l,l)*(c*c) + 2*A(k,l)*c*s + A(k,k)*(s*s);
                    B(k,l) = (A(k,k)-A(l,l))*c*s + A(k,l)*((c*c)-(s*s));
                    B(l,k) = B(k,l);

                    for(int i = 0; i < N; i++){
                        if(i != k && i != l){
                            B(i,i) = A(i,i);
                            B(i,k) = A(i,k)*c - A(i,l)*s;
                            B(i,l) = A(i,l)*c + A(i,k)*s;
                            B(k,i) = A(i,k)*c - A(i,l)*s;
                            B(l,i) = A(i,l)*c + A(i,k)*s;
                        }
                        double r_ik = R(i,k);
                        double r_il = R(i,l);
                        R(i,k) = c*r_ik - s*r_il;
                        R(i,l) = c*r_il + s*r_ik;
                    }

                    A=B;
                    //cout << k << l <<" " << A(l,l) << " " << A(k,k) << " " << A(k,l) << " " << tau << " "<< theta << " \n \n";
                }
            }
        }
        //Check if elements are zero
        int counter_zero_elements = 0;
        int number_of_elements_over_diagonal = (N*N - N)/2;
        for(int u = 0; u < N; u++){
            for(int v = u+1; v < N; v++){
                if(fabs(A(u, v)) < eps){
                    //cout << u << " " << v << "\n \n \n";
                    counter_zero_elements++;
                }
                //cout << number_of_elements_over_diagonal << "   " << counter_zero_elements << " \n\n";
                if(counter_zero_elements == number_of_elements_over_diagonal){
                    //cout << "hei \n \n";
                    goto theEnd;
                }
            }
        }
    }

theEnd:
    //cout << A_original << " \n" << A;
    cout << "";
    //return A;

}

//vec eigenvalues(mat &A, mat &B, int N, int M, double tau, mat A_original, double eps, mat& R, double w_r)
vec eigenvalues(mat A, mat B, int N, int M, double &tau, mat &A_original, mat &R, double &eps){

    Jacobi_rotation(A, B, N, tau, A_original, eps, R);
    vec eig(N);
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            if(i == j){
                eig(i) = A(i,j);
            }
        }
    }
    return eig;
}



int main(int argc, char* argv[]) {

    int dumb = 0;
    char* dumb1[1];
    int result = Catch::Session().run(dumb, dumb1);

    //Oppgave d

    vec w_r = {0.01, 0.5, 1., 5.};  // can change
    int N = 4;                      // can change
    int M = 100;
    int n = w_r.size();
    double rho_max = 2, tau, eps;
    mat eig_ = zeros<mat>(n,n), A_original;
    string potential = "coulomb";
    mat A, B, R;

    for (int i = 0; i<n; i++){
        //setup(N, 100, 2, 1e-7, potential, w_r(i));
        setup(N, M, rho_max, potential, w_r(i));
        //eig = eigenvalues(A, B, N, M, tau, A_original, eps, R, w_r(i));
        vec eig = eigenvalues(A, B, N, M, tau, A_original, R, eps) ;
        for (int j=0; j<N; j++){
            eig_(j,i) = eig(j);
    }
        //Jacobi_rotation(A, B, N, k, l, tau, t, c, s, i, counter_zero_elements, number_of_elements_over_diagonal, u, v, A_original, eps, j, R);
        //cout << R;
    }
    cout << eig_;
    //mat A = Jacobi_rotation(A, B, N, k, l, tau, t, c, s, i, counter_zero_elements, number_of_elements_over_diagonal, u, v, A_original, eps, j, R);
    cout << R;




}/*
TEST_CASE( "Eigenvalues are correct", "[eigenvalues]" ) {
    double w_r = 0;
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic", w_r)(0) - 2.7561) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic", w_r)(1) - 7.5639) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic", w_r)(2) - 15.3575) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic", w_r)(3) - 26.3174) < 1e-4);
}*/

/*
TEST_CASE("hei") {
    mat A;
    mat R;
    setup(4,10,A,R);
    int rotations = rotate(1e-8, A, R);
    vec lowestMeigenvalues = lowestEigenvaluesEigenvectors(3,A,R);
    REQUIRE(arma::dot(R(:,0), R(:,0) == Approx(1)));
    REQUIRE(arma::dot(R(:,0), R(:,1) == Approx(0)));
    REQUIRE(lowestMeigenvalues(0) == Approx(3.938768375));
}
*/
