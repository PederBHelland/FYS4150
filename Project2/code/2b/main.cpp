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

mat Jacobi_rotation(mat A, mat B, int N, int k, int l, double tau, double t, double c, double s, int i, int counter_zero_elements, int number_of_elements_over_diagonal, int u, int v, mat A_original, double eps, int j, mat &R)
{
    R = eye(N,N);
    while(1){
        for(k = 0; k < N; k++){
            for(l = k+1; l < N; l++){
                if(fabs(A(k, l)) > eps){
                    //cout << fabs(A(k,l));
                    tau = (A(l,l) - A(k, k))/(2.*A(k, l));
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

                    for(i = 0; i < N; i++){
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
                    //cout << "hei \n \n";
                    goto theEnd;
                }
            }
        }
    }

theEnd:
    //cout << A_original << " \n" << A;
    cout << "";
    return A;

}

vec eigenvalues(mat A, mat B, int N, int k, int l, double tau, double t, double c, double s, int i, int counter_zero_elements, int number_of_elements_over_diagonal, int u, int v, mat A_original, double eps, int j, mat R)
{
    A = Jacobi_rotation(A, B, N, k, l, tau, t, c, s, i, counter_zero_elements, number_of_elements_over_diagonal, u, v, A_original, eps, j, R);
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

vec setup(int N, int M, double rho_max, double eps, string potential, double w_r=0) {
    int i, j, k, l, x, counter_zero_elements, u, v;
    int number_of_elements_over_diagonal;
    double tau;
    double theta;
    double rho_min = 0;

    mat A = zeros<mat>(N,N);
    mat B = zeros<mat>(N,N);
    mat A_original = zeros<mat>(N,N);
    mat R = zeros<mat>(N,N);
    vec rho(N), V(N), eig(N), eig_vec(N);

    double t, c, s; //tan, cos, sin

    double h = (rho_max-rho_min)/N;

    for(i = 1; i<N+1; i++){
        rho(i-1) = rho_min + i*h;
        if(potential == "harmonic"){
            V(i-1) = rho(i-1)*rho(i-1);
        }
        if(potential == "coulomb"){
            V(i-1) = w_r*w_r*rho(i-1)*rho(i-1) + 1./rho(i-1);
        }
    }

    for(i = 1; i< N-1; i++){
        A(i,i-1) = -1./(h*h);
        A(i,i) = 2./(h*h) + V(i);
        A(i,i+1) = -1./(h*h);
    }

    A(0,0) = 2./(h*h) + V(0);
    A(N-1,N-1) = 2./(h*h) + V(N-1);
    A(0,1) = A(N-1,N-2) = -1./(h*h);
    A_original = A;

    B=A;
    eig = eigenvalues(A, B, N, k, l, tau, t, c, s, i, counter_zero_elements, number_of_elements_over_diagonal, u, v, A_original, eps, j, R);

    return eig;
}


// int main og test eigenvalue skal staa her


void setup(int N, double rho_max, mat& A, mat& R) {
    // construct A, R
}

int rotate(double eps, mat& A, mat& R) {
    // rotate
    // A is diagonal
    // R is eigenvec
    int num_iterations;
    return num_iterations;
}

vec lowestEigenvaluesEigenvectors(int M, mat& A, mat& R) {
    // sort A
    // sort R accoding to A
    //bruk sort(A) og sort_index(A). Sort_index(A) gir array som svarer til posisjonene til tilsvarende egenvektorer i matrisen R
    vec lowestMeigenvalues; // vector of the M lowest eigenvalues
    return lowestMeigenvalues;
}
/*
int main() {
    mat A;
    mat R;
    setup(10,5,A,R);
    int rotations = rotate(1e-8, A, R);
    vec lowestMeigenvalues = lowestEigenvaluesEigenvectors(3,A,R);
    cout << lowestMeigenvalues << endl;
}
*/
int main(int argc, char* argv[]) {
    int dumb = 0;
    char* dumb1[1];
    int result = Catch::Session().run(dumb, dumb1);

    //Oppgave d
    double w_r = 0.01;
    string potential = "coulomb";
    vec eig = setup(3, 100, 5, 1e-7, potential, w_r);
    //if(potential == "coulomb"){
    //    cout << "The eigenvalues in the interaction case are:" << endl << eig << endl;
    //}

    //Legger til den andre main


    mat A;
    mat R;
    //cout << R;
    setup(10,5,A,R);
    int rotations = rotate(1e-8, A, R);
    vec lowestMeigenvalues = lowestEigenvaluesEigenvectors(3,A,R);
    cout << lowestMeigenvalues << endl;
    cout << R << endl;

}

TEST_CASE( "Eigenvalues are correct", "[eigenvalues]" ) {
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(0) - 2.7561) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(1) - 7.5639) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(2) - 15.3575) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(3) - 26.3174) < 1e-4);
}

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
