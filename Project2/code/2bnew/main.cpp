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

void setup(int N, double rho_max, mat& A, mat& R, string potential, vec& rho, double w_r=0) {
    double rho_min = 0;
    double h = (rho_max-rho_min)/(N+1);

    vec V(N), eig(N), eig_vec(N);

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
    //B=A;
}

void Jacobi_rotation(double eps, int N, mat& A, mat& R)
{
    double t, tau, c, s, r_ik, r_il;
    while(1){
        for(int k = 0; k < N; k++){
            for(int l = k+1; l < N; l++){
                if(fabs(A(k, l)) > eps){
                    tau = (A(l,l) - A(k, k))/(2.*A(k, l));
                    if(tau > 0){ //choose smallest absolute value of t
                        t = 1.0/(tau+sqrt(1.0+tau*tau));
                    }
                    else {
                        t = -1.0/(-tau+sqrt(1.0+tau*tau));
                    }
                    c = 1./(sqrt(1+(t*t)));
                    s = t*c;

                    //Change the elements
                    double A_kk = A(k,k);
                    double A_ll = A(l,l);
                    A(k,k) = A(k,k)*(c*c) - 2*A(k,l)*c*s + A(l,l)*(s*s);
                    A(l,l) = A(l,l)*(c*c) + 2*A(k,l)*c*s + A_kk*(s*s);
                    A(k,l) = (A_kk-A_ll)*c*s + A(k,l)*((c*c)-(s*s));
                    A(l,k) = A(k,l);

                    for(int i = 0; i < N; i++){
                        if(i != k && i != l){
                            double A_ik = A(i,k);
                            A(i,i) = A(i,i);
                            A(i,k) = A_ik*c - A(i,l)*s;
                            A(i,l) = A(i,l)*c + A_ik*s;
                            A(k,i) = A(i,k);
                            A(l,i) = A(i,l);
                        }
                        r_ik = R(i,k);
                        r_il = R(i,l);
                        R(i,k) = c*r_ik - s*r_il;
                        R(i,l) = c*r_il + s*r_ik;
                    }
                }
            }
        }

        //Check if elements are zero
        int counter_zero_elements = 0;
        int number_of_elements_over_diagonal = (N*N - N)/2;
        for(int u = 0; u < N; u++){
            for(int v = u+1; v < N; v++){
                if(fabs(A(u, v)) < eps){
                    counter_zero_elements++;
                }
                if(counter_zero_elements == number_of_elements_over_diagonal){
                    goto theEnd;
                }
            }
        }
}
theEnd:
    cout << "";
}

vec eigenvalues(int N, mat& A, mat& R, vec& eig_sorted, uvec& eigvec_sorted)
{
    //cout << A << endl;
    vec eig(N);
    //cout << A << endl;
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            if(i == j){
                eig(i) = A(i,j);
            }
        }
    }
    //cout << eig << endl;
    eig_sorted = sort(eig);
    //cout << eig << endl;
    eigvec_sorted = sort_index(eig);
    //cout << eig_sorted << endl;
    return eig_sorted;
}

void eigenvalues_matrix(int N, mat& A, mat& R, string potential, double eps, vec& eig_sorted, uvec& eigvec_sorted, vec& rho, vec w_r){
    //Oppgave d
    int n = w_r.size();
    mat eig_ = zeros<mat>(N,n);
    for (int i = 0; i<n; i++){
        setup(N, 5, A, R, potential, rho, w_r(i));
        Jacobi_rotation(eps, N, A, R);
        vec eig = eigenvalues(N, A, R, eig_sorted, eigvec_sorted);
        for (int j=0; j<N; j++){
            eig_(j,i) = eig(j);
        }
    }
}

void print_eigenvectors_to_file(int N, vec& rho, mat& R, uvec& eigvec_sorted, double w_r){
    //Print to output file
    string path = string("/uio/hume/student-u85/monande/FYS4150/Project2/results_w") + to_string(w_r) + ".txt";
    ofstream myfile(path);
    myfile << N << endl;

    for(int j=0; j<3; j++){
        for(int i=0; i<N; i++){
            myfile << rho(i) << " " << R(i, eigvec_sorted(j)) << "\n";
        }
    }

    myfile.close();
}

int main(int argc, char* argv[]) {

    int dumb = 0;
    char* dumb1[1];
    int result = Catch::Session().run(dumb, dumb1);

    int N = 200;
    mat A = zeros<mat>(N,N);
    mat R = eye(N,N);

    vec w_r = {0.01, 0.5, 1.0, 5.0};
    vec eig_sorted = zeros<vec>(N);
    vec rho =zeros<vec>(N);
    uvec eigvec_sorted = zeros<uvec>(N);

    int rho_max = 50;
    double eps  = 1e-8;

    //setup(N, rho_max, A, R, "harmonic", rho);
    //Jacobi_rotation(eps, N, A, R);
    //cout << eigenvalues(N, A, R, eig_sorted, eigvec_sorted);

    //eigenvalues_matrix(N, A, R, "coulomb", eps, eig_sorted, eigvec_sorted, rho, w_r);

    for(int i = 0; i < w_r.size(); i++){
        setup(N, rho_max, A, R, "coulomb", rho, w_r(i));
        Jacobi_rotation(eps, N, A, R);
        eigenvalues(N, A, R, eig_sorted, eigvec_sorted);
        print_eigenvectors_to_file(N, rho, R, eigvec_sorted, w_r(i));
        R = eye(N,N);
        A = zeros<mat>(N,N);
        eig_sorted = zeros<vec>(N);
        rho =zeros<vec>(N);
        eigvec_sorted = zeros<uvec>(N);
    }


}

TEST_CASE( "Eigenvalues are correct", "[eigenvalues]" ) {
    int N = 4;
    int rho_max = 5;
    double eps  = 1e-8;
    mat A = zeros<mat>(N,N);
    mat R = eye(N,N);
    vec rho =zeros<vec>(N);
    vec eig_sorted = zeros<vec>(N);
    uvec eigvec_sorted = zeros<uvec>(N);

    setup(N, rho_max, A, R, "harmonic", rho);
    Jacobi_rotation(eps, N, A, R);

    REQUIRE(fabs(eigenvalues(N,A,R, eig_sorted,eigvec_sorted)(0) - 2.6867) < 1e-4);
    REQUIRE(fabs(eigenvalues(N,A,R, eig_sorted,eigvec_sorted)(1) - 6.1130) < 1e-4);
    REQUIRE(fabs(eigenvalues(N,A,R, eig_sorted,eigvec_sorted)(2) - 11.0586) < 1e-4);
    REQUIRE(fabs(eigenvalues(N,A,R, eig_sorted,eigvec_sorted)(3) - 18.1417) < 1e-4);
}

TEST_CASE( "The eigenvectors are orthonormal", "[eigenvectors]" ) {
    int N = 100;
    int rho_max = 5;
    double eps  = 1e-8;
    mat A = zeros<mat>(N,N);
    mat R = eye(N,N);
    vec rho =zeros<vec>(N);
    vec eig_sorted = zeros<vec>(N);
    uvec eigvec_sorted = zeros<uvec>(N);

    setup(N, rho_max, A, R, "harmonic", rho);
    Jacobi_rotation(eps, N, A, R);
    mat R_T = trans(R);
    mat I = eye(N,N);
    mat r =( R*R_T - I);
    bool s = all(vectorise(r) < 1e-10);
    REQUIRE(s);

}


// int main og test eigenvalue skal staa her

/*
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
*/
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


    /*
    //Oppgave d
    vec w_r = {0.01, 0.5, 1., 2.};  // can change The last element should be 5 not 2
    int N = 4;                      // can change
    int n = w_r.size();
    mat eig_ = zeros<mat>(n,n);
    string potential = "coulomb";
    for (int i = 0; i<n; i++){
        vec eig = setup(N, 100, 5, 1e-7, potential, w_r(i));
        for (int j=0; j<N; j++){
            eig_(j,i) = eig(j);
    }}
    cout << eig_;
    */
    //cout << eig_;

    //if(potential == "coulomb"){
    //    cout << "The eigenvalues in the interaction case are:" << endl << eig << endl;
    //}

    //Legger til den andre main


    /*mat A;
    mat R;
    //cout << R;
    setup(10,5,A,R);
    int rotations = rotate(1e-8, A, R);
    vec lowestMeigenvalues = lowestEigenvaluesEigenvectors(3,A,R);
    cout << lowestMeigenvalues << endl;
*/

/*
TEST_CASE( "Eigenvalues are correct", "[eigenvalues]" ) {
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(0) - 2.7561) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(1) - 7.5639) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(2) - 15.3575) < 1e-4);
    REQUIRE(fabs(setup(4, 100, 5, 1e-7, "harmonic")(3) - 26.3174) < 1e-4);
}
*/
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
