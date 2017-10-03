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
}

int Jacobi_rotation(double eps, int N, mat& A, mat& R)
{
    double t, tau, c, s, r_ik, r_il;
    int number_of_rotations = 0;
    while(1){
        number_of_rotations++;
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
                    //All elements are zero
                    goto theEnd;
                }
            }
        }
}
theEnd:
    return number_of_rotations;
}

vec eigenvalues(int N, mat& A, mat& R, vec& eig_sorted, uvec& eigvec_sorted)
{
    vec eig(N);
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            if(i == j){
                eig(i) = A(i,j);
            }
        }
    }
    eig_sorted = sort(eig);
    eigvec_sorted = sort_index(eig);
    return eig_sorted;
}

void eigenvalues_matrix(int N, mat& A, mat& R, string potential, double eps, vec& eig_sorted, uvec& eigvec_sorted, vec& rho, vec w_r){
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

void print_eigenvectors_to_file_without_repulsion(int N, vec& rho, mat& R, uvec& eigvec_sorted){
    string path = string("/uio/hume/student-u69/pederbh/FYS4150/Project2/results_without_repulsion.txt");
    ofstream myfile(path);
    myfile << N << " " << N << endl;

    for(int j=0; j<3; j++){
        for(int i=0; i<N; i++){
            myfile << rho(i) << " " << R(i, eigvec_sorted(j)) << "\n";
        }
    }

    myfile.close();
}

void print_eigenvectors_to_file(int N, vec& rho, mat& R, uvec& eigvec_sorted, double w_r){
    //uio/hume/student-u85/monande/FYS4150/
    string path = string("/uio/hume/student-u69/pederbh/FYS4150/Project2/results_w") + to_string(w_r) + ".txt";
    ofstream myfile(path);
    cout << w_r << endl;
    myfile << N << " " << w_r <<endl;


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

    int N = 400;
    mat A = zeros<mat>(N,N);
    mat R = eye(N,N);
    int rho_max_ = 5;
    vec eig_sorted = zeros<vec>(N);
    vec rho = zeros<vec>(N);
    uvec eigvec_sorted = zeros<uvec>(N);
    mat eigvec_arma;
    vec eig_arma;
    double eps  = 1e-8;
    int w_r_ = 0;
    setup(N, rho_max_, A, R, "harmonic", rho, w_r_);
    Jacobi_rotation(eps, N, A, R);
    print_eigenvectors_to_file_without_repulsion(N, rho, R, eigvec_sorted);

    vec w_r = {0.01, 0.5, 1.0, 5.0}; // For comparison with the article: w_r = 0.25, 0.05
    vec rho_max = {55, 8, 5.5, 2.5}; // 10, 40



    for(int i = 0; i < w_r.size(); i++){  

        //Jacobi rotation
        setup(N, rho_max(i), A, R, "harmonic", rho, w_r(i));
        clock_t t;
        t = clock();
        int number_of_rotations = Jacobi_rotation(eps, N, A, R);       
        vec eig = eigenvalues(N, A, R, eig_sorted, eigvec_sorted);
        t = clock()-t;
        double sec = ((double)t)/CLOCKS_PER_SEC;
        cout << "Number of rotations: " << number_of_rotations << endl;
        cout << "Time Jacobi rotation: " << setprecision(30) << sec << endl;

        print_eigenvectors_to_file(N, rho, R, eigvec_sorted, w_r(i));
        //cout << "The calculated eigenvalues using Jacobi's rotation method: " << endl;
        //cout << eig << endl;

        //Armadillo
        setup(N, rho_max(i), A, R, "harmonic", rho, w_r(i));
        clock_t t2;
        t2 = clock();
        eig_sym(eig_arma, eigvec_arma, A);
        t2 = clock()-t2;
        double sec2 = ((double)t2)/CLOCKS_PER_SEC;
        cout << "Time Armadillo: " << setprecision(30) << sec2 << endl;

        //cout << "The calculated eigenvalues using Armadillo's eigsym function: " << endl;
        //cout << eig_arma << endl;

        //Reset values
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
