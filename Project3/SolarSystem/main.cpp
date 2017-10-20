#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"
#include <cstdio>

using namespace std;

void PrintInitialValues(int, double, double, double *, double *, int);
void PrintFinalValues(int, double *, double *);


int main()
{
    int IntegrationPoints;  // No. of integration points
    int FinalTime;       // End time of calculation
    int Dimension;           // No. of spatial dimensions

    remove("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetPositionEulerResults.txt");
    remove("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetEnergiesEulerResults.txt");
    remove("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetPositionVerletResults.txt");
    remove("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetEnergiesVerletResults.txt");

        cout << "Earth-Sun binary system" << endl;
        Dimension = 2;

        FinalTime = 1;
        IntegrationPoints = 10000*FinalTime;

        double epsilon = 0.0;

        double TimeStep = FinalTime/IntegrationPoints;
        double x[3],v[3];  // positions and velocities
        // initial position x = 1AU, y = z = 0, vx = 2pi, vy=0, vz=0
        //planet planet1(0.000003,1.,0.0,0.0,0.0,6.3,0.); // Earth: (mass,x,y,z,vx,vy,vz)
        planet planet1(1.,0.,0.,0.,0.,0.,0.);           // Sun: (mass,x,y,z,vx,vy,vz)
        planet planet2(0.000003,1.,0.0,0.0,0.0,2*M_PI,0.); //Earth, testing to find circular orbit
//        planet planet3(0.95e-3, 5.2, 0., 0., 0., 2*M_PI/sqrt(5.2), 0. ); //Jupiter 2.76
//        planet planet4(3.3e-7, 1.52, 0., 0., 0., 2*M_PI/sqrt(1.52), 0. ); //Mars
//        planet planet5(2.45e-6, 0.72, 0., 0., 0., 2*M_PI/sqrt(0.72), 0. ); //Venus
//        planet planet6(2.75e-4, 9.54, 0., 0., 0., 2*M_PI/sqrt(9.54), 0. ); //Saturn
//        planet planet7(1.65e-7, 0.39, 0., 0., 0., 2*M_PI/sqrt(0.39), 0. ); //Merkur
//        planet planet8(4.4e-5, 19.19, 0., 0., 0., 2*M_PI/sqrt(19.19), 0. );//Uranus
//        planet planet9(0.515e-4, 30.06, 0., 0., 0., 2*M_PI/sqrt(30.06), 0. ); //Neptun
//        planet planet10(0.655e-8, 39.53, 0., 0., 0., 2*M_PI/sqrt(39.53), 0. ); //Pluto

        solver binary_vv(5.0);
        binary_vv.setFileWriting(true,100);

        binary_vv.add(planet1);
        binary_vv.add(planet2);
//        binary_vv.add(planet3);
//        binary_vv.add(planet4);
//        binary_vv.add(planet5);
//        binary_vv.add(planet6);
//        binary_vv.add(planet7);
//        binary_vv.add(planet8);
//        binary_vv.add(planet9);
//        binary_vv.add(planet10);

        PrintInitialValues(Dimension,TimeStep,FinalTime,x,v,IntegrationPoints);

        //cout << "INTIANEOPITJNEOAI: " << planet x1 << endl;

        cout << "Velocity Verlet results for the Sun-Earth system:" << endl;
        binary_vv.VelocityVerlet(Dimension,IntegrationPoints,FinalTime,epsilon);

        //binary_vv.Euler(Dimension, IntegrationPoints, FinalTime, epsilon);

        for(int j = 0; j < Dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }
        PrintFinalValues(Dimension,x,v);
    return 0;
}



void PrintInitialValues(int Dimension,double TimeStep, double FinalTime,double *x_initial,double *v_initial, int N){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << TimeStep << "; final time = " << FinalTime << "; integration points = " << N << endl;

    cout << "Initial position = ";
    for(int j=0;j<Dimension;j++) cout << x_initial[j] << " ";
    cout << endl;

    cout << "Initial velocity = ";
    for(int j=0;j<Dimension;j++) cout << v_initial[j] << " ";
    cout << endl;
}

void PrintFinalValues(int Dimension,double *x_final,double *v_final){
    // A function that prints out the final results of the calculation

    cout << "Final position = ";
    for(int j=0; j<Dimension; j++) cout << x_final[j] << " ";
    cout << endl;

    cout << "Final velocity = ";
    for(int j=0; j<Dimension; j++) cout << v_final[j] << " ";
    cout << endl;
}
