#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include "time.h"
#include "armadillo"

using namespace std;
using namespace arma;

solver::solver()
{
    total_planets = 0;
    radius = 100;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
}

solver::solver(double radi)
{
    total_planets = 0;
    radius = radi;
    total_mass = 0;
    G = 4*M_PI*M_PI;
    totalKinetic = 0;
    totalPotential = 0;
}

void solver::add(planet newplanet)
{
    total_planets += 1;
    total_mass += newplanet.mass;
    all_planets.push_back(newplanet);
}

void solver::addM(planet newplanet)
{
    total_planets +=1;
    all_planets.push_back(newplanet);
}

void solver::GravitationalConstant()
{
    //G = (4*M_PI*M_PI/32)*radius*radius*radius/total_mass;
    G = 4*M_PI*M_PI; //(AU)³/(solarMass*yr²)
}

void solver::print_position(std::ofstream &output, int dimension, double time,int number)
{   // Writes mass, position and velocity to a file "output"
    if(dimension > 3 || dimension <= 0) dimension = 3;
    else{
        for(int i=0;i<number;i++){
            planet &current = all_planets[i];
            output << time << "\t" << i+1 << "\t" << current.mass;
            for(int j=0;j<dimension;j++) output << "\t" << current.position[j];
            for(int j=0;j<dimension;j++) output << "\t" << current.velocity[j];
            output << std::endl;
        }
    }
}





void solver::print_energy(std::ofstream &output, double time,double epsilon)
{   // Writes energies to a file "output"

    this->KineticEnergySystem();
    this->PotentialEnergySystem(epsilon);
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        output << time << "\t" << nr << "\t";
        output << Current.kinetic << "\t" << Current.potential << std::endl;
    }
}

void::solver::Euler(int dimension, int integrationPoints, int finalTime, double epsilon)
{ /*Euler method to solve two coupled ODEs.*/

    // Define time step
    double timeStep = finalTime/(double)integrationPoints;
    double time = 0.0;
    double loss = 0.; // Possible energy loss
    double gravitationalConstant = 39.4784; // (AU)³/solarMass*yr²

    // Set up arrays
    double **acceleration = setup_matrix(total_planets,3);
    double **acceleration_new = setup_matrix(total_planets,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension
    int n = 0; //number of iterations

    if(doFileWriting == true){
        writeInformationToFile("Euler", integrationPoints, dimension);
    }

    // Set up clock to measure the time usage
    double startClock = clock();
    while(time < finalTime){
        // Loop over all planets
        for(int nr1=0; nr1<total_planets; nr1++){
            planet &current = all_planets[nr1]; // Current planet we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Calculate forces in each dimension
            for(int nr2=0; nr2<total_planets; nr2++){
                if(nr1!=nr2) {
                    planet &other = all_planets[nr2];
                    GravitationalForce(current,other,Fx,Fy,Fz,epsilon);
                }
            }

            // Acceleration in each dimension for current planet
            acceleration[nr1][0] = Fx/current.mass;
            acceleration[nr1][1] = Fy/current.mass;
            acceleration[nr1][2] = Fz/current.mass;

            // Calculate new position for current planet
            for(int j=0; j<dimension; j++) {
                double new_position = current.position[j] + current.velocity[j]*timeStep;
                current.position[j] = new_position;
            }

            // Calculate new velocity for current planet
            for(int j=0; j<dimension; j++){
                double new_velocity = current.velocity[j] + acceleration[nr1][j]*timeStep;
                current.velocity[j] = new_velocity;
            }

            double kineticEnergy = current.KineticEnergy();
            double potentialEnergy = 0;
            double angularMomentum = 0;
            for(int nr2=nr1+1; nr2<total_planets; nr2++){
                planet &other = all_planets[nr2];
                potentialEnergy += current.PotentialEnergy(other, gravitationalConstant, epsilon);
                angularMomentum += current.AngularMomentum(other);
            }

            if(n%freq == 0 && doFileWriting == true){
                writeToFile("Euler", current, time+timeStep, kineticEnergy, potentialEnergy, angularMomentum);
            }
        }

    time += timeStep;
    n++;
    }
    double elapsedTime = clock() - startClock;
    std::cout << "Total time Euler = " << "\t" << ((float)(elapsedTime)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
}


void solver::VelocityVerlet(int dimension, int integrationPoints, int finalTime, double epsilon)//(int dimension, int integrationPoints, double final_time, double epsilon)
{
    /*Euler method to solve two coupled ODEs.*/

    // Define time step
    double timeStep = finalTime/(double)integrationPoints;
    double time = 0.0;
    double loss = 0.; // Possible energy loss
    double gravitationalConstant = 39.4784; // (AU)³/solarMass*yr²

    // Set up arrays
    double **acceleration = setup_matrix(total_planets,3);
    double **acceleration_new = setup_matrix(total_planets,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension
    int n = 0; //number of iterations

    if(doFileWriting == true){
        writeInformationToFile("Euler", integrationPoints, dimension);
    }

    // Set up clock to measure the time usage
    double startClock = clock();
    while(time < finalTime){
        // Loop over all planets
        for(int nr1=0; nr1<total_planets; nr1++){
            planet &current = all_planets[nr1]; // Current planet we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Calculate forces in each dimension
            for(int nr2=0; nr2<total_planets; nr2++){
                if(nr1!=nr2) {
                    planet &other = all_planets[nr2];
                    GravitationalForce(current,other,Fx,Fy,Fz,epsilon);
                }
            }

            // Acceleration in each dimension for current planet
            acceleration[nr1][0] = Fx/current.mass;
            acceleration[nr1][1] = Fy/current.mass;
            acceleration[nr1][2] = Fz/current.mass;

            // Calculate new position for current planet
            for(int j=0; j<dimension; j++) {
                double new_position = current.position[j] + current.velocity[j]*timeStep;
                current.position[j] = new_position;
            }

            // Calculate new velocity for current planet
            for(int j=0; j<dimension; j++){
                double new_velocity = current.velocity[j] + acceleration[nr1][j]*timeStep;
                current.velocity[j] = new_velocity;
            }

            double kineticEnergy = current.KineticEnergy();
            double potentialEnergy = 0;
            double angularMomentum = 0;
            for(int nr2=nr1+1; nr2<total_planets; nr2++){
                planet &other = all_planets[nr2];
                potentialEnergy += current.PotentialEnergy(other, gravitationalConstant, epsilon);
                angularMomentum += current.AngularMomentum(other);
            }

            if(n%freq == 0 && doFileWriting == true){
                writeToFile("Euler", current, time+timeStep, kineticEnergy, potentialEnergy, angularMomentum);
            }
        }

    time += timeStep;
    n++;
    }
    double elapsedTime = clock() - startClock;
    std::cout << "Total time Euler = " << "\t" << ((float)(elapsedTime)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time

}


double ** solver::setup_matrix(int height,int width)
{   // Function to set up a 2D array

    // Set up matrix
    double **matrix;
    matrix = new double*[height];

    // Allocate memory
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];

    // Set values to zero
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

void solver::delete_matrix(double **matrix)
{   // Function to deallocate memory of a 2D array

    for (int i=0; i<total_planets; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void solver::GravitationalForce(planet &current,planet &other,double &Fx,double &Fy,double &Fz,double epsilon)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    double smoothing = epsilon*epsilon*epsilon;

    // Calculate the forces in each direction
    Fx -= this->G*current.mass*other.mass*relative_distance[0]/((r*r*r) + smoothing);
    Fy -= this->G*current.mass*other.mass*relative_distance[1]/((r*r*r) + smoothing);
    Fz -= this->G*current.mass*other.mass*relative_distance[2]/((r*r*r) + smoothing);
}

void solver::GravitationalForce_RK(double x_rel,double y_rel,double z_rel,double &Fx,double &Fy,double &Fz,double mass1,double mass2)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double r = sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel);

    // Calculate the forces in each direction
    Fx -= this->G*mass1*mass2*x_rel/(r*r*r);
    Fy -= this->G*mass1*mass2*y_rel/(r*r*r);
    Fz -= this->G*mass1*mass2*z_rel/(r*r*r);
}

void solver::KineticEnergySystem()
{
    totalKinetic = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.kinetic = Current.KineticEnergy();
    }
}

void solver::PotentialEnergySystem(double epsilon)
{
    totalPotential = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        Current.potential = 0;
    }
    for(int nr1=0;nr1<total_planets;nr1++){
        planet &Current = all_planets[nr1];
        for(int nr2=nr1+1;nr2<total_planets;nr2++){
            planet &Other = all_planets[nr2];
            Current.potential += Current.PotentialEnergy(Other,G,epsilon);
            Other.potential += Other.PotentialEnergy(Current,G,epsilon);
        }
    }
}

bool solver::Bound(planet OnePlanet)
{
    return ((OnePlanet.kinetic + OnePlanet.potential) < 0.0);
}


double solver::EnergyLoss()
{
    bool bound;
    vector<int> indices;
    double EnergyLoss = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        bound = this->Bound(Current);
        if(!bound){
            indices.push_back(nr);
        }
    }
    for(int i=0;i<indices.size();i++) EnergyLoss += all_planets[indices[i]].KineticEnergy();
    return EnergyLoss;

}

void::solver::writeInformationToFile(string type, int integrationPoints, int dim){
    string planetPositionPath= string("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetPosition") + type + "Results.txt";
    ofstream myPlanetPositionFile;
    myPlanetPositionFile.open(planetPositionPath,std::ios::app);

    string planetEnergiesPath= string("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetEnergies") + type + "Results.txt";
    ofstream myPlanetEnergiesFile(planetEnergiesPath);

    myPlanetPositionFile << integrationPoints << endl;
    myPlanetPositionFile << dim << endl;
    myPlanetPositionFile << total_planets << endl;

    myPlanetEnergiesFile << integrationPoints << endl;
    myPlanetEnergiesFile << dim << endl;
    myPlanetEnergiesFile << total_planets << endl;

    myPlanetPositionFile.close();
    myPlanetEnergiesFile.close();
}

void solver::writeToFile(string type, planet current, double time, double kineticEnergy, double potentialEnergy, double angularMomentum) {
    string planetPositionPath= string("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetPosition") + type + "Results.txt";
    ofstream myPlanetPositionFile;
    myPlanetPositionFile.open(planetPositionPath,std::ios::app);

    string planetEnergiesPath= string("/uio/hume/student-u85/monande/FYS4150/Project3/PlanetEnergies") + type + "Results.txt";
    ofstream myPlanetEnergiesFile(planetEnergiesPath);

    myPlanetPositionFile << current.position[0] << " " << current.position[1] << " " << current.position[2] << endl;
    myPlanetEnergiesFile << time << " " << kineticEnergy << " " << potentialEnergy << " " << angularMomentum << endl;

    myPlanetPositionFile.close();
    myPlanetEnergiesFile.close();

}
