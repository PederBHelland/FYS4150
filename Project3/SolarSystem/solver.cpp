#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include <iomanip>
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

            double sunPosition_x;
            double sunPosition_y;
            double sunPosition_z;
            if (current.mass == 1.){
               sunPosition_x = current.position[0];
               sunPosition_y = current.position[1];
               sunPosition_z = current.position[2];
            }


            current.position[0] -= sunPosition_x;
            current.position[1] -= sunPosition_y;
            current.position[2] -= sunPosition_z;

            double kineticEnergy = current.KineticEnergy();
            double potentialEnergy = 0;
            double angularMomentum = 0;
            for(int nr2=0; nr2<total_planets; nr2++){
                if(nr2 != nr1){
                    planet &other = all_planets[nr2];
                    potentialEnergy += current.PotentialEnergy(other, gravitationalConstant, epsilon);
                    angularMomentum += current.AngularMomentum(other);
                }

            }

            if(n%freq == 0 && doFileWriting == true){
                writeToFile("Euler", current, time+timeStep, integrationPoints, kineticEnergy, potentialEnergy, angularMomentum);
            }
        }

    time += timeStep;
    n++;
    }
    double elapsedTime = clock() - startClock;
    std::cout << "Total time Euler = " << "\t" << ((float)(elapsedTime)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
}

void solver::VelocityVerlet(int dimension, int integrationPoints, double final_time, double epsilon)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    // Define time step
    double timeStep = final_time/((double) integrationPoints);
    double time = 0.0;
    double loss = 0.; // Possible energy loss
    //double gravitationalConstant = 39.4784; // (AU)³/solarMass*yr²
    int lostPlanets[integrationPoints];
    double gravitationalConstant = 39.4784; // (AU)³/solarMass*yr²

    // Set up arrays
    double **acceleration = setup_matrix(total_planets,3);
    double **acceleration_new = setup_matrix(total_planets,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension

    int n = 0;
    lostPlanets[n] = 0;
    //output_lost << time << "\t" << lostPlanets[n] << std::endl;
    n+=1;

    // Set up clock to measure the time usage
    clock_t planet_VV,finish_VV;
    planet_VV = clock();

    // PLANET CALCULATIONS
    // Loop over time
    time += timeStep;

    if(doFileWriting == true){
        writeInformationToFile("Verlet", integrationPoints, dimension);
    }

    while(time < final_time){
        if ((n % 1000) == 0) cout << time/((double)final_time)*100.0 << endl;
        lostPlanets[n] = 0;
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
                current.position[j] += current.velocity[j]*timeStep + 0.5*timeStep*timeStep*acceleration[nr1][j];
            }

//            double sunPosition_x;
//            double sunPosition_y;
//            double sunPosition_z;
//            if (current.mass == 1.){
//               sunPosition_x = current.position[0];
//               sunPosition_y = current.position[1];
//               sunPosition_z = current.position[2];
//            }
//            current.position[0] -= sunPosition_x;
//            current.position[1] -= sunPosition_y;
//            current.position[2] -= sunPosition_z;


            // Update position so that the center of mass is always in the origin
            vec centerOfMassVector = centerOfMass(dimension);
            current.position[0] -= centerOfMassVector[0];
            current.position[1] -= centerOfMassVector[1];
            current.position[2] -= centerOfMassVector[2];


            if (current.mass == 1.65e-7){
                double mercuryDistance = sqrt(current.position[0]*current.position[0] + current.position[1]*current.position[1] + current.position[2]*current.position[2]);
                if(n%freq == 0 && doFileWriting == true) {
//                    writeMercuryToFile(mercuryDistance, time, integrationPoints);
                }

            }


            // Loop over all other planets
                for(int nr2=0; nr2<total_planets; nr2++){
                    if(nr1!=nr2) {
                        planet &other = all_planets[nr2];
                        GravitationalForce(current,other,Fxnew,Fynew,Fznew,epsilon);
                    }
                }

            // Acceleration each dimension exerted for current planet
            acceleration_new[nr1][0] = Fxnew/current.mass;
            acceleration_new[nr1][1] = Fynew/current.mass;
            acceleration_new[nr1][2] = Fznew/current.mass;

            // Calculate new velocity for current planet
            for(int j=0; j<dimension; j++){
                current.velocity[j] += 0.5*timeStep*(acceleration[nr1][j] + acceleration_new[nr1][j]);
            }


            double kineticEnergy = current.KineticEnergy();
            double potentialEnergy = 0;
            double angularMomentum = 0;
            for(int nr2=0; nr2<total_planets; nr2++){
                if(nr2 != nr1){
                    planet &other = all_planets[nr2];
                    potentialEnergy += current.PotentialEnergy(other, gravitationalConstant, epsilon);
                    angularMomentum += current.AngularMomentum(other);

                }
            }

            if(n%freq == 0 && doFileWriting == true){
                writeToFile("Verlet", current, time+timeStep, integrationPoints, kineticEnergy, potentialEnergy, angularMomentum);
            }
        }

//        loss += EnergyLoss();

//        for(int nr=0;nr<total_planets;nr++){
//            planet &Current = all_planets[nr];
//            if(!(this->Bound(Current))){
//                lostPlanets[n] += 1;
//            }
//        }

        n++;
        time += timeStep;
    }

    //Stop clock and print out time usage
    finish_VV = clock();
    std::cout << "Total time Verlet = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time

    std::cout << "One time step = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC)/integrationPoints << " seconds" << std::endl; // print elapsed time

    //loss = EnergyLoss();
    std::cout << "Total energyloss due to unbound planets: " << loss << std::endl;

    double boundPlanets = 0;
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        if(this->Bound(Current)){
            //output_bound << nr << std::endl;
            boundPlanets += 1;
        }
    }
    std::cout << "There are " << boundPlanets << " bound planets at the end of the run" << std::endl;

    //Clear memory
    delete_matrix(acceleration);
    delete_matrix(acceleration_new);
}

void solver::VelocityVerletMercury(int dimension, int integrationPoints, double final_time, double epsilon, planet &Sun, planet &Mercury)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    // Define time step
    double timeStep = final_time/((double) integrationPoints);
    //cout << timeStep << endl;
    double time = 0.0;
    double gravitationalConstant = 39.4784; // (AU)³/solarMass*yr²

    // Set up arrays
    vec acceleration(3);
    vec acceleration_new(3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension

    int n = 0;
    n+=1;

    // Set up clock to measure the time usage
    clock_t planet_VV,finish_VV;
    planet_VV = clock();

    // PLANET CALCULATIONS
    // Loop over time
    time += timeStep;

    while(time < final_time){
        //cout << time << endl;
        if ((n % 1000) == 0) cout << time/((double)final_time)*100.0 << endl;

        Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

        GravitationalForceSunMercury(Mercury,Sun,Fx,Fy,Fz,epsilon);

        // Acceleration in each dimension for current planet
        acceleration[0] = Fx/Mercury.mass;
        acceleration[1] = Fy/Mercury.mass;
        acceleration[2] = Fz/Mercury.mass;

        // Calculate new position for current planet
        for(int j=0; j<dimension; j++) {
            //cout << acceleration[j] << endl;
            Mercury.position[j] += Mercury.velocity[j]*timeStep + 0.5*timeStep*timeStep*acceleration[j];
        }

        // Update position so that the center of mass is always in the origin
/*        vec centerOfMassVector = centerOfMass(dimension);
        Mercury.position[0] -= centerOfMassVector[0];
        Mercury.position[1] -= centerOfMassVector[1];
        Mercury.position[2] -= centerOfMassVector[2];*/

//      double mercuryDistance = sqrt(Mercury.position[0]*Mercury.position[0] + Mercury.position[1]*Mercury.position[1] + Mercury.position[2]*Mercury.position[2]);

        //double mercuryDistance = Mercury.distance(Sun);
        GravitationalForceSunMercury(Mercury, Sun, Fxnew, Fynew, Fznew, epsilon);


        // Acceleration each dimension exerted for current planet
        acceleration_new[0] = Fxnew/Mercury.mass;
        acceleration_new[1] = Fynew/Mercury.mass;
        acceleration_new[2] = Fznew/Mercury.mass;

        // Calculate new velocity for current planet
        for(int j=0; j<dimension; j++){
            Mercury.velocity[j] += 0.5*timeStep*(acceleration[j] + acceleration_new[j]);
        }

        //REPEAT FOR THE SUN

        Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

        GravitationalForceSunMercury(Sun,Mercury,Fx,Fy,Fz,epsilon);

        // Acceleration in each dimension for current planet
        acceleration[0] = Fx/Sun.mass;
        acceleration[1] = Fy/Sun.mass;
        acceleration[2] = Fz/Sun.mass;

        // Calculate new position for current planet
        for(int j=0; j<dimension; j++) {
            //cout << acceleration[j] << endl;
            Sun.position[j] += Sun.velocity[j]*timeStep + 0.5*timeStep*timeStep*acceleration[j];
        }

        // Update position so that the center of mass is always in the origin
  /*      centerOfMassVector = centerOfMass(dimension);
        Sun.position[0] -= centerOfMassVector[0];
        Sun.position[1] -= centerOfMassVector[1];
        Sun.position[2] -= centerOfMassVector[2];*/
        //cout << centerOfMassVector << endl;

//      double mercuryDistance = sqrt(Mercury.position[0]*Mercury.position[0] + Mercury.position[1]*Mercury.position[1] + Mercury.position[2]*Mercury.position[2]);

        double mercuryDistance = Mercury.distance(Sun);
        GravitationalForceSunMercury(Sun, Mercury, Fxnew, Fynew, Fznew, epsilon);

        // Acceleration each dimension exerted for current planet
        acceleration_new[0] = Fxnew/Sun.mass;
        acceleration_new[1] = Fynew/Sun.mass;
        acceleration_new[2] = Fznew/Sun.mass;

        // Calculate new velocity for current planet
        for(int j=0; j<dimension; j++){
            Sun.velocity[j] += 0.5*timeStep*(acceleration[j] + acceleration_new[j]);
        }

        // Write to file
        if(n%freq == 0 && doFileWriting == true){
            writeMercuryToFile(mercuryDistance, Sun, Mercury, time, integrationPoints);
        }
        n++;
        time += timeStep;
        cout << time << endl;
    }


    //Stop clock and print out time usage
    finish_VV = clock();
    std::cout << "Total time Verlet = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time

    std::cout << "One time step = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC)/integrationPoints << " seconds" << std::endl; // print elapsed time


    //Clear memory
    //delete_matrix(acceleration);
    //delete_matrix(acceleration_new);
}

vec solver::centerOfMass(int dimension){
    vec centerOfMass(3);
    centerOfMass[0] = centerOfMass[1] = centerOfMass[2] = 0;
    for(int i=0; i<total_planets; i++){
        planet &current = all_planets[i];
        for(int j=0; j<dimension; j++){
            centerOfMass[j] += (current.mass*current.position[j])/total_mass;
        }
    }
    return centerOfMass;
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

void solver::GravitationalForceSunMercury(planet &current,planet &other,double &Fx,double &Fy,double &Fz,double epsilon)
{   // Function that calculates the gravitational force between two objects, component by component.

    // Calculate relative distance between current planet and all other planets
    double relative_distance[3];

    for(int j = 0; j < 3; j++) relative_distance[j] = current.position[j]-other.position[j];
    double r = current.distance(other);
    double smoothing = epsilon*epsilon*epsilon;
    double l = current.AngularMomentum(other);
    //double c = 63239.7263; // AU/Year
    double c = 1.;
    // Calculate the forces in each direction
    Fx -= (this->G*current.mass*other.mass*relative_distance[0]/((r*r*r))*(1+ 3*l*l/(r*r*c*c)) + smoothing);
    Fy -= (this->G*current.mass*other.mass*relative_distance[1]/((r*r*r))*(1+ 3*l*l/(r*r*c*c)) + smoothing);
    Fz -= (this->G*current.mass*other.mass*relative_distance[2]/((r*r*r))*(1+ 3*l*l/(r*r*c*c)) + smoothing);
}

//void solver::GravitationalForce_RK(double x_rel,double y_rel,double z_rel,double &Fx,double &Fy,double &Fz,double mass1,double mass2)
//{   // Function that calculates the gravitational force between two objects, component by component.

//    // Calculate relative distance between current planet and all other planets
//    double r = sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel);

//    // Calculate the forces in each direction
//    Fx -= this->G*mass1*mass2*x_rel/(r*r*r);
//    Fy -= this->G*mass1*mass2*y_rel/(r*r*r);
//    Fz -= this->G*mass1*mass2*z_rel/(r*r*r);
//}

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
        for(int nr2=0;nr2<total_planets;nr2++){
            if (nr2!= nr1){
                planet &Other = all_planets[nr2];
                Current.potential += Current.PotentialEnergy(Other,G,epsilon);
                Other.potential += Other.PotentialEnergy(Current,G,epsilon);

        }

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
    string planetPositionPath= string("/uio/hume/student-u69/pederbh/FYS4150/Project3/PlanetPosition") + type + "Results.txt";
    ofstream myPlanetPositionFile;
    myPlanetPositionFile.open(planetPositionPath,std::ios::app);

    string planetEnergiesPath= string("/uio/hume/student-u69/pederbh/FYS4150/Project3/PlanetEnergies") + type + "Results.txt";
    ofstream myPlanetEnergiesFile(planetEnergiesPath);

    string mercuryPositionPath= string("/uio/hume/student-u69/pederbh/FYS4150/Project3/MercuryPositionResults.txt");
    ofstream myMercuryPositionFile;
    myMercuryPositionFile.open(mercuryPositionPath,std::ios::app);

    myPlanetPositionFile << integrationPoints << endl;
    myPlanetPositionFile << dim << endl;
    myPlanetPositionFile << total_planets << endl;

    myPlanetEnergiesFile << integrationPoints << endl;
    myPlanetEnergiesFile << dim << endl;
    myPlanetEnergiesFile << total_planets << endl;

    myMercuryPositionFile << integrationPoints << endl;
    myMercuryPositionFile << dim << endl;


    myPlanetPositionFile.close();
    myPlanetEnergiesFile.close();
    myMercuryPositionFile.close();
}

void solver::writeToFile(string type, planet current, double time, int integrationPoints, double kineticEnergy, double potentialEnergy, double angularMomentum) {
    string planetPositionPath= string("/uio/hume/student-u69/pederbh/FYS4150/Project3/PlanetPosition") + type + "Results.txt";
    ofstream myPlanetPositionFile;
    myPlanetPositionFile.open(planetPositionPath,std::ios::app);

    string planetEnergiesPath= string("/uio/hume/student-u69/pederbh/FYS4150/Project3/PlanetEnergies") + type + "Results.txt";
    ofstream myPlanetEnergiesFile;
    myPlanetEnergiesFile.open(planetEnergiesPath,std::ios::app);

    myPlanetPositionFile << current.position[0] << " " << current.position[1] << " " << current.position[2] << endl;
    myPlanetEnergiesFile << time << " " << kineticEnergy << " " << potentialEnergy << " " << angularMomentum << endl;

    myPlanetPositionFile.close();
    myPlanetEnergiesFile.close();
}

void solver::writeMercuryToFile(double mercuryDistance, planet &Sun, planet &Mercury, double time, int integrationPoints) {

    string mercuryPositionPath= string("/uio/hume/student-u69/pederbh/FYS4150/Project3/MercuryPositionResults.txt");

    if (! myMercuryPositionFile.is_open()) {
        myMercuryPositionFile.open(mercuryPositionPath,std::ios::out);
    }


    myMercuryPositionFile << setprecision(16) << time << " " << mercuryDistance << " " << Sun.position[0] << " " << Sun.position[1] << " " << Sun.position[2] << " " << Mercury.position[0] << " " << Mercury.position[1] << " " << Mercury.position[2] << endl;


    //myMercuryPositionFile.close();
}
