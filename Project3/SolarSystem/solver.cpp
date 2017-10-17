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

void::solver::EulerEarthSun(int integrationPoints, double finalTime, double initialEarthPosition_x, double intialEarthPosition_y, double initialVelocity_x, double initialVelocity_y)
{ /*Euler method to solve two coupled ODEs.*/

    // Define time step
    double timeStep = finalTime/((double) integrationPoints);
    double time = 0.0;
    double gravitationalConstant = 39.4784; // (AU)³/solarMass*yr²
    double massEarth = 3e-6; //solar mass
    double massSun = 1.0;
    double epsilon = 0.0;

    mat accelerationEarth(integrationPoints,2), velocityEarth(integrationPoints,2), positionEarth(integrationPoints,2), forceOnEarth(integrationPoints,2);

    //Initial conditions
    positionEarth(0,0) = initialEarthPosition_x;
    positionEarth(0,1) = intialEarthPosition_y;

    double initialDistance = sqrt(positionEarth(0,0)*positionEarth(0,0) + positionEarth(0,1)*positionEarth(0,1));
    //double initialVelocity = sqrt(gravitationalConstant*massSun/initialDistance);

    velocityEarth(0,0) = initialVelocity_x;
    velocityEarth(0,1) = initialVelocity_y;

    //Write kinetic and potential energy to file
    string path1="/uio/hume/student-u69/pederbh/FYS4150/Project3/TestingEnergiesEulerEarthSunResults.txt";
    ofstream myfile(path1);
    myfile << integrationPoints << endl;
    planet planet2(1.,0.,0.,0.,0.,0.,0.);           // Sun: (mass,x,y,z,vx,vy,vz)

    double startClock = clock();
    for (int i = 0; i < integrationPoints-1; i++){
        //Distance between the earth and the sun
        double distance = sqrt(positionEarth(i,0)*positionEarth(i,0) + positionEarth(i,1)*positionEarth(i,1));

        //Gravitational force in x and y direction
        double F_ = -(gravitationalConstant*massSun*massEarth)/(distance*distance*distance); // forceOnEarth % distance
        forceOnEarth(i,0) = positionEarth(i,0)*F_;
        forceOnEarth(i,1) = positionEarth(i,1)*F_;

        //Motion of the earth
        accelerationEarth(i,0) = forceOnEarth(i,0)/massEarth;
        accelerationEarth(i,1) = forceOnEarth(i,1)/massEarth;
        velocityEarth(i+1,0) = velocityEarth(i,0) + accelerationEarth(i,0)*timeStep;
        velocityEarth(i+1,1) = velocityEarth(i,1) + accelerationEarth(i,1)*timeStep;
        positionEarth(i+1,0) = positionEarth(i,0) + velocityEarth(i+1,0)*timeStep;
        positionEarth(i+1,1) = positionEarth(i,1) + velocityEarth(i+1,1)*timeStep;

        planet planet1(massEarth,positionEarth(i+1,0),positionEarth(i+1,1),0.0,velocityEarth(i+1,0),velocityEarth(i+1,1),0.); //Earth, testing to find circular orbit
        double kineticEnergyEarth = planet1.KineticEnergy();
        double potentialEnergyEarth = planet1.PotentialEnergy(planet2, gravitationalConstant, epsilon);
        double angularMomentumEarth = planet1.AngularMomentum(planet2);

        time += timeStep;

        myfile << time << " " << kineticEnergyEarth << " " << potentialEnergyEarth << " " << angularMomentumEarth << endl;
    }
    double elapsedTime = clock() - startClock;
    std::cout << "Total time Euler = " << "\t" << ((float)(elapsedTime)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time


    myfile.close();
    //Print to output file
    string path2="/uio/hume/student-u69/pederbh/FYS4150/Project3/TestingEulerEarthSunResults.txt";
    ofstream myfile2(path2);
    myfile2 << integrationPoints << endl;
    for(int i=0; i<integrationPoints; i++){
        myfile2 << positionEarth(i,0) << " " << positionEarth(i,1) << endl;
    }
    myfile2.close();
}

void solver::VelocityVerlet(int dimension, int integration_points, double final_time, int print_number, double epsilon)
{   /*  Velocity-Verlet solver for two coupeled ODEs in a given number of dimensions.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/


    // Define time step
    double time_step = final_time/((double) integration_points);
    double time = 0.0;
    double loss = 0.; // Possible energy loss
    int lostPlanets[integration_points];

//    // Create files for data storage
//    char *filename = new char[1000];
//    char *filenameE = new char[1000];
//    char *filenameB = new char[1000];
//    char *filenameLost = new char[1000];
//        sprintf(filename, "PlanetsVV_%d_%.3f.txt",total_planets,time_step);
//        sprintf(filenameE, "PlanetsVV_energy_%d_%.3f.txt",total_planets,time_step);
//        sprintf(filenameB,"Planetsbound_%d_%.3f.txt",total_planets,time_step);
//        sprintf(filenameLost,"Planetslost_%d_%.3f.txt",total_planets,time_step);
//    std::ofstream output_file(filename);
//    std::ofstream output_energy(filenameE);
//    std::ofstream output_bound(filenameB);
//    std::ofstream output_lost(filenameLost);

    // Set up arrays
    double **acceleration = setup_matrix(total_planets,3);
    double **acceleration_new = setup_matrix(total_planets,3);

    // Initialize forces
    double Fx,Fy,Fz,Fxnew,Fynew,Fznew; // Forces in each dimension

    // Write initial values to file
//    print_position(output_file,dimension,time,print_number);
//    print_energy(output_energy,time,epsilon);

    int n = 0;
    lostPlanets[n] = 0;
//    output_lost << time << "\t" << lostPlanets[n] << std::endl;
//    n+=1;

    // Set up clock to measure the time usage
    clock_t planet_VV,finish_VV;
    planet_VV = clock();

    // PLANET CALCULATIONS
    // Loop over time
    time += time_step;

//    vector<ofstream> a;
//    a.push_back();

    string earthPath=string("/uio/hume/student-u69/pederbh/FYS4150/Project3/EarthPositionVerletResults.txt");
    ofstream myEarthFile(earthPath);
    myEarthFile << integration_points << endl;

    while(time < final_time){
        lostPlanets[n] = 0;
        // Loop over all planets
        for(int nr1=0; nr1<total_planets; nr1++){
            planet &current = all_planets[nr1]; // Current planet we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Calculate forces in each dimension
                for(int nr2=nr1+1; nr2<total_planets; nr2++){
                    planet &other = all_planets[nr2];
                    GravitationalForce(current,other,Fx,Fy,Fz,epsilon);
                }

            // Acceleration in each dimension for current planet
            acceleration[nr1][0] = Fx/current.mass;
            acceleration[nr1][1] = Fy/current.mass;
            acceleration[nr1][2] = Fz/current.mass;

            // Calculate new position for current planet
            for(int j=0; j<dimension; j++) {
                current.position[j] += current.velocity[j]*time_step + 0.5*time_step*time_step*acceleration[nr1][j];
            }

            // Loop over all other planets
                for(int nr2=nr1+1; nr2<total_planets; nr2++){
                    planet &other = all_planets[nr2];
                    GravitationalForce(current,other,Fxnew,Fynew,Fznew,epsilon);
                }

            // Acceleration each dimension exerted for current planet
            acceleration_new[nr1][0] = Fxnew/current.mass;
            acceleration_new[nr1][1] = Fynew/current.mass;
            acceleration_new[nr1][2] = Fznew/current.mass;

            // Calculate new velocity for current planet
            for(int j=0; j<dimension; j++) current.velocity[j] += 0.5*time_step*(acceleration[nr1][j] + acceleration_new[nr1][j]);

            //Print to output file
            if(nr1==0){
                myEarthFile << current.position[0] << " " << current.position[1] << endl;
            }
        }



        // Energy conservation

        // Write current values to file and increase time
//        print_position(output_file,dimension,time,print_number);
//        print_energy(output_energy,time,epsilon);

        loss += EnergyLoss();

        for(int nr=0;nr<total_planets;nr++){
            planet &Current = all_planets[nr];
            if(!(this->Bound(Current))){
                lostPlanets[n] += 1;
            }
        }
//        output_lost << time << "\t" << lostPlanets[n] << std::endl;
        n += 1;
        time += time_step;
    }
    // Stop clock and print out time usage
    finish_VV = clock();
    std::cout << "Total time Verlet = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time
    myEarthFile.close();

//    std::cout << "One time step = " << "\t" << ((float)(finish_VV - planet_VV)/CLOCKS_PER_SEC)/integration_points << " seconds" << std::endl; // print elapsed time

//    //loss = EnergyLoss();
//    std::cout << "Total energyloss due to unbound planets: " << loss << std::endl;

//    double boundPlanets = 0;
//    for(int nr=0;nr<total_planets;nr++){
//        planet &Current = all_planets[nr];
//        if(this->Bound(Current)){
//            output_bound << nr << std::endl;
//            boundPlanets += 1;
//        }
//    }
//    std::cout << "There are " << boundPlanets << " bound planets at the end of the run" << std::endl;

//    // Close files
//    output_file.close();
//    output_energy.close();
//    output_bound.close();
//    output_lost.close();

    // Clear memory
    delete_matrix(acceleration);
    delete_matrix(acceleration_new);
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
