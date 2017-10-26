#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>
#include "armadillo"
using std::vector;
using namespace arma;

class solver
{
public:
    friend class planet;

    // properties
    double radius,total_mass,G;
    int total_planets;
    vector<planet> all_planets;
    double totalKinetic;
    double totalPotential;
    bool doFileWriting = true;
    int freq           = 10;

    std::ofstream myMercuryPositionFile;

    void setFileWriting(bool w, int f) {
        doFileWriting = w;
        freq = f;
    }
    // constants

    // initializers
    solver();
    solver(double radi);

    // functions
    void add(planet newplanet);
    void addM(planet newplanet);
    void GravitationalConstant();
    void print_position(std::ofstream &output, int dimension, double time, int number);
    void print_energy(std::ofstream &output, double time, double epsilon);
    void VelocityVerlet(int dimension, int integration_points, double final_time, double epsilon);//(int dimension, int integrationPoints, int finalTime, double epsilon);//
    void VelocityVerletMercury(int dimension, int integration_points, double final_time, double epsilon, planet &Sun, planet &Mercury);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz);
    void GravitationalForceRelativistic(planet &current, planet &other, double &Fx, double &Fy, double &Fz);
    void GravitationalForce_RK(double x_rel, double y_rel, double z_rel, double &Fx, double &Fy, double &Fz, double mass1, double mass2);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    void Euler(int dimension, int integrationPoints, int finalTime, double epsilon);
    void writeInformationToFile(std::string type, int integrationPoints, int dim);
    void writeMercuryToFile(double mercuryDistance, planet &Sun, planet &Mercury, double time, int integrationPoints);
    void writeToFile(std::string type, planet& current, double time, int integrationPoints, double kineticEnergy, double potentialEnergy, double angularMomentum);
    vec centerOfMass(int dimension);
};

#endif // SOLVER_H
