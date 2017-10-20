
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
        lostPlanets[n] = 0;
        // Loop over all planets
        for(int nr1=0; nr1<total_planets; nr1++){
            planet &current = all_planets[nr1]; // Current planet we are looking at

            Fx = Fy = Fz = Fxnew = Fynew = Fznew = 0.0; // Reset forces before each run

            // Calculate forces in each dimension
                for(int nr2=nr1+1; nr2<total_planets; nr2++){
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
                current.position[j] += current.velocity[j]*timeStep ;//+ 0.5*timeStep*timeStep*acceleration[nr1][j];
            }

//            // Loop over all other planets
//                for(int nr2=nr1+1; nr2<total_planets; nr2++){
//                    if(nr1!=nr2) {
//                        planet &other = all_planets[nr2];
//                        GravitationalForce(current,other,Fxnew,Fynew,Fznew,epsilon);
//                    }
//                }

//            // Acceleration each dimension exerted for current planet
//            acceleration_new[nr1][0] = Fxnew/current.mass;
//            acceleration_new[nr1][1] = Fynew/current.mass;
//            acceleration_new[nr1][2] = Fznew/current.mass;

            // Calculate new velocity for current planet
            for(int j=0; j<dimension; j++){
                current.velocity[j] += timeStep*acceleration[nr1][j];//0.5*timeStep*(acceleration[nr1][j] + acceleration_new[nr1][j]);
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
                writeToFile("Verlet", current, time+timeStep, kineticEnergy, potentialEnergy, angularMomentum);
            }
        }

        loss += EnergyLoss();

        for(int nr=0;nr<total_planets;nr++){
            planet &Current = all_planets[nr];
            if(!(this->Bound(Current))){
                lostPlanets[n] += 1;
            }
        }

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
