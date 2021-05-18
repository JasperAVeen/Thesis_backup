/*    Copyright (c) 2010-2019, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <string>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "G:/Thesis/applicationOutput.h"

#include "Classes/mgatrajectory.h"
#include "Problems/applicationOutput.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Classes/tt_functions.h"
#include "Classes/maketraj.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Classes/lagrangepoints.h"
#include "Classes/analysisb_functions.h"



#include <Tudat/InputOutput/basicInputOutput.h>
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"



//Toggles

//! Execute  main
int main( )
{

    std::vector<std::string> alternativeKernel {"thesis_Lagrange_kernels.bsp"};
    spice_interface::loadStandardSpiceKernels( alternativeKernel );
    std::string outputPath = tudat_applications::getOutputPath( ) + "Tester/Rendezvous/";

    std::string system = "Em";
    std::string LPO_type = "halo";
    std::string LP = "L1";
    std::pair<int, int> orbitID_bounds = std::make_pair(0, 2500);

    std::unordered_map<int, std::map<double, Eigen::Vector6d>> allLPOs = importLPO_cr3b(system, LPO_type, LP, orbitID_bounds);

    //Example variables
    int orbitID = 1500;
    double orbitFraction = 0.5;
    double arrivalTime = 1200 * physical_constants::JULIAN_DAY;
    double stayTime = 20;
    double departureTime = arrivalTime + (stayTime * physical_constants::JULIAN_DAY);

    Eigen::Vector6d LPOcart_Arrstate = getLPO_ArrState(system, allLPOs, orbitID, orbitFraction, arrivalTime);
    Eigen::Vector6d LPOcart_Depstate = getLPO_DepState(system, allLPOs, orbitID, orbitFraction, departureTime, stayTime);

    std::map<double, Eigen::Vector6d> LPO_orbit = getLPO_cart(system, allLPOs, orbitID, arrivalTime);


    std::map<int, Eigen::Vector6d> boundaryStates {};
    boundaryStates[1] = LPOcart_Arrstate;
    boundaryStates[2] = LPOcart_Depstate;

    input_output::writeDataMapToTextFile(boundaryStates, "gateway_rendezvous.dat", outputPath );
    input_output::writeDataMapToTextFile(LPO_orbit, "LPO_map.dat", outputPath );



    /*
    for (int orbitID {0}; orbitID < 1177; orbitID++) {
    std::map<double, Eigen::Vector6d> Map2 = getLPO_cart ("Em", "horizontal", "L1",  orbitID, 1000, true) ;
    }
    */





    /*
    std::vector<std::string> alternativeKernel {"thesis_Lagrange_kernels.bsp"};
    spice_interface::loadStandardSpiceKernels( alternativeKernel );
    std::cout << "\nSuccesfully loaded kernels!\n";
    std::string outputPath = tudat_applications::getOutputPath( ) + "Tester/LibPoints/Final/";

    //initializing storage
    std::map<double, Eigen::Vector6d> SE_L1;
    std::map<double, Eigen::Vector6d> SE_L2;
    std::map<double, Eigen::Vector6d> SE_L3;
    std::map<double, Eigen::Vector6d> SE_L4;
    std::map<double, Eigen::Vector6d> SE_L5;

    std::map<double, Eigen::Vector6d> Em_L1;
    std::map<double, Eigen::Vector6d> Em_L2;
    std::map<double, Eigen::Vector6d> Em_L3;
    std::map<double, Eigen::Vector6d> Em_L4;
    std::map<double, Eigen::Vector6d> Em_L5;

    std::map<double, Eigen::Vector6d> Em_L1_E;
    std::map<double, Eigen::Vector6d> Em_L2_E;
    std::map<double, Eigen::Vector6d> Em_L3_E;
    std::map<double, Eigen::Vector6d> Em_L4_E;
    std::map<double, Eigen::Vector6d> Em_L5_E;

    std::map<double, Eigen::Vector6d> Sun_EB;
    std::map<double, Eigen::Vector6d> Sun_Earth;
    std::map<double, Eigen::Vector6d> Earth_Moon;



    //creating Ephemeris vecs
    std::vector<int> SunEB = {3};
    std::vector<int> SunMB = {4};
    std::vector<int> SunE = {3};
    std::vector<int> vecEm = {0};




    double initial_time{1000 * physical_constants::JULIAN_DAY};
    double end_time{initial_time + 5 * physical_constants::JULIAN_YEAR};
    double time_step{physical_constants::JULIAN_DAY};

    for(double time {initial_time}; time < end_time; time += time_step) {
        SE_L1[time] << getLibrationState("SE", 1, time, false);
        SE_L2[time] << getLibrationState("SE", 2, time, false);
        SE_L3[time] << getLibrationState("SE", 3, time, false);
        SE_L4[time] << getLibrationState("SE", 4, time, false);
        SE_L5[time] << getLibrationState("SE", 5, time, false);


        Em_L1[time] << getLibrationState("Em", 1, time, true);
        Em_L2[time] << getLibrationState("Em", 2, time, true);
        Em_L3[time] << getLibrationState("Em", 3, time, true);
        Em_L4[time] << getLibrationState("Em", 4, time, true);
        Em_L5[time] << getLibrationState("Em", 5, time, true);


        Em_L1_E[time] << getLibrationState("Em", 1, time, false);
        Em_L2_E[time] << getLibrationState("Em", 2, time, false);
        Em_L3_E[time] << getLibrationState("Em", 3, time, false);
        Em_L4_E[time] << getLibrationState("Em", 4, time, false);
        Em_L5_E[time] << getLibrationState("Em", 5, time, false);

        Sun_EB[time] << createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true)[0]->getCartesianState( time );
        Sun_Earth[time] << createEphemerisVec(SunE, "SUN", "ECLIPJ2000", false)[0]->getCartesianState( time );
        Earth_Moon[time] << createEphemerisVec(vecEm, "EARTH", "ECLIPJ2000", true)[0]->getCartesianState( time );
    }


    input_output::writeDataMapToTextFile(SE_L1, "SE_L1.dat", outputPath );
    input_output::writeDataMapToTextFile(SE_L2, "SE_L2.dat", outputPath );
    input_output::writeDataMapToTextFile(SE_L3, "SE_L3.dat", outputPath );
    input_output::writeDataMapToTextFile(SE_L4, "SE_L4.dat", outputPath );
    input_output::writeDataMapToTextFile(SE_L5, "SE_L5.dat", outputPath );

    input_output::writeDataMapToTextFile(Em_L1, "Em_L1.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L2, "Em_L2.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L3, "Em_L3.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L4, "Em_L4.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L5, "Em_L5.dat", outputPath );

    input_output::writeDataMapToTextFile(Em_L1_E, "Em_L1_E.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L2_E, "Em_L2_E.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L3_E, "Em_L3_E.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L4_E, "Em_L4_E.dat", outputPath );
    input_output::writeDataMapToTextFile(Em_L5_E, "Em_L5_E.dat", outputPath );

    input_output::writeDataMapToTextFile(Sun_EB, "Sun_EB.dat", outputPath );
    input_output::writeDataMapToTextFile(Sun_Earth, "Sun_Earth.dat", outputPath );
    input_output::writeDataMapToTextFile(Earth_Moon, "Earth_Moon.dat", outputPath );















    double sunGravitationalParameter = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
    double earthGravitationalParameter = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
    double EBGravitationalParameter = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;


    std::vector<int> SunEB = {3}; std::vector<std::shared_ptr< ephemerides::Ephemeris >> Sun_EB_Ephemeris = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true);


    circular_restricted_three_body_problem::LagrangePoint SE ( sunGravitationalParameter, EBGravitationalParameter,
                                                                  std::make_shared< root_finders::NewtonRaphson >( 1.0e-14, 1000 ) );



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Creating positions of Lagrange points in Corotating frame
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map<int, Eigen::Vector6d> relposMap;

    for (int i {1}; i < 6; i++) {
    SE.computeLocationOfLibrationPoint( i );
    Eigen::Vector3d relPosVec = SE.getLocationOfLagrangeLibrationPoint( );
    relposMap[i] << relPosVec, 0, 0, 0;
}
    input_output::writeDataMapToTextFile(relposMap, "relposMap.dat", outputPath );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Creating positions of Lagrange points in Cartesian
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map<double, Eigen::Vector6d> cartposMap_L1;
    std::map<double, Eigen::Vector6d> cartposMap_L2;
    std::map<double, Eigen::Vector6d> cartposMap_L3;
    std::map<double, Eigen::Vector6d> cartposMap_L4;
    std::map<double, Eigen::Vector6d> cartposMap_L5;
    std::map<double, Eigen::Vector6d> cartposMap_EB;



    double initial_time{1000 * physical_constants::JULIAN_DAY};
    double end_time{initial_time + 5 * physical_constants::JULIAN_YEAR};
    double time_step{physical_constants::JULIAN_DAY};

    for(double time {initial_time}; time < end_time; time += time_step) {
        Eigen::Vector6d stateEB_Sun = Sun_EB_Ephemeris[0]->getCartesianState( time );

        cartposMap_EB [time] << stateEB_Sun;

        //angles
        double theta = atan2(stateEB_Sun(1), stateEB_Sun(0));
        double alfa{};
        if(stateEB_Sun(2)<0){
            alfa = std::acos(sqrt(pow(stateEB_Sun(0),2) + pow(stateEB_Sun(1),2))/stateEB_Sun.segment(0,3).norm());
        }
        else{
            alfa = - std::acos(sqrt(pow(stateEB_Sun(0),2) + pow(stateEB_Sun(1),2))/stateEB_Sun.segment(0,3).norm());
        }

        double beta{};
        if(stateEB_Sun(5)>0){
            beta = std::acos(sqrt(pow(stateEB_Sun(3),2) + pow(stateEB_Sun(4),2))/stateEB_Sun.segment(3,3).norm());
        }
        else{
            beta = - std::acos(sqrt(pow(stateEB_Sun(3),2) + pow(stateEB_Sun(4),2))/stateEB_Sun.segment(3,3).norm());
        }

        Eigen::Vector6d Cartesian_L1 = circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis( sunGravitationalParameter, earthGravitationalParameter,
                 stateEB_Sun.norm(), relposMap[1], theta, alfa, beta);
        cartposMap_L1[time] << Cartesian_L1;

        Eigen::Vector6d Cartesian_L2 = circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis( sunGravitationalParameter, earthGravitationalParameter,
                 stateEB_Sun.norm(), relposMap[2], theta, alfa, beta);
        cartposMap_L2[time] << Cartesian_L2;

        Eigen::Vector6d Cartesian_L3 = circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis( sunGravitationalParameter, earthGravitationalParameter,
                 stateEB_Sun.norm(), relposMap[3], theta, alfa, beta);
        cartposMap_L3[time] << Cartesian_L3;

        Eigen::Vector6d Cartesian_L4 = circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis( sunGravitationalParameter, earthGravitationalParameter,
                 stateEB_Sun.norm(), relposMap[4], theta, alfa, beta);
        cartposMap_L4[time] << Cartesian_L4;

        Eigen::Vector6d Cartesian_L5 = circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis( sunGravitationalParameter, earthGravitationalParameter,
                 stateEB_Sun.norm(), relposMap[5], theta, alfa, beta);
        cartposMap_L5[time] << Cartesian_L5;
        }

    input_output::writeDataMapToTextFile(cartposMap_L1, "cartposMap_L1.dat", outputPath );
    input_output::writeDataMapToTextFile(cartposMap_L2, "cartposMap_L2.dat", outputPath );
    input_output::writeDataMapToTextFile(cartposMap_L3, "cartposMap_L3.dat", outputPath );
    input_output::writeDataMapToTextFile(cartposMap_L4, "cartposMap_L4.dat", outputPath );
    input_output::writeDataMapToTextFile(cartposMap_L5, "cartposMap_L5.dat", outputPath );
    input_output::writeDataMapToTextFile(cartposMap_EB, "cartposMap_EB.dat", outputPath );


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Creating positions of Lagrange points from Spice
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map<double, Eigen::Vector6d> SELib1_spice;
    std::map<double, Eigen::Vector6d> SELib2_spice;
    std::map<double, Eigen::Vector6d> SELib4_spice;
    std::map<double, Eigen::Vector6d> SELib5_spice;

    int timeincr2{0};

    for(double time{initial_time}; time < end_time; time +=time_step){
    SELib1_spice[time] <<  spice_interface::getBodyCartesianStateAtEpoch("391", "EARTH_BARYCENTER", "ECLIPJ2000", "NONE", time);
    SELib2_spice[time] << spice_interface::getBodyCartesianStateAtEpoch("392", "EARTH_BARYCENTER", "ECLIPJ2000", "NONE", time);
    SELib4_spice[time] << spice_interface::getBodyCartesianStateAtEpoch("394", "10", "ECLIPJ2000", "NONE", time);
    SELib5_spice[time] << spice_interface::getBodyCartesianStateAtEpoch("395", "10", "ECLIPJ2000", "NONE", time);
    std::cout << timeincr2 <<std::endl;
    timeincr2++;
    }

    input_output::writeDataMapToTextFile(SELib1_spice, "SELib1_spice.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib2_spice, "SELib2_spice.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib4_spice, "SELib4_spice.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib5_spice, "SELib5_spice.dat", outputPath );


    //Compute position of Lagrange points through Tudat
    double initial_time{1000 * physical_constants::JULIAN_DAY};
    double end_time{initial_time + 1 * physical_constants::JULIAN_YEAR};
    double time_step{physical_constants::JULIAN_DAY};

    std::map<double, Eigen::Vector6d> SELib1;
    std::map<double, Eigen::Vector6d> SELib2;
    std::map<double, Eigen::Vector6d> SELib4;
    std::map<double, Eigen::Vector6d> SELib5;

    int timeincr{1};

    for(double time{initial_time}; time < end_time; time +=time_step){
    Eigen::Vector6d SEL1 = getLibrationState("SE", 1, time, true);
    Eigen::Vector6d SEL2 = getLibrationState("SE", 2, time, true);
    Eigen::Vector6d SEL4 = getLibrationState("SE", 4, time, true);
    Eigen::Vector6d SEL5 = getLibrationState("SE", 5, time, true);
    SELib1[time] << SEL1;
    SELib2[time] << SEL2;
    SELib4[time] << SEL4;
    SELib5[time] << SEL5;
    std::cout << timeincr <<std::endl;
    timeincr++;
   }

    input_output::writeDataMapToTextFile(SELib1, "SELib1.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib2, "SELib2.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib4, "SELib4.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib5, "SELib5.dat", outputPath );





    std::map<double, Eigen::Vector6d> EmLib1;
    std::map<double, Eigen::Vector6d> EmLib2;
    std::map<double, Eigen::Vector6d> EmLib4;
    std::map<double, Eigen::Vector6d> EmLib5;

    int timeincr{1};

    for(double time{initial_time}; time < end_time; time +=time_step){
    Eigen::Vector6d EmL1 = getLibrationState("Em", 1, time, true);
    Eigen::Vector6d EmL2 = getLibrationState("Em", 2, time, true);
    Eigen::Vector6d EmL4 = getLibrationState("Em", 4, time, true);
    Eigen::Vector6d EmL5 = getLibrationState("Em", 5, time, true);
    EmLib1[time] << EmL1;
    EmLib2[time] << EmL2;
    EmLib4[time] << EmL4;
    EmLib5[time] << EmL5;
    std::cout << timeincr <<std::endl;
    timeincr++;
   }

    input_output::writeDataMapToTextFile(EmLib1, "EmLib1.dat", outputPath );
    input_output::writeDataMapToTextFile(EmLib2, "EmLib2.dat", outputPath );
    input_output::writeDataMapToTextFile(EmLib4, "EmLib4.dat", outputPath );
    input_output::writeDataMapToTextFile(EmLib5, "EmLib5.dat", outputPath );




   //Compute positions of Lagrange points through Spice
    std::map<double, Eigen::Vector6d> SELib1_spice;
    std::map<double, Eigen::Vector6d> SELib2_spice;
    std::map<double, Eigen::Vector6d> SELib4_spice;
    std::map<double, Eigen::Vector6d> SELib5_spice;

    int timeincr2{0};

    for(double time{initial_time}; time < end_time; time +=time_step){
    SELib1_spice[time] <<  spice_interface::getBodyCartesianStateAtEpoch("391", "EARTH_BARYCENTER", "ECLIPJ2000", "NONE", time);
    SELib2_spice[time] << spice_interface::getBodyCartesianStateAtEpoch("392", "EARTH_BARYCENTER", "ECLIPJ2000", "NONE", time);
    SELib4_spice[time] << spice_interface::getBodyCartesianStateAtEpoch("394", "10", "ECLIPJ2000", "NONE", time);
    SELib5_spice[time] << spice_interface::getBodyCartesianStateAtEpoch("395", "10", "ECLIPJ2000", "NONE", time);
    std::cout << timeincr2 <<std::endl;
    timeincr2++;
    }

    input_output::writeDataMapToTextFile(SELib1_spice, "SELib1_spice.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib2_spice, "SELib2_spice.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib4_spice, "SELib4_spice.dat", outputPath );
    input_output::writeDataMapToTextFile(SELib5_spice, "SELib5_spice.dat", outputPath );



    //Planets
    std::vector<int> SunEB = {3};
    std::map<double, Eigen::Vector6d> EB_SSB;
    std::map<double, Eigen::Vector6d> EB_SUN;

    std::map<double, Eigen::Vector6d> SUN_SSB;


    for(double time{initial_time}; time < end_time; time +=time_step){
    EB_SSB[time] <<  spice_interface::getBodyCartesianStateAtEpoch("EARTH_BARYCENTER", "SSB", "ECLIPJ2000", "NONE", time);
    EB_SUN[time] <<  createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true)[0]->getCartesianState( time );

    SUN_SSB[time] <<  spice_interface::getBodyCartesianStateAtEpoch("SUN", "SSB", "ECLIPJ2000", "NONE", time);
    }
    input_output::writeDataMapToTextFile(EB_SUN, "EB_SUN.dat", outputPath );
    input_output::writeDataMapToTextFile(EB_SSB, "EB_SSB.dat", outputPath );

    input_output::writeDataMapToTextFile(SUN_SSB, "SUN_SSB.dat", outputPath );

*/
    std::cout << "\n\n------------------------\nProgram ended succesfully!\n--------------------------\n\n";

    return 0;






}
