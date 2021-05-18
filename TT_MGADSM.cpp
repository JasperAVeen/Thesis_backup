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

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "applicationOutput.h"

#include "Classes/mgatrajectory.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Classes/tt_functions.h"
#include "Classes/maketraj.h"


#include <Tudat/InputOutput/basicInputOutput.h>
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"


//Toggles
//! Execute  main
int main( )
{   //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////Toggles and buttons//////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    //What trajectory should be optimized?
    std::string flyby_seq = "EdVdVdEdJdS";


    //What are the bounds wrt departure time and TOF?
    double start_date {-1000}; //MJD2000
    double launch_window {1000};   //days
    double end_date {1000000}; //MJD2000
    const double TOF_constraint = 30*365.25; //days

    //What are departure and capture conditions?
    double R_e {6371e3}; double R_GEO {42164e3}; double h_parking {185e3};
    double R_n{24766e3}; double h_na {488000e3}; double R_np {4000e3};

    double DepEcc{0};//{(R_GEO-(h_parking+R_e))/(R_GEO+(h_parking+R_e))};
    double DepSMA{std::numeric_limits< double >::infinity( )};//{((R_GEO)+(h_parking+R_e))/2};
    double CapEcc{0};//{((h_na+R_n)-(R_np+R_n))/((h_na+R_n)+(R_np+R_n))};
    double CapSMA{std::numeric_limits< double >::infinity( )};//{((h_na+R_n)+(R_np+R_n))/2};

    //Set seed:
    std::vector<unsigned int> seedvec = {123, 456, 789, 101112, 131415};
    std::map<unsigned int, double> DVperSeed;

    for (int j{0}; j < seedvec.size(); j++){
    unsigned int seed {seedvec.at(j)};
    std::string seedstring {std::to_string(seed)};

    std::cout << "==========\nSeed: " << seed <<"\n==========\n";


    //Toggles:
    const bool printBest = false;
    const bool useTripTime = false;

    const double sunGravitationalParameter = 1.32712428e20;


    //OutputPath
    std::string pathSuffix = "MGADSM/IndividualTrajectories/" + flyby_seq + "/" + seedstring;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    spice_interface::loadStandardSpiceKernels();


    // Create flyby sequence
    std::vector< int > flybySequence = createFlyBySeq(flyby_seq);

    //Create LegTypeVec
    std::vector< LegTypeTransfer > LegTypeVec =  createLegType (flyby_seq);

    // Create Bounds
    int num_DSMs = std::count(flyby_seq.begin(), flyby_seq.end(), 'd');
    size_t numberOfParameters = flybySequence.size()+(num_DSMs*4);
    std::vector< std::vector< double > > Bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );
    Bounds =  createCassiniBounds (flybySequence , start_date, launch_window, num_DSMs, Bounds, LegTypeVec);

    //create Departure and Capture conditions
    Eigen::VectorXd semiMajorAxes(2);
    semiMajorAxes << DepSMA, CapSMA;
    Eigen::VectorXd eccentricities(2);
    eccentricities << DepEcc, CapEcc;

    //Set seed for reproducible results
    pagmo::random_device::set_seed( seed );

    // Create object to compute the problem fitness
    problem prob{ MGAtrajectory( Bounds, flybySequence, LegTypeVec, flyby_seq, TOF_constraint, semiMajorAxes, eccentricities, end_date ) };


    // Select NSGA2 algorithm for problem
    algorithm algo{de(1u, 0.7, 0.97, 2u, 1e-10, 1e-10 )};

    // Create an island with individuals
    size_t numInd{60};
    island isl{algo, prob, numInd };

    // Evolve for numGen generations
    int numGen {50000};
    for( int i = 0 ; i < numGen+1; i++ )
    {

        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }

        isl.wait_check( ); // Raises errors

        // Write current iteration results to file
        printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga" + std::to_string( i ), false, pathSuffix );
        printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga" + std::to_string( i ), true, pathSuffix );
        if(i%10000==0){
        std::cout<<i<<"  -   ";
        std::cout<<isl.get_population().champion_f().at(0)<<std::endl;
}
    }


//Creating results for best individual

    int numberOfLegs = flybySequence.size( );


    auto TopIndividual = Eigen::VectorXd(numberOfLegs+1+(4*num_DSMs));
    for (size_t i = 0 ; i < isl.get_population().champion_x().size() ; i++ ){
        TopIndividual[i] = isl.get_population().champion_x()[i];
    }
    TopIndividual *= physical_constants::JULIAN_DAY;
    TopIndividual[ numberOfLegs ] = 1;//dummy

    //Write the variable vector for DSM variables below

    if(num_DSMs>0){
        int count {numberOfLegs+1};
        for(size_t i{0}; i < num_DSMs; i++){
            TopIndividual[ count ] = isl.get_population().champion_x()[ count-1 ];
            TopIndividual[ count + 1 ] = isl.get_population().champion_x()[ count ];
            TopIndividual[ count + 2 ] = isl.get_population().champion_x()[ count + 1 ];
            TopIndividual[ count + 3 ] = isl.get_population().champion_x()[ count + 2 ];
            count+=4;
        }
    }

    std::cout << "\n------------------------\nVariable Vector\n------------------------" << TopIndividual << "\n------------------------\n";
    Eigen::VectorXd minimumPericenterRadii{};
    minimumPericenterRadii.resize( flybySequence.size() );

    if(num_DSMs==0) {minimumPericenterRadii=createMinPerRad(flybySequence);}
    else {
        for (int i{0}; i < flybySequence.size(); i++){
            minimumPericenterRadii [i] = TUDAT_NAN;
        }

    }


    makeTraj BestmgaTraj = makeTraj( flybySequence.size(), LegTypeVec, flyby_seq, createEphemerisVec (flybySequence),
                                 createGravParamVec(flybySequence),
                                 TopIndividual,
                                 sunGravitationalParameter,
                                 minimumPericenterRadii,
                                 semiMajorAxes,
                                 eccentricities
                                 );



    // Vectors for the specific maneuvers and the total delta v
    std::vector< Eigen::Vector3d > positionVectorTop;
    std::vector< double > timeVectorTop;
    std::vector< double > deltaVVectorTop;
    double resultingDeltaVTop;

    // Calculate the orbits
    BestmgaTraj.calculateTrajectory( resultingDeltaVTop );

    BestmgaTraj.maneuvers( positionVectorTop, timeVectorTop, deltaVVectorTop );


    // Define vectors to calculate intermediate points
    std::vector< Eigen::Vector3d > interPositionVectorTop;
    std::vector< double > interTimeVectorTop;

    // Calculate intermediate points
    BestmgaTraj.intermediatePoints( 1000.0 , interPositionVectorTop, interTimeVectorTop );

    // Define vectors to calculate manoeuvres
    std::vector< Eigen::Vector3d > manPositionVectorTop;
    std::vector< double > manTimeVectorTop;
    std::vector< double > manDeltaVVectorTop;

    // Calculate maneuvers
    BestmgaTraj.maneuvers( manPositionVectorTop, manTimeVectorTop, manDeltaVVectorTop );


//Saving results for best individual

    // Writing intermediate points to file
    std::string outputFileTraj = tudat_applications::getOutputPath( ) + pathSuffix + "/TopTrajectory.dat";
    writeTrajectoryToFile( interPositionVectorTop, interTimeVectorTop, outputFileTraj );

    // Writing intermediate points to file
    std::string outputFileMan = tudat_applications::getOutputPath( ) + pathSuffix + "/TopManeuvers.dat";
    writeTrajectoryToFile( manPositionVectorTop, manTimeVectorTop, outputFileMan );



    std::ofstream file(tudat_applications::getOutputPath( ) + pathSuffix +  "/TopDeltaV.dat");
    for (size_t i {0}; i < deltaVVectorTop.size(); i++){
        file << i << " " <<deltaVVectorTop.at(i) << std::endl;
    }
    file.close();

    std::cout << "=================================================\n" << "Best individual: \n" << "=================================================\n";
    std::cout << " Total Delta V: " << resultingDeltaVTop <<std::endl << " Total TOF: " << (timeVectorTop[ timeVectorTop.size() -1 ] - timeVectorTop[ 0 ])/physical_constants::JULIAN_DAY << " days"
              << std::endl << " TOF in years: " << (timeVectorTop[ timeVectorTop.size() -1 ] - timeVectorTop[ 0 ])/physical_constants::JULIAN_YEAR << std::endl <<  "=================================================\n";


    for(int i{0}; i<timeVectorTop.size(); i++){
        if(i>0) std::cout << "Planet: " << i << " --- " << "Time: " << timeVectorTop[i]/physical_constants::JULIAN_DAY << " --- " << "TOF: " << (timeVectorTop[i]-timeVectorTop[i-1])/physical_constants::JULIAN_DAY <<  " --- " << "Delta V: " << deltaVVectorTop[i] << std::endl;
        else std::cout << "Planet: " << i << " --- " << "Time: " << timeVectorTop[i]/physical_constants::JULIAN_DAY << " --- " << "TOF: " << "Launch" <<  " --- " << "Delta V: " << deltaVVectorTop[i] << std::endl;

    }
    std::cout << "=================================================\n=================================================\n";


    //add best for every seed to map
    DVperSeed[seed] = resultingDeltaVTop;

}
    std::string outputFileDVperSeed = tudat_applications::getOutputPath( ) + "MGADSM/IndividualTrajectories/" + flyby_seq;
    input_output::writeDataMapToTextFile(DVperSeed, "DVperSeed.dat", outputFileDVperSeed );

    return 0;

}
