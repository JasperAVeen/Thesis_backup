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
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Classes/tt_functions.h"
#include "Classes/maketraj.h"
#include "Classes/trajectoryoptions.h"

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
    spice_interface::loadStandardSpiceKernels();

    //What trajectory should be optimized?
    std::vector<std::string> vecGateway = {"G1", "G2", "G3", "G4", "G5", "G6"};
    std::vector<std::string> vecSegment = {"EG", "GM", "MG", "GE"};
    bool useDCT = true;
    bool useDSM = true;
    bool useMGA = true;
    bool useMGADSM = false; //Always false! For MGADSM simulations use different script!

    //Set Seed:
    unsigned int seed {123}; 

    //What are the bounds wrt departure time and TOF?
    double start_date {10958}; //MJD2000
    double launch_window {5*365.25};   //days
    double end_date {14610}; //MJD2000
    const double TOF_constraint = 1*365.25; //days

    //Toggles:
    const bool printBest = false;
    const bool useTripTime = false;

    const double sunGravitationalParameter = 1.32712428e20;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    for(size_t pointerGW{0}; pointerGW < vecGateway.size(); pointerGW++){
        for(size_t pointerSeg{0}; pointerSeg < vecSegment.size(); pointerSeg++){


    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    std::string gateway = vecGateway.at(pointerGW);
    std::string segment = vecSegment.at(pointerSeg);

    TrajectoryOptions TrajOptions = TrajectoryOptions(gateway, segment, useDCT, useDSM, useMGA, useMGADSM );
    std::vector<std::string> flyby_vec=  TrajOptions.makeStringVector();

    std::map<std::string, Eigen::Vector2d> segmentMap{};

    for(size_t pointerFB{0}; pointerFB < flyby_vec.size(); pointerFB++) {

    std::string flyby_seq = flyby_vec.at(pointerFB);
    std::string pathSuffix = "MGADSM/"+ gateway + "/" + segment + "/" + flyby_seq;

    // Create flyby sequence
    std::vector< int > flybySequence = createFlyBySeq(flyby_seq);

    //Create LegTypeVec
    std::vector< LegTypeTransfer > LegTypeVec =  createLegType (flyby_seq);

    // Create Bounds
    int num_DSMs = std::count(flyby_seq.begin(), flyby_seq.end(), 'd');
    size_t numberOfParameters = flybySequence.size()+(num_DSMs*4);
    std::vector< std::vector< double > > Bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );
    Bounds =  createBounds (flybySequence , start_date, launch_window, num_DSMs, Bounds, LegTypeVec);

    //create Departure and Capture conditions
    Eigen::VectorXd semiMajorAxes = getSMA(flyby_seq);
    Eigen::VectorXd eccentricities = getEccentricities(flyby_seq);

    //Set seed for reproducible results
    pagmo::random_device::set_seed( seed );

    // Create object to compute the problem fitness
    problem prob{ MGAtrajectory( Bounds, flybySequence, LegTypeVec, flyby_seq, TOF_constraint, semiMajorAxes, eccentricities, end_date ) };

    // Select NSGA2 algorithm for problem
    algorithm algo{de(1u, 0.7, 0.97, 2u, 1e-4, 1e-4 )};

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Optimization parameters

    size_t numInd{120};
    int numGen {10};

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////


    // Create an island with individuals
    island isl{algo, prob, numInd };

    // Evolve for numGen generations
    for( int i = 0 ; i < numGen; i++ )
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
//        if(i%1==0){
//        std::cout<<i<<std::endl;
//        std::cout<<"Delta V: "<<isl.get_population().champion_f()[0]<<std::endl;
//        }

    }

    std::cout << "\nCompleted:  " << flyby_seq;


    int numberOfLegs = flybySequence.size( );


    // Create variable vector.
    Eigen::VectorXd TopIndividual(numberOfLegs+1+(4*num_DSMs));

    double TOF = 0;
    for(int i = 0; i < numberOfLegs ; i++){
        TopIndividual[ i ] = isl.get_population().champion_x()[i];
            }

    TopIndividual *= physical_constants::JULIAN_DAY;
    TopIndividual[ numberOfLegs ] = 1;//dummy

    //Write the variable vector for DSM variables below

    if(num_DSMs>0){
        int count {numberOfLegs+1};
        for(size_t i{0}; i<num_DSMs; i++){
            TopIndividual[ count ] = isl.get_population().champion_x()[ count-1 ];
            TopIndividual[ count + 1 ] = isl.get_population().champion_x()[ count ];
            TopIndividual[ count + 2 ] = isl.get_population().champion_x()[ count + 1 ];
            TopIndividual[ count + 3 ] = isl.get_population().champion_x()[ count + 2 ];
            count+=4;
        }
    }

    Eigen::VectorXd minimumPericenterRadii{};
    minimumPericenterRadii.resize( flybySequence.size() );

    if(num_DSMs==0) {minimumPericenterRadii=createMinPerRad(flybySequence);}
    else {
        for (int i{0}; i < numberOfLegs; i++){
            minimumPericenterRadii [i] = TUDAT_NAN;
        }

    }



    makeTraj BestmgaTraj( numberOfLegs, LegTypeVec, flyby_seq, createEphemerisVec (flybySequence),
                          createGravParamVec(flybySequence), TopIndividual, sunGravitationalParameter,
                          minimumPericenterRadii, semiMajorAxes, eccentricities );


    // Vectors for the specific maneuvers and the total delta v
    std::vector< Eigen::Vector3d > positionVectorTop;
    std::vector< double > timeVectorTop;
    std::vector< double > deltaVVectorTop;
    double resultingDeltaVTop;

    // Calculate the orbits
    BestmgaTraj.calculateTrajectory( resultingDeltaVTop );

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

    double BestDeltaV = isl.get_population().champion_f()[0];

    double BestTOF = 0;
    for(int i = 1; i < flybySequence.size() ; i++){
            BestTOF += TopIndividual[i];
    }
    BestTOF = BestTOF/physical_constants::JULIAN_DAY;


    std::cout << "\nGateway: " << gateway << "     Segment: " << segment << "      flyby sequence: " << flyby_seq << std::endl << "Delta V: " << BestDeltaV << "\nTOF: " <<BestTOF;


        Eigen::Vector2d vectorMap;
        vectorMap << BestDeltaV, BestTOF;
        segmentMap[flyby_seq] << vectorMap;

        if(pointerFB == flyby_vec.size()-1){
        std::string outputPathsegmentMap = tudat_applications::getOutputPath( ) + "MGADSM/"+ gateway + "/" + segment;
        input_output::writeDataMapToTextFile(segmentMap, "segmentMap.dat", outputPathsegmentMap );
    }
    
            } //end flybysequence loop
        } // end sequence loop
    } // end GW loop


    return 0;

}
