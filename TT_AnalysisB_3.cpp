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
#include <cmath>
#include <Eigen/Core>


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
#include "Classes/analysisb.h"
#include "Classes/analysisb_functions.h"


#include <Tudat/InputOutput/basicInputOutput.h>
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

//! Execute  main
int main( )
{   //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////Toggles and buttons//////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    //Optimize using two objectives? (DV & TOF)
    bool multiObj = false;


    std::vector<std::string> Gateway_options = {"G1"};//, "G7", "G8", "G9"};
    std::string LPO_type = "halo";

    for(size_t GW_iter{0}; GW_iter<Gateway_options.size(); GW_iter++){
        std::string Gateway = Gateway_options.at(GW_iter);


    std::vector<std::string> EG_trajectoryOptions = getTrajectoryOptions(Gateway, "EG");//{"EdG1"};//, "EG1", "EdG1", "EmG1", "EdmdG1"}; //folows from Analysis A
    std::vector<std::string> GM_trajectoryOptions = getTrajectoryOptions(Gateway, "GM");//{"G1dmdM"};//, "G7EmM", "G7EM", "G7mM"};//, "G1M", "G1EmM", "G1mM"}; //folows from Analysis A

    //What are the bounds wrt departure time and TOF?
    double start_date {10958}; //MJD2000
    double launch_window {5*365.25};   //days

    //Min and max stay time at gateway
    std::pair <double, double> Staytime = std::make_pair( 1, 50 ); //MJD

    //create Bounds Map
    int Nsteps = 3;
    std::map<std::string, std::pair<double, double>> boundsMap {};
    boundsMap["SMA"] = std::make_pair(300e3, 50000e3);
    boundsMap["Ecc"] = std::make_pair(0, 0.90);
    boundsMap["AOP"] = std::make_pair(0, 360);
    boundsMap["orbitID"] = getOrbitIDbounds(Nsteps, Gateway, LPO_type);


    //Set seed vector:
    std::vector<unsigned int> seedvec = {123, 234, 345, 456, 567};//, 345, 456, 567};


    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    spice_interface::loadStandardSpiceKernels();

    std::string pathDirectory = "AnalysisB_1805/"+ Gateway + "/" + LPO_type + "/";

    //Import orbits of interest
    std::unordered_map<int, std::map<double, Eigen::VectorXd>> orbitLibrary = createOrbitLibrary(Gateway, boundsMap, LPO_type, Nsteps, pathDirectory);

    //Map to save all results in for each EM_trajectory
    std::map<std::string, Eigen::VectorXd> resultsMap;

    for (size_t EG_iter{0}; EG_iter < EG_trajectoryOptions.size(); EG_iter++) {
        for (size_t GM_iter{0}; GM_iter < GM_trajectoryOptions.size(); GM_iter++) {

            //Define entire trajectory
            std::vector<std::string> EM_trajectory = {EG_trajectoryOptions.at(EG_iter), GM_trajectoryOptions.at(GM_iter) };
            std::string EMtraj_string = EM_trajectory.at(0) + " x " + EM_trajectory.at(1);
            std::string pathSuffix = pathDirectory + EMtraj_string;

            std::cout << "\n********* - Starting optimization of: " << EM_trajectory[0] << "  x  " << EM_trajectory[1] << " - *********\n\n";

            //Create Bounds
            std::vector< std::vector< double > > Bounds =  B_createBounds (Gateway, EM_trajectory, start_date, launch_window, Staytime, boundsMap["orbitID"] );
            std::vector<int> ParameterIndices = getParameterIndices (EM_trajectory);
            size_t numberOfParameters = ParameterIndices.at(0);

            //creating maps to fill during optimization
            std::map<unsigned int, double > SeedmapFitness{}; //Map to save Objectives mapped to seed number
            std::map<double, std::vector<double>> FitnessmapIndividual{}; //Map to save Variable Vector to Objective
            std::vector<double> seedFitnessVec{};

             std::cout << "\n-Bounds were created!";



            //********************************************OPTIMIZATION*******************************************************************************************************************************************************************************
                       //Loop over different seeds
                       for (int pointerSeed{0}; pointerSeed < seedvec.size(); pointerSeed++){ //start seed loop
                       unsigned int seed {seedvec.at(pointerSeed)};
                       std::string seedstring {std::to_string(seed)};

                       //Set seed for reproducible results
                       pagmo::random_device::set_seed( seed );

                       // Create object to compute the problem fitness
                       problem prob{ AnalysisB( EM_trajectory, Bounds, orbitLibrary, multiObj ) };

                       // Select optim algorithm
                       algorithm algo;

                       if(multiObj){
                           algo = nsga2();
                       }
                       else{
                           algo = de(1u, 0.7, 0.97, 2u); // Taken from Musegaas
                       }

                       std::cout << "\n-Optimization problem setup!";


                       //////////////////////////////////////////////////////////////////////////////////////////////
                       //////////////////////////////////////////////////////////////////////////////////////////////
                       // Optimization parameters
                       size_t numInd {10 + 2 * numberOfParameters}; // Taken from Musegaas
                       int numGen {5000};
                       int saveIter {500};

                       //////////////////////////////////////////////////////////////////////////////////////////////
                       //////////////////////////////////////////////////////////////////////////////////////////////

                       // Create an island with individuals
                       island isl{algo, prob, numInd };

                       // Evolve for numGen generations
                       for( int i = 0 ; i < numGen+1; i++ )
                           {

                           isl.evolve( );

                           while( isl.status( ) != pagmo::evolve_status::idle &&
                               isl.status( ) != pagmo::evolve_status::idle_error )
                           {
                               isl.wait( );
                           }



                           //isl.wait_check( ); // Raises errors

                           //std::cout << "\nchecked!";


                           // Write current iteration results to file
                           if(i%saveIter==0){
                           if(!multiObj){
                           printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga" + std::to_string( i ), false, pathSuffix+"/SingleObj/"+seedstring );
                           printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga" + std::to_string( i ), true, pathSuffix+"/SingleObj/"+seedstring );
                           std::cout << "\n============ Iter: " << i << " - finished" << "   -   Fitness: "<<isl.get_population().champion_f()[0]<<std::endl;
                             }
                           if(multiObj){
                           printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga" + std::to_string( i ), false, pathSuffix+"/MultiObj/"+seedstring );
                           printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga" + std::to_string( i ), true, pathSuffix+"/MultiObj/"+seedstring );
                           std::cout << "\n============ Iter: " << i << " - finished" << std::endl;
                           }
                           }
                       }
           //********************************************************************************************************************************************************************************************************************************


                       if(!multiObj){

                       //////////////////////////////////////////////////////////////////////////////////////////////
                       ///////////////////////////////Save parameters for each seed ////////////////////////////////
                       double fitness = isl.get_population().champion_f()[0];
                       std::vector<double> Individual = isl.get_population().champion_x();

                       SeedmapFitness[seed] = fitness;
                       FitnessmapIndividual[fitness] = Individual;
                       seedFitnessVec.push_back(fitness);

                       std::cout << "\n***************\nGateway: " << Gateway << "     Transfer: " << EM_trajectory[0] << "  x  " << EM_trajectory[1] << "      Seed: " << seed <<   "\nFitness: " << fitness << "\n***************\n\n";
                       }
                       else {
                           std::cout << "\n***************\nGateway: " << Gateway << "     Transfer: " << EM_trajectory[0] << "  x  " << EM_trajectory[1] << "      Seed: " << seed  << "\n***************\n\n";

                       }

                        } //end of seed loop


                    if(!multiObj){

                       // Get VV of the seed number with lowest DV outcome
                       double minFitness_seeds = *std::min_element(seedFitnessVec.begin(),seedFitnessVec.end());
                       std::vector<double> BestIndividual =   FitnessmapIndividual[ minFitness_seeds ];  //BestIndividual over all seed values

                       std::cout << "\nBestIndividual:   ";
                       for(int i{0}; i < BestIndividual.size(); i++) {
                           std::cout << BestIndividual.at(i) << ", ";
                       }

                       //Compute all parameters and save all files for this BestIndividual
                       Eigen::VectorXd BestIndividual_Results = getNsaveTrajectories ( EM_trajectory,  BestIndividual,  ParameterIndices,  pathSuffix, orbitLibrary );

                       std::cout << "\nBestIndividual_Results: \n" << BestIndividual_Results;

                        resultsMap [EMtraj_string] = BestIndividual_Results;
                        std::string outputPath_resultsMap = tudat_applications::getOutputPath( ) + "AnalysisB/"+ Gateway + "/" + LPO_type;
                        input_output::writeDataMapToTextFile(resultsMap, "resultsMap.dat", outputPath_resultsMap );
                    }

                    std::cout << "\n********* - Finished optimization of: " << EM_trajectory[0] << "  x  " << EM_trajectory[1] << " - *********\n\n";


        } //end of GM options loop
    } // end of EG options loop

    std::cout << "\n*********\n*********\n*********\n - Finished optimization of Gateway:  " << Gateway << " - \n*********\n*********\n*********\n";


    } //end of Gateway loop
    std::cout << "\n\n------------------------\nProgram ended succesfully!\n--------------------------\n\n";

    return 0;

}
