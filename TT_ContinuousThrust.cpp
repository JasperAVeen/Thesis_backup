#include <iostream>
#include <fstream>
#include <functional>

#include <boost/filesystem.hpp>
#include "applicationOutput.h"

#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Classes/tt_functions.h"
#include "Classes/hodographicshapingproblem.h"


#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/sphericalShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLegSettings.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include "Tudat/SimulationSetup/optimisationSettings.h"

using namespace tudat;
using namespace tudat::shape_based_methods;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;



//! Execute  main
int main( )
{   //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////Toggles and buttons//////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
   //What trajectory should be optimized?
    std::string trajectory {"EM"};

    //What are the bounds wrt departure time and TOF?
    std::pair< double, double > departureTimeBounds =
            std::make_pair( 10958 * physical_constants::JULIAN_DAY, (10958 + 5 * 365.25) * physical_constants::JULIAN_DAY  );
    std::pair< double, double > timeOfFlightBounds =
            std::make_pair( 100.0 * physical_constants::JULIAN_DAY, physical_constants::JULIAN_YEAR );

    //What is the initial mass of the spacecraft?
    double initialMass = 1000.0;

    //What is the specific Impulse of the spacecraft's engines?
    double specificImpulse = 3000.0;

    //What is the maximum number of allowed revolutions around the central body?
    int maxRev = 0;

    //Do we incorporate orbits into our boundary conditions (departure & capture)?
    bool useDepOrbit = false;
    bool useArrOrbit = false;

    //Set seed:
    pagmo::random_device::set_seed( 123456789 );

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    spice_interface::loadStandardSpiceKernels();

    std::vector< int > flybySequence = createFlyBySeq(trajectory);

    // Ephemeris functions of bodies.
    std::vector< ephemerides::EphemerisPointer > ephemerisVector = createEphemerisVec ( flybySequence);

    std::function< Eigen::Vector6d( const double ) > departureStateFunction = [ = ]( const double currentTime )
    { return ephemerisVector[0]->getCartesianState( currentTime ); };
    std::function< Eigen::Vector6d( const double ) > arrivalStateFunction = [ = ]( const double currentTime )
    { return ephemerisVector[1]->getCartesianState( currentTime ); };



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////        GRID SEARCH FOR HODOGRAPHIC SHAPING LOWEST-ORDER SOLUTION            /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::vector< double > > bounds{};

   bounds = createContBounds ( departureTimeBounds, timeOfFlightBounds, useDepOrbit, useArrOrbit);

   for (int revolutions = 0; revolutions < maxRev+1; revolutions++)
   {
       std::cout << "==============================\nRevolutions: " << revolutions << "\n==============================\n\n";

    for( int useMultiObjective = 0; useMultiObjective < 2; useMultiObjective++ )
   {

        ////////////////////////Optimization variables/////////////////////////////
        size_t Population_size = 3000;
        int num_Gen = 300;
        int output_iter = 50;
        //////////////////////////////////////////////////////////////////////////


            // Create object to compute the problem fitness
            problem prob{ HodographicShapingProblem( departureStateFunction, arrivalStateFunction, spice_interface::getBodyGravitationalParameter( "Sun" ),
                            revolutions, std::bind( &getShapingBasisFunctions, std::placeholders::_1, revolutions ), bounds, trajectory, useDepOrbit, useArrOrbit,
                            useMultiObjective, initialMass ) };

            //sade, gaco, sga, de
            algorithm algo;
            if( !useMultiObjective )
            {
                algo = algorithm{ de(1u, 0.7, 0.97, 2u, 1e-4, 1e-4 )};
            }
            else
            {
                algo = algorithm{ nsga2( ) };
            }

            // Create an island with 1000 individuals
            island isl{algo, prob, Population_size };

            // Evolve for 100 generations
            for( int i = 0 ; i < num_Gen+1; i++ )
            {
                isl.evolve( );
                while( isl.status( ) != pagmo::evolve_status::idle &&
                       isl.status( ) != pagmo::evolve_status::idle_error )
                {
                    isl.wait( );
                }

                if( i % output_iter  == 0 )
                {
                    if( !useMultiObjective )
                    {
                        std::cout<<"Iteration: "<<" "<<i<<"; Best Delta V: "<<isl.get_population( ).champion_f( ).at( 0 ); if(useDepOrbit) std::cout << "; SMA: " << isl.get_population( ).champion_x( ).at( 7 ); std::cout << std::endl;

                        printPopulationToFile( isl.get_population( ).get_x( ), "hodograph_single_objective_" + std::to_string( i / output_iter ), false, "ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );
                        printPopulationToFile( isl.get_population( ).get_f( ), "hodograph_single_objective_" + std::to_string( i / output_iter ), true, "ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );
                    }
                    else
                    {
                        std::cout<<"Iteration: "<<" "<<i<<std::endl;

                        printPopulationToFile( isl.get_population( ).get_x( ), "hodograph_multi_objective_" + std::to_string( i / output_iter ), false, "ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );
                        printPopulationToFile( isl.get_population( ).get_f( ), "hodograph_multi_objective_" + std::to_string( i / output_iter ), true, "ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );
                    }
                }

            }

            if( !useMultiObjective )
            {
                //Create bestPopulation
                std::vector< double > bestPopulation = isl.get_population( ).champion_x( );

                Eigen::VectorXd radialFreeParameters = Eigen::VectorXd::Zero( 2 );
                Eigen::VectorXd normalFreeParameters = Eigen::VectorXd::Zero( 0 );
                Eigen::VectorXd axialFreeParameters = Eigen::VectorXd::Zero( 3 );

                radialFreeParameters << bestPopulation.at( 2 ), bestPopulation.at( 3 );
                axialFreeParameters << bestPopulation.at( 4 ), bestPopulation.at( 5 ), bestPopulation.at( 6 );

                double timeOfFlight = bestPopulation.at( 1 );
                double departureTime = bestPopulation.at( 0 );
                double arrivalTime = departureTime + timeOfFlight;
                double smadep=0;
                double smaarr=0;
                if(useDepOrbit) smadep = bestPopulation.at(7);
                if(useArrOrbit && !useDepOrbit) smaarr =  bestPopulation.at(7);
                if(useArrOrbit && useDepOrbit) smaarr =  bestPopulation.at(8);


                double bestDeltaV = isl.get_population( ).champion_f( ).at( 0 );

                Eigen::Vector2d FlightTimes; FlightTimes << departureTime/physical_constants::JULIAN_DAY, timeOfFlight/physical_constants::JULIAN_DAY;

                std::map<double, Eigen::Vector2d> TrajectoryResults;
                TrajectoryResults[bestDeltaV] << FlightTimes;
                input_output::writeDataMapToTextFile(
                            TrajectoryResults, "TrajectoryResults.dat", tudat_pagmo_applications::getOutputPath( )+"ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );


                std::vector< std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > > shapingFunctions = getShapingBasisFunctions(
                            timeOfFlight, revolutions );

                Eigen::Vector6d departureState = getDepState( departureStateFunction, departureTime, trajectory, useDepOrbit, smadep );
                Eigen::Vector6d arrivalState = getArrState( arrivalStateFunction, arrivalTime, trajectory, useArrOrbit, smaarr );
                std::shared_ptr< HodographicShaping > hodographicShaping =
                        std::make_shared< HodographicShaping >(
                            departureState, arrivalState, timeOfFlight,
                            spice_interface::getBodyGravitationalParameter( "Sun" ), revolutions,
                            shapingFunctions.at( 0 ), shapingFunctions.at( 1 ), shapingFunctions.at( 2 ),
                            radialFreeParameters, normalFreeParameters, axialFreeParameters, initialMass );


                // Save results
                std::map<int, Eigen::Vector6d> startendconditions;
                startendconditions[1] = departureStateFunction(departureTime);
                startendconditions[2] = arrivalStateFunction(arrivalTime);
                startendconditions[3] = getDepState( departureStateFunction, departureTime, trajectory, useDepOrbit );
                startendconditions[4] = getArrState( arrivalStateFunction, arrivalTime, trajectory, useArrOrbit );
                input_output::writeDataMapToTextFile(
                            startendconditions, "startendconditions.dat", tudat_pagmo_applications::getOutputPath( )+"ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );

                int numberOfSteps = 1000;
                double stepSize = timeOfFlight / static_cast< double >( numberOfSteps );
                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                        std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize );

                std::vector< double > epochsToSaveResults;
                for ( int i = 0 ; i <= numberOfSteps ; i++ )
                {
                    epochsToSaveResults.push_back( i * stepSize );
                }

                std::map< double, Eigen::Vector6d > hodographicShapingTrajectory;
                std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
                std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
                std::map< double, Eigen::VectorXd > hodographicShapingThrustAcceleration;

                hodographicShaping->getTrajectory(
                            epochsToSaveResults, hodographicShapingTrajectory );

                hodographicShaping->getMassProfile(
                            epochsToSaveResults, hodographicShapingMassProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
                hodographicShaping->getThrustForceProfile(
                            epochsToSaveResults, hodographicShapingThrustProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
                hodographicShaping->getCylindricalThrustAccelerationProfile(
                            epochsToSaveResults, hodographicShapingThrustAcceleration );

                input_output::writeDataMapToTextFile(
                            hodographicShapingTrajectory, "hodographicShapingOptimalTrajectory.dat", tudat_pagmo_applications::getOutputPath( )+"ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );

                input_output::writeDataMapToTextFile(
                            hodographicShapingMassProfile, "hodographicShapingOptimalMassProfile.dat", tudat_pagmo_applications::getOutputPath( )+"ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustProfile, "hodographicShapingOptimalThrustProfile.dat", tudat_pagmo_applications::getOutputPath( )+"ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustAcceleration, "hodographicShapingOptimalThrustAcceleration.dat", tudat_pagmo_applications::getOutputPath( )+"ContinuousThrust/"+trajectory+ "/" +  std::to_string(revolutions) );
            }
        }
   }
}
