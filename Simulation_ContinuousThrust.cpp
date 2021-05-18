#include <iostream>
#include <fstream>
#include <functional>

#include <boost/filesystem.hpp>
#include "G:/Thesis/applicationOutput.h"

#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Classes/tt_functions.h"
#include "Classes/hodographicshapingproblem.h"
#include "Classes/trajectoryoptions.h"



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
    spice_interface::loadStandardSpiceKernels();

    //What trajectory should be optimized?
    std::vector<std::string> vecGateway = { "G1", "G2", "G3", "G4", "G5", "G6"};
    std::vector<std::string> vecSegment = { "EG", "GM"};
    bool useDCT = true;
    bool useDSM = false;
    bool useMGA = false;
    bool useMGADSM = false;

    //What are the bounds wrt departure time and TOF?
    std::pair< double, double > departureTimeBounds =
            std::make_pair( 10958 * physical_constants::JULIAN_DAY, (10958 + 5 * 365.25) * physical_constants::JULIAN_DAY  );

    //What is the initial mass of the spacecraft?
    double initialMass = 1000.0;

    //What is the specific Impulse of the spacecraft's engines?
    double specificImpulse = 3000.0;

    //What is the maximum number of allowed revolutions around the central body?
    int maxRev = 0;

    //Set seed vector:
    std::vector<unsigned int> seedvec = {123, 456, 789, 987, 654}; //terugveranderen

    // MultiObjecttive? Yes -- 2, No -- 1
    int MultiObjective {1};

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    for(size_t pointerGW{0}; pointerGW < vecGateway.size(); pointerGW++){
        for(size_t pointerSeg{0}; pointerSeg < vecSegment.size(); pointerSeg++){


    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    std::string gateway = vecGateway.at(pointerGW);
    std::string segment = vecSegment.at(pointerSeg);

    std::pair< double, double > timeOfFlightBounds =
                    std::make_pair( 1 * physical_constants::JULIAN_DAY, 1.0 * physical_constants::JULIAN_YEAR );

    if((pointerGW < 4) && (pointerSeg == 0 || pointerSeg == 3)) //Em gateways- segments EG and GE
    {  timeOfFlightBounds = std::make_pair( 1.0 * physical_constants::JULIAN_DAY, 25.0 * physical_constants::JULIAN_DAY );}
    if((pointerGW == 4) && (pointerSeg == 0 || pointerSeg == 3)) //SE gateways- segments EG and GE
    { timeOfFlightBounds = std::make_pair( 1 * physical_constants::JULIAN_DAY, 100.0 * physical_constants::JULIAN_DAY );}
    if((pointerGW == 5) && (pointerSeg == 1 || pointerSeg == 2)) //SM gateways- segments GM and MG
    {timeOfFlightBounds = std::make_pair( 1 * physical_constants::JULIAN_DAY, 100.0 * physical_constants::JULIAN_DAY );}



    TrajectoryOptions TrajOptions = TrajectoryOptions(gateway, segment, useDCT, useDSM, useMGA, useMGADSM );
    std::vector<std::string> flyby_vec=  TrajOptions.makeStringVector();

    for(size_t pointerFB{0}; pointerFB < flyby_vec.size(); pointerFB++) {

    std::string trajectory = flyby_vec.at(pointerFB);

    std::string pathSuffix = "ContinuousThrust_0512/"+ gateway + "/" + segment + "/";

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

    //Do we need to use departure and arrival orbits?
    bool useDepOrbit = false;
    bool useArrOrbit = false;

    if(pointerSeg == 0 || pointerSeg == 2){useDepOrbit = true; useArrOrbit = false;}
    if(pointerSeg == 1 || pointerSeg == 3){useDepOrbit = false; useArrOrbit = true;}

    std::vector< std::vector< double > > bounds{};

   bounds = createContBounds ( departureTimeBounds, timeOfFlightBounds, useDepOrbit, useArrOrbit);

   size_t numParams = bounds[0].size();

   for (int revolutions = 0; revolutions < maxRev+1; revolutions++)
   {

    for( int useMultiObjective = 0; useMultiObjective < MultiObjective; useMultiObjective++ )
   {
        std::cout << "\n----------------------------------------------------------------------------\n Gateway:  " <<  gateway << "  Segment:  " << segment << "  Trajectory:  " << trajectory << std::endl;


        ////////////////////////Optimization variables/////////////////////////////
        size_t Population_size = 10 * numParams;
        size_t num_Gen = 500 * numParams;
        int output_iter = 100;
        //////////////////////////////////////////////////////////////////////////

        std::map<int, Eigen::Vector3d> seedResults;


        //Loop over different seeds
        for (int pointerSeed{0}; pointerSeed < seedvec.size(); pointerSeed++){ //start seed loop
        unsigned int seed {seedvec.at(pointerSeed)};
        std::string seedstring {std::to_string(seed)};

        //Set seed for reproducible results
        pagmo::random_device::set_seed( seed );


            // Create object to compute the problem fitness
            problem prob{ HodographicShapingProblem( departureStateFunction, arrivalStateFunction, spice_interface::getBodyGravitationalParameter( "Sun" ),
                            revolutions, std::bind( &getShapingBasisFunctions, std::placeholders::_1, revolutions ), bounds, trajectory, useDepOrbit, useArrOrbit,
                            useMultiObjective, initialMass ) };

            //sade, gaco, sga, de
            algorithm algo;
            if( !useMultiObjective )
            {
               algo = algorithm {de(1u, 0.7, 0.97, 2u)}; // Taken from Musegaas
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
                        Eigen::Vector3d Vinf_vec; Vinf_vec << isl.get_population( ).champion_x( ).at( 7 ) , isl.get_population( ).champion_x( ).at( 8 ) , isl.get_population( ).champion_x( ).at( 9 );
                        //std::cout<<"Iteration: "<<" "<<i<<"; Best Total Delta V: "<<isl.get_population( ).champion_f( ).at( 0 ) <<std::endl; if(useDepOrbit) std::cout << "; Vinf_vec norm: " << Vinf_vec.norm() << "\n DepDeltaV" << getDepDeltaV(Vinf_vec, trajectory) << std::endl;

                        printPopulationToFile( isl.get_population( ).get_x( ), "hodograph_single_objective_" + std::to_string( i / output_iter ), false, pathSuffix + seedstring + "/" );
                        printPopulationToFile( isl.get_population( ).get_f( ), "hodograph_single_objective_" + std::to_string( i / output_iter ), true, pathSuffix + seedstring + "/" );
                    }
                    else
                    {
                        std::cout<<"Iteration: "<<" "<<i<<std::endl;

                        printPopulationToFile( isl.get_population( ).get_x( ), "hodograph_multi_objective_" + std::to_string( i / output_iter ), false, pathSuffix  + seedstring + "/");
                        printPopulationToFile( isl.get_population( ).get_f( ), "hodograph_multi_objective_" + std::to_string( i / output_iter ), true, pathSuffix  + seedstring + "/");
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
                Eigen::Vector3d Vinf_vec_dep = Eigen::Vector3d::Zero();
                Eigen::Vector3d Vinf_vec_arr = Eigen::Vector3d::Zero();

                if(useDepOrbit && !useArrOrbit) { Vinf_vec_dep <<  bestPopulation.at( 7 ),  bestPopulation.at( 8 ),  bestPopulation.at( 9 );}
                if(!useDepOrbit && useArrOrbit) { Vinf_vec_arr <<  bestPopulation.at( 7 ),  bestPopulation.at( 8 ),  bestPopulation.at( 9 );}
                if(useDepOrbit && useArrOrbit) { Vinf_vec_dep <<  bestPopulation.at( 7 ),  bestPopulation.at( 8 ),  bestPopulation.at( 9 );
                                                 Vinf_vec_arr <<  bestPopulation.at( 10 ),  bestPopulation.at( 11 ),  bestPopulation.at( 12 );}

                double bestDeltaV = isl.get_population( ).champion_f( ).at( 0 );

                Eigen::Vector2d FlightTimes; FlightTimes << timeOfFlight/physical_constants::JULIAN_DAY, departureTime/physical_constants::JULIAN_DAY;

                std::map<double, Eigen::Vector2d> TrajectoryResults;
                TrajectoryResults[bestDeltaV] << FlightTimes;
                input_output::writeDataMapToTextFile(
                            TrajectoryResults, "TrajectoryResults.dat", tudat_applications::getOutputPath( )+pathSuffix   + seedstring + "/" );

                //Display most important output
                std::cout << " -- Seed: " << seed << "   -   Departure Time:  " << departureTime/physical_constants::JULIAN_DAY << "   -   TOF:  " << timeOfFlight/physical_constants::JULIAN_DAY << "   -   DeltaV:  " <<bestDeltaV << "\n";


                std::vector< std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > > shapingFunctions = getShapingBasisFunctions(
                            timeOfFlight, revolutions );

                Eigen::Vector6d departureState = getDepartureState( departureStateFunction, departureTime, trajectory, Vinf_vec_dep );
                Eigen::Vector6d arrivalState = getArrivalState( arrivalStateFunction, arrivalTime, trajectory, Vinf_vec_arr );

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
                startendconditions[3] = departureState;
                startendconditions[4] = arrivalState;
                input_output::writeDataMapToTextFile(
                            startendconditions, "startendconditions.dat", tudat_applications::getOutputPath( )+pathSuffix   + seedstring + "/");

                int numberOfSteps = timeOfFlight/1000;
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
                            hodographicShapingTrajectory, "hodographicShapingOptimalTrajectory.dat", tudat_applications::getOutputPath( )+pathSuffix  + seedstring + "/");

                input_output::writeDataMapToTextFile(
                            hodographicShapingMassProfile, "hodographicShapingOptimalMassProfile.dat", tudat_applications::getOutputPath( )+pathSuffix  + seedstring + "/");

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustProfile, "hodographicShapingOptimalThrustProfile.dat", tudat_applications::getOutputPath( )+pathSuffix  + seedstring + "/");

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustAcceleration, "hodographicShapingOptimalThrustAcceleration.dat", tudat_applications::getOutputPath( )+pathSuffix  + seedstring + "/");


                seedResults [seed] << bestDeltaV, timeOfFlight/physical_constants::JULIAN_DAY, departureTime/physical_constants::JULIAN_DAY;

                }

            input_output::writeDataMapToTextFile(
                        seedResults, "seedResults.dat", tudat_applications::getOutputPath( )+pathSuffix );

            } //end of seed loop


        }
   }
}
        }
    }
}
