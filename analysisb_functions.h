#ifndef ANALYSISB_FUNCTIONS_H
#define ANALYSISB_FUNCTIONS_H

#include <Tudat/InputOutput/basicInputOutput.h>
#include <fstream>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <string.h>


#include "lagrangepoints.h"
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "tt_functions.h"
#include "G:/Thesis/applicationOutput.h"
#include "b_maketraj.h"



#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLegSettings.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShapingOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"



using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
using namespace tudat;
using namespace tudat::shape_based_methods;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;


Eigen::Vector6d B_getDepartureState (double, std::string, std::vector< ephemerides::EphemerisPointer >, std::map<double, Eigen::VectorXd>, double, double, Eigen::Vector3d = Eigen::Vector3d::Zero() );

Eigen::Vector6d B_getArrivalState (double, std::string, std::vector< ephemerides::EphemerisPointer >, std::map<double, Eigen::VectorXd>, double, Eigen::Vector3d = Eigen::Vector3d::Zero() );

std::pair< double, double > getTOFbounds(std::string, std::string);

std::vector< int > getParameterIndices (std::vector<std::string>);

//Eigen::VectorXd getVVfromInd (std::vector<std::string>, std::vector<double>, std::vector<int>);

Eigen::VectorXd getNsaveTrajectories (std::vector<std::string>, std::vector<double>, std::vector<int>, std::string, std::unordered_map<int, std::map<double, Eigen::VectorXd>>  );

std::vector< std::vector< double > > B_createBounds (std::string, std::vector<std::string>, double, double, std::pair<double, double>, std::pair<int, int>,  bool returnFlight = false );

std::map<double, Eigen::VectorXd> getLPO_cr3b (std::string, std::string, std::string, int, bool = false);

//std::map<double, Eigen::Vector6d> getLPO_cart (std::string,  std::map<double, Eigen::Vector6d>, double, int = 3);

//double getLPO_period (std::string, std::map<double, Eigen::Vector6d>, double);

//Eigen::Vector6d getLPO_ArrState (std::string, std::map<double, Eigen::Vector6d>, double, double);

//Eigen::Vector6d getLPO_DepState (std::string, std::map<double, Eigen::Vector6d>,  double, double, double);

std::vector< B_LegTypeTransfer > createLegType_B (std::string);

std::unordered_map<int, std::map<double, Eigen::VectorXd>> importLPO_cr3b (std::string system, std::string LPO_type, std::string LP, std::pair<int, int> orbitID_bounds);

Eigen::Vector6d getLPGateway_ArrState (std::string, std::map<double, Eigen::VectorXd>, double, double, bool = true);

Eigen::Vector6d getLPGateway_DepState (std::string, std::map<double, Eigen::VectorXd>, double, double, double, bool = true);

std::map<double,Eigen::Vector6d> getLPGateway_orbit (std::string, std::map<double, Eigen::VectorXd>, double, double, double, bool = true );

std::map<double,Eigen::Vector6d> getBody_orbit (std::string, double, std::string = "Sun", int = 1);

std::unordered_map<int, std::map<double, Eigen::VectorXd>> importLPO_cr3b (std::string system, std::string LPO_type, std::string LP, std::pair<int, int> orbitID_bounds);

std::unordered_map<int, std::map<double, Eigen::VectorXd>> generateCentralOrbits (std::string planet, std::pair<double, double> SMA_bounds, std::pair<double, double> eccentricity_bounds, std::pair<double, double> AOP_bounds, int Nsteps, std::string pathDirectory );

Eigen::Vector6d getCentralGateway_ArrState (std::string planet, bool useRefFrame, std::map<double, Eigen::VectorXd> centralOrbit, double arrivalFraction, double arrivalTime);

Eigen::Vector6d getCentralGateway_DepState (std::string planet, bool useRefFrame, std::map<double, Eigen::VectorXd> centralOrbit, double arrivalFraction, double departurTime, double stayTime);

std::map<double,Eigen::Vector6d> getCentralGateway_orbit (std::string planet, std::map<double, Eigen::VectorXd> centralOrbit, double arrivalFraction, double arrivalTime, double stayTime, bool useRefFrame);

Eigen::Vector6d convertXdtoSixd (Eigen::VectorXd Xd_vector);

std::unordered_map<int, std::map<double, Eigen::VectorXd>> createOrbitLibrary (std::string Gateway, std::map<std::string, std::pair<double, double>> boundsMap, std::string LPO_type, int Nsteps, std::string pathDirectory);

std::pair<double, double> getOrbitIDbounds( int, std::string, std::string, int = 0, int = 0);

std::vector<std::string> getTrajectoryOptions (std::string Gateway, std::string segment);



typedef std::vector< std::shared_ptr< tudat::shape_based_methods::BaseFunctionHodographicShaping > > BaseFunctionVector;


#endif // ANALYSISB_FUNCTIONS_H
