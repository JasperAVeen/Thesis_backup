#ifndef TT_FUNCTIONS_H
#define TT_FUNCTIONS_H

#include <Tudat/InputOutput/basicInputOutput.h>
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

#include "lagrangepoints.h"
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "maketraj.h"


#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLegSettings.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShapingOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

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

//Declaration of all functions:

std::vector< int > createFlyBySeq (std::string);

Eigen::VectorXd createMinPerRad (std::vector< int > );

Eigen::VectorXd createMaxPerRad (std::vector< int > );

std::vector< std::vector< double > > createBounds (std::vector< int >, double, double, int, std::vector< LegTypeTransfer >);

std::vector< std::vector< double > > createCassiniBounds (std::vector< int >, double, double, int, std::vector< std::vector< double > >, std::vector< LegTypeTransfer >);

std::vector< std::vector< double > > createMessengerBounds (std::vector< int >, double, double, int,  std::vector< LegTypeTransfer >);

std::vector< std::vector< double > > createContBounds ( std::pair< double, double >,  std::pair< double, double >, bool, bool);

std::vector< LegTypeTransfer > createLegType (std::string);

std::vector<std::shared_ptr< ephemerides::Ephemeris >> createEphemerisVec (std::vector< int >, std::string = "SUN",  std::string = "ECLIPJ2000", bool = true);

Eigen::VectorXd createGravParamVec (std::vector< int >);

std::vector< std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > > getShapingBasisFunctions( const double , const int );

Eigen::Vector2d getEccentricities ( std::string );

Eigen::Vector2d getSMA ( std::string );

Eigen::Vector6d getDepState( const std::function< Eigen::Vector6d( const double )>, double, std::string, bool, double =  NULL);

Eigen::Vector6d getDepartureState (const std::function< Eigen::Vector6d( const double )>, double, std::string, Eigen::Vector3d );

Eigen::Vector6d getArrivalState (const std::function< Eigen::Vector6d( const double )>, double, std::string, Eigen::Vector3d );

Eigen::Vector6d getArrState( const std::function< Eigen::Vector6d( const double )>, double, std::string, bool, double = NULL );

Eigen::Vector6d getLibrationState (std::string, int, double, bool = false );

Eigen::Vector6d convertSphericalToCartesianState( const Eigen::Vector6d& );

double getDepDeltaV( Eigen::Vector3d, std::string );

double getArrDeltaV( Eigen::Vector3d, std::string );

Eigen::Matrix3d makeSkewed (Eigen::Vector3d);

//Eigen::VectorXd getVariableVector (std::vector< LegTypeTransfer >, std::vector< double >);


//Eigen::Vector3d getLagrangePosVec(std::string, std::string, int);



#endif // TT_FUNCTIONS_H
