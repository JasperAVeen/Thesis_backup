#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

#include "src/applyDifferentialCorrection.h"
#include "src/createInitialConditions.h"
#include "src/computeDifferentialCorrection.h"
#include "src/checkEigenvalues.h"
#include "src/computeManifolds.h"
#include "src/connectManifoldsAtTheta.h"
#include "src/propagateOrbit.h"
#include "src/richardsonThirdOrderApproximation.h"
#include "src/stateDerivativeModel.h"

const double SUN_mu { 1.3271244004194e+020};
const double EB_mu {4.0350323550226e+014};
const double Earth_mu {3.986004354361e+014};
const double Moon_mu {4902800066163.8};
const double Mars_mu {42828373620699};


double massParameter = tudat::circular_restricted_three_body_problem::computeMassParameter(SUN_mu , Mars_mu );


int main() {

    //int librationPointNr {1};
    //std::string orbitType {"horizontal"};

    //createInitialConditions(1, "horizontal" );
    //createInitialConditions(2, "horizontal" );

    //createInitialConditions(1, "vertical" );
    //createInitialConditions(2, "vertical" );

    //createInitialConditions(1, "halo" );
    createInitialConditions(1, "halo" );

/*
    int orbitIdOne {5500};
    double desiredJacobiEnergy {3.05};

    Eigen::VectorXd selectedInitialConditions = readInitialConditionsFromFile(librationPointNr, orbitType, orbitIdOne, orbitIdOne + 1, massParameter);
    Eigen::VectorXd refinedJacobiEnergyResult = refineOrbitJacobiEnergy(librationPointNr, orbitType, desiredJacobiEnergy,
                                                                        selectedInitialConditions.segment(1, 6),
                                                                        selectedInitialConditions(0),
                                                                        selectedInitialConditions.segment(8, 6),
                                                                        selectedInitialConditions(7), massParameter);
    Eigen::VectorXd initialStateVector = refinedJacobiEnergyResult.segment(0, 6);
    double orbitalPeriod               = refinedJacobiEnergyResult(6);

    Eigen::MatrixXd fullInitialState = getFullInitialState( initialStateVector );
    std::map< double, Eigen::Vector6d > stateHistory;
    std::pair< Eigen::MatrixXd, double > endState = propagateOrbitToFinalCondition( fullInitialState, massParameter, orbitalPeriod, 1, stateHistory, 100, 0.0 );

    writeStateHistoryToFile( stateHistory, orbitIdOne, orbitType, librationPointNr, 1000, false );
*/

    std::cout << "\n=======\nProgram ended succesfully! \n=======\n";








    return 0;
}
