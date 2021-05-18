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
#include "Problems/applicationOutput.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/exportTrajectory.h"
#include "Classes/tt_functions.h"

#include <Tudat/InputOutput/basicInputOutput.h>
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/planetTrajectory.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

//Toggles



//! Execute  main
int main( )
{
    const double sunGravitationalParameter = 1.32712428e20;

    // Calculate trajectories of the planets and output to file
    std::vector< Eigen::Vector3d > positionVectorEarth;
    std::vector< double > timeVectorEarth;
    std::vector< Eigen::Vector3d > positionVectorVenus;
    std::vector< double > timeVectorVenus;
    std::vector< Eigen::Vector3d > positionVectorMercury;
    std::vector< double > timeVectorMercury;
    std::vector< Eigen::Vector3d > positionVectorMars;
    std::vector< double > timeVectorMars;
    std::vector< Eigen::Vector3d > positionVectorJupiter;
    std::vector< double > timeVectorJupiter;
    std::vector< Eigen::Vector3d > positionVectorSaturn;
    std::vector< double > timeVectorSaturn;
    std::vector< Eigen::Vector3d > positionVectorUranus;
    std::vector< double > timeVectorUranus;
    std::vector< Eigen::Vector3d > positionVectorNeptune;
    std::vector< double > timeVectorNeptune;

    // Earth
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorEarth,
                timeVectorEarth );

    // Venus
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorVenus,
                timeVectorVenus );

    // Mercury
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorMercury,
                timeVectorMercury );
    // Mars
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorMars,
                timeVectorMars );
    // Jupiter
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorJupiter,
                timeVectorJupiter );
    // Saturn
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorSaturn,
                timeVectorSaturn );
    // Uranus
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::uranus ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorUranus,
                timeVectorUranus );
    // Neptune
    returnSingleRevolutionPlanetTrajectory(
                std::make_shared< ephemerides::ApproximatePlanetPositions >(
                    ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::neptune ),
                sunGravitationalParameter,
                1171.64503236,
                1000.0,
                positionVectorNeptune,
                timeVectorNeptune );



    std::string outputMap = "PlanetOrbits";
    std::string outputFilePlanetE = "C:/tudatBundle/tudatApplications/Thesis/SimulationOutput/PlanetOrbits/earthTrajectory.dat";
    std::cout << outputFilePlanetE;
    writeTrajectoryToFile( positionVectorEarth, timeVectorEarth, outputFilePlanetE );

    std::string outputFilePlanetV = tudat_applications::getOutputPath(  ) + outputMap+ "/venusTrajectory.dat";
    writeTrajectoryToFile( positionVectorVenus, timeVectorVenus, outputFilePlanetV );

    std::string outputFilePlanetM = tudat_applications::getOutputPath(  )+ outputMap + "/mercuryTrajectory.dat";
    writeTrajectoryToFile( positionVectorMercury, timeVectorMercury, outputFilePlanetM );

    std::string outputFilePlanetMars = tudat_applications::getOutputPath(  )+ outputMap + "/marsTrajectory.dat";
    writeTrajectoryToFile( positionVectorMars, timeVectorMars, outputFilePlanetMars );

    std::string outputFilePlanetJ = tudat_applications::getOutputPath(  )+ outputMap + "/jupiterTrajectory.dat";
    writeTrajectoryToFile( positionVectorJupiter, timeVectorJupiter, outputFilePlanetJ );

    std::string outputFilePlanetS = tudat_applications::getOutputPath(  )+ outputMap + "/saturnTrajectory.dat";
    writeTrajectoryToFile( positionVectorSaturn, timeVectorSaturn, outputFilePlanetS );

    std::string outputFilePlanetU = tudat_applications::getOutputPath(  )+ outputMap + "/uranusTrajectory.dat";
    writeTrajectoryToFile( positionVectorUranus, timeVectorUranus, outputFilePlanetU );

    std::string outputFilePlanetN = tudat_applications::getOutputPath(  )+ outputMap + "/neptuneTrajectory.dat";
    writeTrajectoryToFile( positionVectorNeptune, timeVectorNeptune, outputFilePlanetN );


}
