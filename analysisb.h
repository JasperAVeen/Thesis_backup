#ifndef ANALYSISB_H
#define ANALYSISB_H


#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/InputOutput/basicInputOutput.h>

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"

#include <random>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"
#include <pagmo/rng.hpp>
#include "tt_functions.h"
#include "maketraj.h"
#include "analysisb_functions.h"
#include "b_maketraj.h"




#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
using namespace tudat;
using namespace pagmo;

struct AnalysisB
{
    AnalysisB( ) {};


    AnalysisB( std::vector<std::string>, std::vector< std::vector< double > >, std::unordered_map<int, std::map<double, Eigen::VectorXd>>, bool );

    // Calculates the fitness
    vector_double fitness(const vector_double &) const;

    std::pair<vector_double, vector_double> get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        if(multiObj_ )
        {
            return 2u;
        }
        else
        {
            return 1u;
        }
    }

    typedef std::vector< std::shared_ptr< tudat::shape_based_methods::BaseFunctionHodographicShaping > > BaseFunctionVector;


private:

    //general parameters
    mutable double totalFitnessValue1{0}; //DV
    mutable double totalFitnessValue2{0}; //TOF in hours

    bool comments_;
    bool multiObj_;
    std::string EG_trajectory_;
    std::string GM_trajectory_;
    bool EG_isCont_;
    bool GM_isCont_;
    std::vector< std::vector< double > > problemBounds_;
    int numberOfRevolutions_{0};
    double initialMass_{1000.0};
    const double sunGravitationalParameter_ = 1.32712428e20;
    std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction_ = std::bind( &getShapingBasisFunctions, std::placeholders::_1, numberOfRevolutions_ );
    mutable double EG_fitnessvalue_ {0};
    mutable double EG_trajectoryDV_ {0};
    mutable double EG_escapeDV_ {0};
    mutable double GM_fitnessvalue_ {0};
    mutable double GM_trajectoryDV_ {0};
    mutable double GM_captureDV_ {0};
    mutable Eigen::Vector3d Vinf_vec_dep_ = Eigen::Vector3d::Zero();
    mutable Eigen::Vector3d Vinf_vec_arr_ = Eigen::Vector3d::Zero();
    int V_inf_vars_ = 3;
    double departureDeltaV_ = 0;
    double arrivalDeltaV_ = 0;

    //Gateway orbit parameters
    mutable std::unordered_map<int, std::map<double, Eigen::VectorXd>> orbitLibrary_;
    mutable std::map<double, Eigen::VectorXd> specificOrbit_;

    mutable int orbitID_;
    mutable double orbitFraction_;
    mutable double stayTime_{0};
    mutable double stayTimeSec_{0};



    //EG parameters
    std::vector< B_LegTypeTransfer > EG_LegTypeVector_;
    std::vector< int > EG_flybysequence_;
    int EG_numLegs_;
    int EG_numDSM_;
    std::vector< ephemerides::EphemerisPointer > EG_ephemerisVector_;
    Eigen::VectorXd EG_gravitationalParameterVector_;
    Eigen::Vector2d EG_semiMajorAxes_;
    Eigen::Vector2d EG_eccentricities_;
    Eigen::VectorXd EG_minimumPericenterRadii_;
    mutable Eigen::VectorXd EG_variableVector_;
    mutable double EG_departureTime_;
    mutable double EG_TOF_;
    mutable double EG_arrivalTime_;
    mutable Eigen::Vector6d EG_departureState_;
    mutable Eigen::Vector6d EG_arrivalState_;
    mutable BaseFunctionVector EG_radialVelocityFunctionComponents_;
    mutable BaseFunctionVector EG_normalVelocityFunctionComponents_;
    mutable BaseFunctionVector EG_axialVelocityFunctionComponents_;
    mutable Eigen::VectorXd EG_freeCoefficientsRadialVelocityFunction_;
    mutable Eigen::VectorXd EG_freeCoefficientsNormalVelocityFunction_;
    mutable Eigen::VectorXd EG_freeCoefficientsAxialVelocityFunction_;
    mutable int numEGvars_;
    mutable int numGMvars_;

    //GM parameters
    std::vector< B_LegTypeTransfer > GM_LegTypeVector_;
    std::vector< int > GM_flybysequence_;
    int GM_numLegs_;
    int GM_numDSM_;
    std::vector< ephemerides::EphemerisPointer > GM_ephemerisVector_;
    Eigen::VectorXd GM_gravitationalParameterVector_;
    Eigen::Vector2d GM_semiMajorAxes_;
    Eigen::Vector2d GM_eccentricities_;
    Eigen::VectorXd GM_minimumPericenterRadii_;
    mutable Eigen::VectorXd GM_variableVector_;
    mutable double GM_departureTime_;
    mutable double GM_TOF_;
    mutable double GM_arrivalTime_;
    mutable Eigen::Vector6d GM_departureState_;
    mutable Eigen::Vector6d GM_arrivalState_;
    mutable BaseFunctionVector GM_radialVelocityFunctionComponents_;
    mutable BaseFunctionVector GM_normalVelocityFunctionComponents_;
    mutable BaseFunctionVector GM_axialVelocityFunctionComponents_;
    mutable Eigen::VectorXd GM_freeCoefficientsRadialVelocityFunction_;
    mutable Eigen::VectorXd GM_freeCoefficientsNormalVelocityFunction_;
    mutable Eigen::VectorXd GM_freeCoefficientsAxialVelocityFunction_;



};
#endif // ANALYSISB_H
