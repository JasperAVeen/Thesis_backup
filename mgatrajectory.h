#ifndef MGATRAJECTORY_H
#define MGATRAJECTORY_H

#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include <Tudat/InputOutput/basicInputOutput.h>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"

#include <random>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"
#include <pagmo/rng.hpp>
#include "tt_functions.h"
#include "maketraj.h"


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

struct MGAtrajectory
{
    MGAtrajectory( const bool useTripTime = false ): useTripTime_( useTripTime ){ }

    MGAtrajectory( std::vector< std::vector< double > > &bounds,
                           std::vector< int > flybySequence, std::vector< LegTypeTransfer > LegTypeVec, std::string flyby_seq, double TOF_constraint,  Eigen::VectorXd semimajoraxes,
                   Eigen::VectorXd eccentricities, double end_date_constraint=10838.5, const bool useTripTime = false );

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    std::string get_name( ) const;

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        if(useTripTime_ )
        {
            return 2u;
        }
        else
        {
            return 1u;
        }

    }

private:

    const std::vector< std::vector< double > > problemBounds_;
    const double TOFconst_{};
    const double end_date_constraint_{};

    bool useTripTime_;

    int numberOfLegs_;
    std::vector< LegTypeTransfer > legTypeVector_;
    std::vector< std::string > bodyNamesVector_;
    std::string flybyseq_;
    std::vector< ephemerides::EphemerisPointer > ephemerisVector_;
    Eigen::VectorXd gravitationalParameterVector_;
    Eigen::VectorXd semiMajorAxes_;
    Eigen::VectorXd eccentricities_;
    Eigen::VectorXd minimumPericenterRadii_;
};


#endif // MGATRAJECTORY_H
