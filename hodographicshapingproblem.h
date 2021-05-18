#ifndef HODOGRAPHICSHAPINGPROBLEM_H
#define HODOGRAPHICSHAPINGPROBLEM_H

#include <vector>
#include <utility>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "lagrangepoints.h"
#include "tt_functions.h"

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
#include <Eigen/Core>

typedef Eigen::Matrix< double, 6, 1 > StateType;

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace pagmo;


//! Test function for a new low-thrust trajectory class in Tudat
struct HodographicShapingProblem
{
    typedef std::vector< std::shared_ptr< tudat::shape_based_methods::BaseFunctionHodographicShaping > > BaseFunctionVector;

    HodographicShapingProblem( ){ }

    HodographicShapingProblem (
            const std::function< Eigen::Vector6d( const double ) >& initialStateFunction,
            const std::function< Eigen::Vector6d( const double ) >& finalStateFunction,
            const double centralBodyGravitationalParameter,
            const int numberOfRevolutions,
            const std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction,
            const std::vector< std::vector< double > >& freeCoefficientsBounds,
            std::string flybySequence,
            bool useDepOrbit,
            bool useArrOrbit,
            const bool minimizeMaximumThrust = false,
            const double initialMass = TUDAT_NAN
            ):
        initialStateFunction_( initialStateFunction ),
        finalStateFunction_( finalStateFunction ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        numberOfRevolutions_( numberOfRevolutions ),
        basisFunctionsFunction_( basisFunctionsFunction ),
        problemBounds_( freeCoefficientsBounds ),
        flybySequence_(flybySequence),
        useDepOrbit_(useDepOrbit),
        useArrOrbit_(useArrOrbit),
        minimizeMaximumThrust_( minimizeMaximumThrust ),
        initialMass_( initialMass )
    {}

    // Calculates the fitness
    std::vector< double > fitness( const std::vector< double > &x ) const;

    std::pair< std::vector< double >, std::vector< double > > get_bounds() const
    {
        return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
    }

    template <typename Archive>
    void serialize(Archive &ar)
    {
        ar(problemBounds_);
    }

    vector_double::size_type get_nobj() const
    {
        return minimizeMaximumThrust_ ? 2u : 1u;
    }

protected:

private:

    const std::function< Eigen::Vector6d( const double ) > initialStateFunction_;

    const std::function< Eigen::Vector6d( const double ) > finalStateFunction_;

    double centralBodyGravitationalParameter_;

    int numberOfRevolutions_;

    std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction_;

    const std::vector< std::vector< double > > problemBounds_;

    std::string flybySequence_;

    bool useDepOrbit_;

    bool useArrOrbit_;

    bool minimizeMaximumThrust_;

    double initialMass_;

};

#endif // HODOGRAPHICSHAPINGPROBLEM_H
