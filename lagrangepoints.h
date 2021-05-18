#ifndef LAGRANGEPOINTS_H
#define LAGRANGEPOINTS_H

#include <cmath>
#include <memory>
#include <Eigen/Core>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"


namespace tudat
{

namespace circular_restricted_three_body_problem
{

//! Compute mass parameter.
/*!
 * Computes dimensionless mass parameter for the CRTBP.
 * \param primaryGravitationalParameter Gravitational parameter of primary body.
 * \param secondaryGravitationalParameter Gravitational parameter of secondary body.
 * \return Dimensionless mass parameter.
 */
inline double computeMassParameter( const double primaryGravitationalParameter,
                                    const double secondaryGravitationalParameter )
{
    return secondaryGravitationalParameter
            / ( primaryGravitationalParameter + secondaryGravitationalParameter );
}

//! Libration point class.
/*!
 * This class includes functions to compute the location of Lagrange libration points in the CRTBP.
 */
class LagrangePoint
{

public:

    //! Lagrange libration points.
    enum LagrangeLibrationPoints { l1, l2, l3, l4, l5 };

    //! Default constructor.
    LagrangePoint( const double aPrimaryGravitationalParameter,
                    const double aSecondaryGravitationalParameter,
                    const root_finders::RootFinderPointer aRootFinder )
        : massParameter( computeMassParameter( aPrimaryGravitationalParameter,
                                               aSecondaryGravitationalParameter ) ), primGravPar(aPrimaryGravitationalParameter), secGravPar(aSecondaryGravitationalParameter),
          rootFinder( aRootFinder )
    { }

    //! Default constructor.
    /*!
     * Default constructor.
     * \param aMassParameter Dimensionless mass parameter of the smaller of the massive bodies in the CRTBP.
     * \param aRootFinder Shared pointer to the rootfinder which is used for finding the L1, L2 and
     * L3 libration points. The rootfinder contains termination conditions inside.
     */

    LagrangePoint( const double aMassParameter,
                    const root_finders::RootFinderPointer aRootFinder )
        : massParameter( aMassParameter ),
          rootFinder( aRootFinder )
    { }


    //! Get dimensionless mass parameter.
    /*!
     * Returns the dimensionless mass parameter based on the gravitational
     * parameters of the primary and secondary bodies.
     * \return Dimensionless mass parameter.
     */
    double getMassParameter( ) { return massParameter; }

    //! Get location of Lagrange libration point.
    /*!
     * Returns the position vector in Cartesian elements of a Lagrange libration point. The
     * libration point is selected when the computeLocationOfLibrationPoint( ) function is
     * called.
     * \return Cartesian position elements of Lagrange libration point.
     */
    Eigen::Vector3d getLocationOfLagrangeLibrationPoint( )
    {
        return positionOfLibrationPoint_;
    }
    const double massParameter;


    //! Compute location of Lagrange libration point.
    /*!
     * Computes the location as a position vector in Cartesian elements of a given Lagrange
     * libration point.
     * \param lagrangeLibrationPoint Lagrange libration point.
     */
    void computeLocationOfLibrationPoint( int LibNum );

    Eigen::Vector6d getLagrangeStateVec(LagrangePoint, int, Eigen::Vector6d);

protected:

private:

    //! Dimensionless mass parameter.
    /*!
    * Dimensionless mass parameter of the smaller of the massive bodies in the CRTBP.
    */

     double primGravPar;

     double secGravPar;

    //! Shared pointer to the rootfinder.
    /*!
     * Shared pointer to the rootfinder. The rootfinder contains termination conditions inside.
     */
    const root_finders::RootFinderPointer rootFinder;

    //! Cartesian position elements of a Lagrange libration point.
    /*!
     * Cartesian position elements of a Lagrange libration point.
     */
    Eigen::Vector3d positionOfLibrationPoint_;

    //! Compute equation of motion for location of L1 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L1 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L1 libration point.
     * \return Value of location function.
     */
    double computeL1LocationFunction( const double xLocationEstimate )
    {
        return xLocationEstimate -
                ( 1.0 - massParameter ) / std::pow( massParameter + xLocationEstimate, 2.0 )
                + massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L1 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L1 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L1 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL1FirstDerivativeLocationFunction( const double xLocationEstimate )
    {
        return 1.0 + 2.0 * ( 1.0 - massParameter )
                / std::pow( massParameter + xLocationEstimate, 3.0 )
                + 2.0 * massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 3.0 );
    }

    //! Compute equation of motion for location of L2 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L2 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L2 libration point.
     * \return Value of location function.
     */
    double computeL2LocationFunction( const double xLocationEstimate )
    {
        return xLocationEstimate
                - ( 1.0 - massParameter ) / std::pow( massParameter + xLocationEstimate, 2.0 )
                - massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L2 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L2 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L2 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL2FirstDerivativeLocationFunction( const double xLocationEstimate )
    {
        return 1.0 + 2.0 * ( 1.0 - massParameter )
                / std::pow( massParameter + xLocationEstimate, 3.0 )
                - 2.0 * massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 3.0 );
    }

    //! Compute equation of motion for location of L3 libration point.
    /*!
     * Computes equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L3 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L3 libration point.
     * \return Value of location function.
     */
    double computeL3LocationFunction( const double xLocationEstimate )
    {
        return xLocationEstimate
                + ( 1.0 - massParameter ) / std::pow( massParameter + xLocationEstimate, 2.0 )
                - massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 2.0 );
    }

    //! Compute derivative of equation of motion for location of L3 libration point.
    /*!
     * Computes derivative of equation of motion, whose root is the x-location in dimensionless
     * coordinates of colinear libration point L3 in the CRTBP (van der Ham, 2012).
     * \param xLocationEstimate Estimate of x-location of L3 libration point.
     * \return Value of first derivative of location function.
     */
    double computeL3FirstDerivativeLocationFunction( const double xLocationEstimate )
    {
        return 1.0 - 2.0 * ( 1.0 - massParameter )
                / std::pow( massParameter + xLocationEstimate, 3.0 )
                -2.0 * massParameter / std::pow( 1.0 - massParameter - xLocationEstimate, 3.0 );
    }
};

// Typedef for shared-pointer to LibrationPoint object.
typedef std::shared_ptr< LagrangePoint > LibrationPointPointer;

} // namespace circular_restricted_three_body_problem
} // namespace tudat


#endif // LAGRANGEPOINTS_H
