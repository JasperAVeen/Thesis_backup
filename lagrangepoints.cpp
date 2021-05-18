#include "lagrangepoints.h"

#include <cmath>
#include <Eigen/Core>
#include "Tudat/Basics/basicTypedefs.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h>

#include <boost/filesystem.hpp>
#include <iostream>
#include <stdexcept>
#include <boost/bind.hpp>
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"

namespace tudat
{

namespace circular_restricted_three_body_problem
{

using namespace root_finders;
using namespace basic_mathematics;

//! Compute location of Lagrange libration point.
void LagrangePoint::computeLocationOfLibrationPoint(
        int LibNum )
{
    LagrangeLibrationPoints LagrangeLibrationPoint{};
    switch (LibNum){
    case(1): LagrangeLibrationPoint = l1; break;
    case(2): LagrangeLibrationPoint = l2; break;
    case(3): LagrangeLibrationPoint = l3; break;
    case(4): LagrangeLibrationPoint = l4; break;
    case(5): LagrangeLibrationPoint = l5; break;
    }
    using std::pow;
    using std::sqrt;

    // Set functions for Newton-Raphson based on collinear libration point passed as input
    // parameter, or computed locations directly of equilateral libration points.
    switch( LagrangeLibrationPoint )
    {
    case l1:
    {
        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = std::make_shared< UnivariateProxy >(
                    std::bind( &LagrangePoint::computeL1LocationFunction, this, std::placeholders::_1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, std::bind( &LagrangePoint::
                computeL1FirstDerivativeLocationFunction, this, std::placeholders::_1 ) );

        // Set position vector of L1 in Cartesian elements based on result of Newton-Raphson
        // root-finding algorithm.
        positionOfLibrationPoint_ << rootFinder->execute( rootFunction, 1.0 ), 0.0, 0.0;
    }
        break;

    case l2:
    {
        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = std::make_shared< UnivariateProxy >(
                    std::bind( &LagrangePoint::computeL2LocationFunction, this, std::placeholders::_1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, std::bind( &LagrangePoint::
                computeL2FirstDerivativeLocationFunction, this, std::placeholders::_1 ) );

        // Set position vector of L1 in Cartesian elements based on result of Newton-Raphson
        // root-finding algorithm.
        positionOfLibrationPoint_ << rootFinder->execute( rootFunction, 1.0 ), 0.0, 0.0;
    }
        break;

    case l3:
    {
        // Create an object containing the function of which we whish to obtain the root from.
        UnivariateProxyPointer rootFunction = std::make_shared< UnivariateProxy >(
                    std::bind( &LagrangePoint::computeL3LocationFunction, this, std::placeholders::_1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, std::bind( &LagrangePoint::
                computeL3FirstDerivativeLocationFunction, this, std::placeholders::_1 ) );

        // Set position vector of L1 in Cartesian elements based on result of Newton-Raphson
        // root-finding algorithm.
        positionOfLibrationPoint_ << rootFinder->execute( rootFunction, -1.0 ), 0.0, 0.0;
    }
        break;

    case l4:

        // Set position vector of L4 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = 0.5 - massParameter;
        positionOfLibrationPoint_.y( ) = 0.5 * sqrt( 3.0 );
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    case l5:

        // Set position vector of L5 in Cartesian elements.
        positionOfLibrationPoint_.x( ) = 0.5 - massParameter;
        positionOfLibrationPoint_.y( ) = -0.5 * sqrt( 3.0 );
        positionOfLibrationPoint_.z( ) = 0.0;

        break;

    default:{
        throw std::runtime_error(
                            "The Lagrange libration point requested does not exist." );}
    }
}

Eigen::Vector6d LagrangePoint::getLagrangeStateVec(LagrangePoint LibPoint, int LibNum, Eigen::Vector6d P1P2state){


    //creating state in corotating normalized frame
    LibPoint.computeLocationOfLibrationPoint( LibNum );
    Eigen::Vector3d relPosVec = LibPoint.getLocationOfLagrangeLibrationPoint( );

    Eigen::Vector3d relVelVec; relVelVec << 0, 0, 0 ;
    Eigen::Vector6d relStateVec; relStateVec <<  relPosVec, relVelVec;


    //Frame transformation into Cartesian frame
    Eigen::Vector6d normalizedCartesian = circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis( primGravPar, secGravPar,
                    relStateVec, P1P2state);


    return  normalizedCartesian;

}





} // namespace circular_restricted_three_body_problem

} // namespace tudat


/*
 double primsecDis{secPosVec.norm()};
   std::cout << primsecDis;
   double theta{atan(secPosVec(1)/secPosVec(0))};
   double phi{atan(secPosVec(2)/secPosVec(0))};
   Eigen::Vector3d LagPosVec(3);
   if(secPosVec(0)>=0 && secPosVec(1)>=0){
       LagPosVec << (cos(theta)*relPosVec(0) - sin(theta) * relPosVec(1) - sin(phi)*relPosVec(2))*primsecDis, (sin(theta)*relPosVec(0) + cos(theta) * relPosVec(1) + cos(phi)*relPosVec(2))*primsecDis, tan(phi)*(relPosVec(0));
   }
   if(secPosVec(0)<=0 && secPosVec(1)>=0){
       LagPosVec << (cos(theta)*relPosVec(0) - sin(theta) * relPosVec(1) - sin(phi)*relPosVec(2))*-1*primsecDis, (sin(theta)*relPosVec(0) + cos(theta) * relPosVec(1) + cos(phi)*relPosVec(2))*-1*primsecDis, tan(phi)*(relPosVec(0));
   }
   if(secPosVec(0)<=0 && secPosVec(1)<=0){
       LagPosVec << (cos(theta)*relPosVec(0) - sin(theta) * relPosVec(1) - sin(phi)*relPosVec(2))*-1*primsecDis, (sin(theta)*relPosVec(0) + cos(theta) * relPosVec(1) + cos(phi)*relPosVec(2))*-1*primsecDis, tan(phi)*(relPosVec(0));
   }
   if(secPosVec(0)>=0 && secPosVec(1)<=0){
       LagPosVec << (cos(theta)*relPosVec(0) - sin(theta) * relPosVec(1) - sin(phi)*relPosVec(2))*primsecDis, (sin(theta)*relPosVec(0) + cos(theta) * relPosVec(1) + cos(phi)*relPosVec(2))*primsecDis, tan(phi)*(relPosVec(0));
   }
 */
