#include "hodographicshapingproblem.h"

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::low_thrust_trajectories;
using namespace pagmo;


// Calculates the fitness
std::vector< double > HodographicShapingProblem::fitness( const std::vector< double > &x ) const
{
    double departureTime = x.at( 0 );
    double timeOfFlight = x.at( 1 );
    double arrivalTime = departureTime + timeOfFlight;

    std::vector< BaseFunctionVector > basisFunctions = basisFunctionsFunction_( timeOfFlight );

    BaseFunctionVector radialVelocityFunctionComponents = basisFunctions.at( 0 );
    BaseFunctionVector normalVelocityFunctionComponents = basisFunctions.at( 1 );
    BaseFunctionVector axialVelocityFunctionComponents = basisFunctions.at( 2 );

    int numberFreeCoefficientsRadialFunction = radialVelocityFunctionComponents.size( ) - 3;
    int numberFreeCoefficientsNormalFunction = normalVelocityFunctionComponents.size( ) - 3;
    int numberFreeCoefficientsAxialFunction = axialVelocityFunctionComponents.size( ) - 3;

    /*
    if( numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 2 !=
            static_cast< int >( x.size( ) ) )
    {
        throw std::runtime_error( "Error, size of design variables vector unconsistent with number of base function components"
                                  "when making a hodographic shaping optimisation problem." );
    }
    */

    Eigen::VectorXd freeCoefficientsRadialVelocityFunction( numberFreeCoefficientsRadialFunction );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction( numberFreeCoefficientsNormalFunction );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction( numberFreeCoefficientsAxialFunction );

    for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
    {
        freeCoefficientsRadialVelocityFunction[ i ] = x[ i + 2 ];
    }
    for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
    {
        freeCoefficientsNormalVelocityFunction[ i ] = x[ i + 2 + numberFreeCoefficientsRadialFunction ];
    }
    for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
    {
        freeCoefficientsAxialVelocityFunction[ i ] = x[ i + 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
    }

    Eigen::Vector3d Vinf_vec_dep = Eigen::Vector3d::Zero();
    Eigen::Vector3d Vinf_vec_arr = Eigen::Vector3d::Zero();
    double departureDeltaV = 0;
    double arrivalDeltaV = 0;
    Eigen::Vector6d departureState = getDepartureState(initialStateFunction_, departureTime, flybySequence_, Vinf_vec_dep);
    Eigen::Vector6d arrivalState = getArrivalState(finalStateFunction_, arrivalTime, flybySequence_, Vinf_vec_arr);

    if(useDepOrbit_ && !useArrOrbit_){
    Vinf_vec_dep <<  x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction ),
                                           x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 1 ),
                                           x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 2 );

    departureState = getDepartureState(initialStateFunction_, departureTime, flybySequence_, Vinf_vec_dep);
    departureDeltaV = getDepDeltaV(Vinf_vec_dep, flybySequence_);
    }

    if(useArrOrbit_ && !useDepOrbit_){
    Vinf_vec_arr <<  x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction ),
                                           x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 1 ),
                                           x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 2 );

     arrivalState = getArrivalState(finalStateFunction_, arrivalTime, flybySequence_, Vinf_vec_arr);
     arrivalDeltaV = getArrDeltaV(Vinf_vec_arr, flybySequence_);

    }
    if(useDepOrbit_ && useArrOrbit_){
        Vinf_vec_dep <<  x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction ),
                                               x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 1 ),
                                               x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 2 );

        departureState = getDepartureState(initialStateFunction_, departureTime, flybySequence_, Vinf_vec_dep);
        departureDeltaV = getDepDeltaV(Vinf_vec_dep, flybySequence_);


        Vinf_vec_arr <<  x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction +3 ),
                                               x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 4 ),
                                               x.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 5 );

         arrivalState = getArrivalState(finalStateFunction_, arrivalTime, flybySequence_, Vinf_vec_arr);
         arrivalDeltaV = getArrDeltaV(Vinf_vec_arr, flybySequence_);


    }





        tudat::shape_based_methods::HodographicShaping hodographicShaping = tudat::shape_based_methods::HodographicShaping(
                    departureState, arrivalState,
                    timeOfFlight, centralBodyGravitationalParameter_,
                    numberOfRevolutions_, radialVelocityFunctionComponents,
                    normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                    freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
                    freeCoefficientsAxialVelocityFunction, initialMass_ );

        std::vector< double > fitnessVector;

        fitnessVector.push_back( hodographicShaping.computeDeltaV( ) + departureDeltaV + arrivalDeltaV );

        if( minimizeMaximumThrust_ )
        {
            Eigen::VectorXd epochsVector = Eigen::VectorXd::LinSpaced( 100, 0, timeOfFlight );

            double maximumAcceleration = 0.0;
            for( int i = 0; i < epochsVector.rows( ); i++ )
            {
                double currentAcceleration =
                        hodographicShaping.computeCurrentThrustAccelerationMagnitude( epochsVector( i ) );
                if( maximumAcceleration < currentAcceleration )
                {
                    maximumAcceleration = currentAcceleration;
                }
            }
            fitnessVector.push_back( maximumAcceleration );

        }

        return fitnessVector;

    }






