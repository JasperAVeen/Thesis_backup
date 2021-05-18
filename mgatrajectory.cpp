#include "mgatrajectory.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"


using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
using namespace tudat;
using namespace pagmo;


MGAtrajectory::MGAtrajectory( std::vector< std::vector< double > > &bounds,
                              std::vector< int > flybySequence, std::vector< LegTypeTransfer > LegTypeVec, std::string flyby_seq, double TOF_constraint, Eigen::VectorXd semimajoraxes,
                      Eigen::VectorXd eccentricities, double end_date_constraint, const bool useTripTime):
    problemBounds_( bounds ), flybyseq_{flyby_seq}, TOFconst_( TOF_constraint ),  useTripTime_( useTripTime ), legTypeVector_(LegTypeVec), semiMajorAxes_{semimajoraxes}, eccentricities_{eccentricities}, end_date_constraint_(end_date_constraint)
{

    // Create the ephemeris, gravitational parameter, and minimum pericentre vector.
    numberOfLegs_ = flybySequence.size( );
    ephemerisVector_.resize( numberOfLegs_ );
    gravitationalParameterVector_.resize( numberOfLegs_ );
    minimumPericenterRadii_.resize( numberOfLegs_ );
    ephemerisVector_ = createEphemerisVec (flybySequence);
    gravitationalParameterVector_=createGravParamVec(flybySequence);

    int num_DSMs = std::count(flyby_seq.begin(), flyby_seq.end(), 'd');


    if(num_DSMs==0) {minimumPericenterRadii_=createMinPerRad(flybySequence);}
    else {
        for (int i{0}; i < numberOfLegs_; i++){
            minimumPericenterRadii_ [i] = TUDAT_NAN;
        }

    }

}

//! Descriptive name of the problem
std::string MGAtrajectory::get_name() const {
    return "MGA transfer trajectory";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > MGAtrajectory::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}

//! Implementation of the fitness function (return delta-v)
std::vector<double> MGAtrajectory::fitness( const std::vector<double> &xv ) const{
    int num_DSMs = std::count(flybyseq_.begin(), flybyseq_.end(), 'd');


    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create variable vector.
    Eigen::VectorXd variableVector(numberOfLegs_+1+(4*num_DSMs));

    double TOF = 0;
    for(int i = 0; i < numberOfLegs_ ; i++){
        variableVector[ i ] = xv[ i ];
        if( i > 0 ){
            TOF += xv[i];
        }
    }

    double end_date = variableVector[ 0 ] + TOF;


    variableVector *= physical_constants::JULIAN_DAY;
    variableVector[ numberOfLegs_ ] = 1;//dummy

    //Write the variable vector for DSM variables below

    if(num_DSMs>0){
        int count {numberOfLegs_+1};
        for(size_t i{0}; i<num_DSMs; i++){
            variableVector[ count ] = xv[ count-1 ];
            variableVector[ count + 1 ] = xv[ count ];
            variableVector[ count + 2 ] = xv[ count + 1 ];
            variableVector[ count + 3 ] = xv[ count + 2 ];
            count+=4;
        }
    }

    //std::cout <<"\nVariable vector: \n" << variableVector << std::endl;


    //Write the variable vector for DSM variables above


    // Create the trajectory problem.
    makeTraj mgaTraj( numberOfLegs_, legTypeVector_, flybyseq_, ephemerisVector_,
                          gravitationalParameterVector_, variableVector, sunGravitationalParameter,
                          minimumPericenterRadii_, semiMajorAxes_, eccentricities_ );

    //std::cout << "\n ---------- Trajectory created";


    // Start the deltaV vector.
    double resultingDeltaV;
    mgaTraj.calculateTrajectory( resultingDeltaV );

    //std::cout << "\n ---------- Trajectory DeltaV computed\n\n";



    if (end_date>end_date_constraint_)
    {
        resultingDeltaV = 1.0E10;
    }
    if (TOF>TOFconst_)
    {
        resultingDeltaV = 1.0E10;
    }

    if ( useTripTime_ ){
        return { resultingDeltaV, TOF };
    }
    else {
        return { resultingDeltaV };
    }



}
