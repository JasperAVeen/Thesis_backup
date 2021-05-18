#include "analysisb.h"
#include "tt_functions.h"

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
using namespace tudat::low_thrust_trajectories;

using namespace tudat;
using namespace pagmo;



AnalysisB::AnalysisB( std::vector<std::string> EM_trajectory, std::vector< std::vector< double > > Bounds,  std::unordered_map<int, std::map<double, Eigen::VectorXd>> orbitLibrary, bool multiObj )
{
    comments_ = false;
    multiObj_ = multiObj;

    orbitLibrary_ = orbitLibrary;

    EG_trajectory_ = EM_trajectory.at(0);
    GM_trajectory_ = EM_trajectory.at(1);
    problemBounds_ = Bounds;
    EG_isCont_ = false;
    GM_isCont_ = false;
    if (EG_trajectory_.find('c') != std::string::npos)
        EG_isCont_ = true;
    if (GM_trajectory_.find('c') != std::string::npos)
        GM_isCont_ = true;

    EG_numDSM_ = std::count(EG_trajectory_.begin(), EG_trajectory_.end(), 'd');
    GM_numDSM_ = std::count(GM_trajectory_.begin(), GM_trajectory_.end(), 'd');


    //If segments are IT, define all needed parameters for trajectory computation
    if(!EG_isCont_) {
        EG_LegTypeVector_ = createLegType_B(EG_trajectory_);
        EG_flybysequence_ = createFlyBySeq(EG_trajectory_);
        EG_numLegs_ = EG_flybysequence_.size();
        EG_ephemerisVector_.resize( EG_numLegs_ );
        EG_gravitationalParameterVector_.resize( EG_numLegs_ );
        EG_minimumPericenterRadii_.resize( EG_numLegs_ );
        EG_ephemerisVector_ = createEphemerisVec (EG_flybysequence_);
        EG_gravitationalParameterVector_=createGravParamVec(EG_flybysequence_);
        if(EG_numDSM_==0) {EG_minimumPericenterRadii_=createMinPerRad(EG_flybysequence_);}
        else {
            for (int i{0}; i < EG_numLegs_; i++){
                EG_minimumPericenterRadii_ [i] = TUDAT_NAN;
            }
        }
        EG_semiMajorAxes_ = getSMA(EG_trajectory_);
        EG_eccentricities_ = getEccentricities(EG_trajectory_);
    }

    if(!GM_isCont_) {
        GM_LegTypeVector_ = createLegType_B(GM_trajectory_);
        GM_flybysequence_ = createFlyBySeq(GM_trajectory_);
        GM_numLegs_ = GM_flybysequence_.size();
        GM_ephemerisVector_.resize( GM_numLegs_ );
        GM_gravitationalParameterVector_.resize( GM_numLegs_ );
        GM_minimumPericenterRadii_.resize( GM_numLegs_ );
        GM_ephemerisVector_ = createEphemerisVec (GM_flybysequence_);
        GM_gravitationalParameterVector_=createGravParamVec(GM_flybysequence_);
        if(GM_numDSM_==0) {GM_minimumPericenterRadii_=createMinPerRad(GM_flybysequence_);}
        else {
            for (int i{0}; i < GM_numLegs_; i++){
                GM_minimumPericenterRadii_ [i] = TUDAT_NAN;
            }
        }
        GM_semiMajorAxes_ = getSMA(GM_trajectory_);
        GM_eccentricities_ = getEccentricities(GM_trajectory_);
    }

    if(EG_isCont_) {
        EG_flybysequence_ = createFlyBySeq(EG_trajectory_);
        EG_ephemerisVector_.resize( 2 );
        EG_ephemerisVector_ = createEphemerisVec ( EG_flybysequence_);
    }
    if(GM_isCont_) {
        GM_flybysequence_ = createFlyBySeq(EG_trajectory_);
        GM_ephemerisVector_.resize( 2 );
        GM_ephemerisVector_ = createEphemerisVec ( GM_flybysequence_);
    }

    if(comments_) std::cout << "\n+ Analysis B initiated";

}

//! Descriptive name of the problem
std::string AnalysisB::get_name() const {
    return "Analysis B";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > AnalysisB::get_bounds() const {    

    if(comments_) std::cout << "\n+ Bounds used";


    return { problemBounds_[0], problemBounds_[1] };
}


//! Implementation of the fitness function (return delta-v)
std::vector<double> AnalysisB::fitness( const std::vector<double> &xv ) const {
    if(comments_) {
    std::cout << "\n======== xv: ========\n";

    for(size_t xv_iter {0}; xv_iter < xv.size(); xv_iter++) {
        std::cout << xv.at(xv_iter) << " ";
    }
}
    if(comments_) std::cout << "\n+ Inside fitness";

    int xv_size = xv.size();
    int lastElement = xv_size - 1;
    int secondToLastElement = xv_size - 2;

    //orbit modelling parameters
    orbitID_ = round(xv.at(secondToLastElement));
    specificOrbit_.clear();
    specificOrbit_ = orbitLibrary_[orbitID_];
    orbitFraction_= xv.at(lastElement);

    //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////    Segment EG     //////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Fill variable vector for each transfer type
    if(comments_) std::cout << "\n+ Inside Segment EG - general";

    if(!EG_isCont_) {
        if(comments_) std::cout << "\n+ Inside Segment EG - IT";

        EG_variableVector_.resize(EG_numLegs_+1+(4*EG_numDSM_));

        EG_TOF_ = 0;
        for(int i = 0; i < EG_numLegs_ ; i++){
            EG_variableVector_[ i ] = xv[ i ];
            if( i > 0 ){
                EG_TOF_ += xv[i];
            }
        }

        EG_TOF_ *= physical_constants::JULIAN_DAY; //seconds
        EG_variableVector_ *= physical_constants::JULIAN_DAY;
        EG_variableVector_[ EG_numLegs_ ] = 1;//dummy


        EG_departureTime_ = EG_variableVector_[ 0 ]; //seconds
        EG_arrivalTime_ = EG_TOF_ + EG_departureTime_; //seconds

        if(EG_numDSM_>0){
            int count {EG_numLegs_+1};
            for(size_t i{0}; i<EG_numDSM_; i++){
                EG_variableVector_[ count ] = xv[ count-1 ];
                EG_variableVector_[ count + 1 ] = xv[ count ];
                EG_variableVector_[ count + 2 ] = xv[ count + 1 ];
                EG_variableVector_[ count + 3 ] = xv[ count + 2 ];
                count+=4;
            }
        }

        numEGvars_ = EG_variableVector_.size() - 1;

        if(comments_) {
        std::cout << "\n======== EG variable vector: ========\n";
        std::cout << EG_variableVector_ << std::endl;
        std::cout << "\nSize: " << numEGvars_;
        }

        if(comments_) std::cout << "\n+ Completed Segment EG - IT";

    }

    if(EG_isCont_) {        
        EG_departureTime_ = xv.at( 0 ) * physical_constants::JULIAN_DAY; // seconds
        EG_TOF_ = xv.at( 1 ) * physical_constants::JULIAN_DAY; // seconds
        EG_arrivalTime_ = EG_departureTime_ + EG_TOF_; //seconds

        std::vector< BaseFunctionVector > basisFunctions = basisFunctionsFunction_( EG_TOF_ );

        EG_radialVelocityFunctionComponents_ = basisFunctions.at( 0 );
        EG_normalVelocityFunctionComponents_ = basisFunctions.at( 1 );
        EG_axialVelocityFunctionComponents_ = basisFunctions.at( 2 );

        int numberFreeCoefficientsRadialFunction = EG_radialVelocityFunctionComponents_.size( ) - 3;
        int numberFreeCoefficientsNormalFunction = EG_normalVelocityFunctionComponents_.size( ) - 3;
        int numberFreeCoefficientsAxialFunction = EG_axialVelocityFunctionComponents_.size( ) - 3;

        EG_freeCoefficientsRadialVelocityFunction_.resize( numberFreeCoefficientsRadialFunction );
        EG_freeCoefficientsNormalVelocityFunction_.resize( numberFreeCoefficientsNormalFunction );
        EG_freeCoefficientsAxialVelocityFunction_.resize( numberFreeCoefficientsAxialFunction );

        for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
        {
            EG_freeCoefficientsRadialVelocityFunction_[ i ] = xv[ i + 2 ];
        }
        for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
        {
            EG_freeCoefficientsNormalVelocityFunction_[ i ] = xv[ i + 2 + numberFreeCoefficientsRadialFunction ];
        }
        for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
        {
            EG_freeCoefficientsAxialVelocityFunction_[ i ] = xv[ i + 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
        }

        Vinf_vec_dep_ = Eigen::Vector3d::Zero();
        Vinf_vec_dep_ <<  xv.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction ),
                                              xv.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 1 ),
                                               xv.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 2 );


        EG_departureState_ = B_getDepartureState(EG_departureTime_, EG_trajectory_, EG_ephemerisVector_, specificOrbit_,  orbitFraction_, stayTime_, Vinf_vec_dep_);
        EG_arrivalState_ =  B_getArrivalState(EG_arrivalTime_, EG_trajectory_, EG_ephemerisVector_,  specificOrbit_,  orbitFraction_);

        numEGvars_ = 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + V_inf_vars_;

    }

    //Compute trajectories for each IT
    if(!EG_isCont_) {
        if(comments_) std::cout << "\n+ Inside Segment EG - computation";

        B_makeTraj EG_IT( EG_numLegs_, EG_LegTypeVector_, EG_trajectory_, EG_ephemerisVector_,
                              EG_gravitationalParameterVector_, EG_variableVector_, sunGravitationalParameter_,
                              EG_minimumPericenterRadii_, EG_semiMajorAxes_, EG_eccentricities_, specificOrbit_,  orbitFraction_ );

        EG_IT.calculateTrajectory( EG_fitnessvalue_ );
        if(comments_) std::cout << "\n+ Completed Segment EG - computation";

    }

    if(EG_isCont_) {
        tudat::shape_based_methods::HodographicShaping EG_CT = tudat::shape_based_methods::HodographicShaping(
                    EG_departureState_, EG_arrivalState_,
                    EG_TOF_, sunGravitationalParameter_,
                    numberOfRevolutions_, EG_radialVelocityFunctionComponents_,
                    EG_normalVelocityFunctionComponents_, EG_axialVelocityFunctionComponents_,
                    EG_freeCoefficientsRadialVelocityFunction_, EG_freeCoefficientsNormalVelocityFunction_,
                    EG_freeCoefficientsAxialVelocityFunction_, initialMass_ );        

       EG_trajectoryDV_ = EG_CT.computeDeltaV( );
       EG_escapeDV_ = getDepDeltaV(Vinf_vec_dep_, EG_trajectory_);
       EG_fitnessvalue_ = EG_trajectoryDV_ + EG_escapeDV_;
    }



    //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////    Segment GM     //////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    if(comments_) std::cout << "\n+ Inside Segment GM - general";

    stayTime_ = xv[numEGvars_];
    stayTimeSec_ = stayTime_ * physical_constants::JULIAN_DAY;


    if(!GM_isCont_) {
        if(comments_) std::cout << "\n+ Inside Segment GM - IT";

        GM_variableVector_.resize(GM_numLegs_+1+(4*GM_numDSM_));

        GM_variableVector_[ 0 ] =  (EG_arrivalTime_ /  physical_constants::JULIAN_DAY) + stayTime_; //Departure time = arrival time EG PLUS stay time (xv[numEGvars_])

        GM_TOF_ = 0;
        for(int i = 1; i < GM_numLegs_ ; i++){
            GM_variableVector_[ i ] = xv[ numEGvars_ + i ];
                GM_TOF_ += xv[numEGvars_ + i];
        }

        GM_TOF_ *= physical_constants::JULIAN_DAY; //seconds
        GM_variableVector_ *= physical_constants::JULIAN_DAY;
        GM_arrivalTime_ = GM_variableVector_[ 0 ] + GM_TOF_ ; //seconds
        GM_variableVector_[ GM_numLegs_ ] = 1;//dummy

        if(GM_numDSM_>0){
            int XVcount {numEGvars_ + GM_numLegs_};
            int GMVarCount {GM_numLegs_ + 1};

            for(size_t i{0}; i<GM_numDSM_; i++){
                GM_variableVector_[ GMVarCount ] = xv[ XVcount ];
                GM_variableVector_[ GMVarCount + 1 ] = xv[ XVcount + 1 ];
                GM_variableVector_[ GMVarCount + 2 ] = xv[ XVcount + 2 ];
                GM_variableVector_[ GMVarCount + 3 ] = xv[ XVcount + 3 ];
                GMVarCount+=4;
                XVcount+=4;
            }
        }

        numGMvars_ = GM_variableVector_.size() - 1;

        if(comments_) {
        std::cout << "\n======== GM variable vector: ========\n";
        std::cout << GM_variableVector_ << std::endl;
        std::cout << "\nSize: " << numGMvars_;
        }

        if(comments_) std::cout << "\n+ Completed Segment GM - IT";

    }

    if(GM_isCont_) {
        GM_departureTime_ = (stayTime_ * physical_constants::JULIAN_DAY  ) + EG_arrivalTime_; //Departure time = arrival time EG PLUS stay time (xv[numEGvars_]) //seconds
        GM_TOF_ = xv.at( numEGvars_ + 1 ) * physical_constants::JULIAN_DAY; //seconds
        GM_arrivalTime_ = GM_departureTime_ + GM_TOF_; //seconds

        std::vector< BaseFunctionVector > basisFunctions = basisFunctionsFunction_( GM_TOF_ );

        GM_radialVelocityFunctionComponents_ = basisFunctions.at( 0 );
        GM_normalVelocityFunctionComponents_ = basisFunctions.at( 1 );
        GM_axialVelocityFunctionComponents_ = basisFunctions.at( 2 );

        int numberFreeCoefficientsRadialFunction = GM_radialVelocityFunctionComponents_.size( ) - 3;
        int numberFreeCoefficientsNormalFunction = GM_normalVelocityFunctionComponents_.size( ) - 3;
        int numberFreeCoefficientsAxialFunction = GM_axialVelocityFunctionComponents_.size( ) - 3;

        GM_freeCoefficientsRadialVelocityFunction_.resize( numberFreeCoefficientsRadialFunction );
        GM_freeCoefficientsNormalVelocityFunction_.resize( numberFreeCoefficientsNormalFunction );
        GM_freeCoefficientsAxialVelocityFunction_.resize( numberFreeCoefficientsAxialFunction );

        for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
        {
            GM_freeCoefficientsRadialVelocityFunction_[ i ] = xv[ numEGvars_ + i + 2 ];
        }
        for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
        {
            GM_freeCoefficientsNormalVelocityFunction_[ i ] = xv[ numEGvars_ +  i + 2 + numberFreeCoefficientsRadialFunction ];
        }
        for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
        {
            GM_freeCoefficientsAxialVelocityFunction_[ i ] = xv[ numEGvars_ + i + 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];

        }

        Vinf_vec_arr_ = Eigen::Vector3d::Zero();
        Vinf_vec_arr_ <<  xv.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction ),
                                               xv.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 1 ),
                                               xv.at( 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 2 );

        GM_departureState_ = B_getDepartureState(GM_departureTime_, GM_trajectory_, GM_ephemerisVector_, specificOrbit_,  orbitFraction_, stayTime_);
        GM_arrivalState_ =  B_getArrivalState(GM_arrivalTime_, GM_trajectory_, GM_ephemerisVector_, specificOrbit_,  orbitFraction_, Vinf_vec_arr_);

        numGMvars_ = 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction;
    }

    //Compute trajectories for each IT
    if(!GM_isCont_) {
        if(comments_) std::cout << "\n+ Inside Segment GM - computation";

        B_makeTraj GM_IT( GM_numLegs_, GM_LegTypeVector_, GM_trajectory_, GM_ephemerisVector_,
                              GM_gravitationalParameterVector_, GM_variableVector_, sunGravitationalParameter_,
                              GM_minimumPericenterRadii_, GM_semiMajorAxes_, GM_eccentricities_, specificOrbit_,  orbitFraction_, stayTime_ );
        GM_IT.calculateTrajectory( GM_fitnessvalue_ );
        if(comments_) std::cout << "\n+ Completed Segment GM - computation";

    }

    if(GM_isCont_) {

        tudat::shape_based_methods::HodographicShaping GM_CT = tudat::shape_based_methods::HodographicShaping(
                    GM_departureState_, GM_arrivalState_,
                    GM_TOF_, sunGravitationalParameter_,
                    numberOfRevolutions_, GM_radialVelocityFunctionComponents_,
                    GM_normalVelocityFunctionComponents_, GM_axialVelocityFunctionComponents_,
                    GM_freeCoefficientsRadialVelocityFunction_, GM_freeCoefficientsNormalVelocityFunction_,
                    GM_freeCoefficientsAxialVelocityFunction_, initialMass_ );

        GM_trajectoryDV_ = GM_CT.computeDeltaV( );
        GM_captureDV_ = getArrDeltaV(Vinf_vec_arr_, GM_trajectory_);
        GM_fitnessvalue_ = GM_trajectoryDV_ + GM_captureDV_;

    }

    if(comments_) std::cout << "\n+ Before fitness vector";

totalFitnessValue1 = GM_fitnessvalue_ + EG_fitnessvalue_; //DV in m/s
totalFitnessValue2 = (EG_TOF_ + GM_TOF_ + stayTimeSec_) / 3600; //TOF in hours

std::vector<double> fitness{};
fitness.push_back(totalFitnessValue1);

if(multiObj_){
    fitness.push_back(totalFitnessValue2);
}

if(comments_) std::cout << "\n+ After fitness vector";


return fitness;
}
