#include "analysisb_functions.h"

Eigen::Vector6d B_getDepartureState ( double Time, std::string trajectory, std::vector< ephemerides::EphemerisPointer > ephemerisVector, std::map<double, Eigen::VectorXd> specificLPO, double orbitFraction_, double stayTime_, Eigen::Vector3d Vinf_vec){

    Eigen::Vector6d boundaryState = ephemerisVector[0]->getCartesianState(Time);

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char bound = trajectory.at(0);
    int switcharg{};
    if(bound=='E')switcharg=13;   if(bound=='m')switcharg=10;   if(bound=='M')switcharg=14;  //Earth & Moon & Mars
    if(bound=='1')switcharg=1;    if(bound=='2')switcharg=2;    if(bound=='3')switcharg=3;   //Lagrange points
    if(bound=='4')switcharg=4;    if(bound=='5')switcharg=5;    if(bound=='6')switcharg=6;   //Lagrange points
    if(bound=='7')switcharg=7;    if(bound=='8')switcharg=8;    if(bound=='9')switcharg=9;   //centralOrbits
    if(bound!='E'&& bound!='m' && bound!='M' && bound!='1'&&bound!='2'&&bound!='3'&&bound!='4'&&bound!='5'&&bound!='6'&&bound!='7'&&bound!='8'&&bound!='9')std::cout << "\nNo switcharg was generated!\n";

    switch (switcharg) {

        case(13):{ //Earth - parking orbit at 250 km
            boundaryState( 3 ) += Vinf_vec[0];
            boundaryState( 4 ) += Vinf_vec[1];
            boundaryState( 5 ) += Vinf_vec[2];
            break;}


        case(10):{ //Moon - parking orbit at 100 km
            boundaryState( 3 ) += Vinf_vec[0];
            boundaryState( 4 ) += Vinf_vec[1];
            boundaryState( 5 ) += Vinf_vec[2];
            break;}

        case(14):{ //Mars - parking orbit at 250 km
            boundaryState( 3 ) += Vinf_vec[0];
            boundaryState( 4 ) += Vinf_vec[1];
            boundaryState( 5 ) += Vinf_vec[2];
            break;}

        case(1):{ //EmL1 - stationary
            boundaryState = getLPGateway_DepState("Em", specificLPO, orbitFraction_, Time, stayTime_);
            break;}

        case(2):{ //EmL2 - stationary
            boundaryState = getLPGateway_DepState("Em", specificLPO, orbitFraction_, Time, stayTime_);
            break;}

        case(3):{ //EmL4 - stationary
            boundaryState = getLibrationState ("Em", 4, Time, true );
            break;}

        case(4):{ //EmL5 - stationary
            boundaryState = getLibrationState ("Em", 5, Time, true );
            break;}

        case(5):{ //SEL2 - stationary
            boundaryState = getLPGateway_DepState("SE", specificLPO, orbitFraction_, Time, stayTime_);
            break;}

        case(6):{ //SML1 - stationary
            boundaryState = getLPGateway_DepState("SM", specificLPO, orbitFraction_, Time, stayTime_);
            break;}

        case(7):{ //Earth orbiting
            boundaryState = getCentralGateway_DepState ("Earth", true, specificLPO, orbitFraction_, Time, stayTime_);
            break;}

        case(8):{ //Moon orbiting
            boundaryState = getCentralGateway_DepState ("Moon", true, specificLPO, orbitFraction_, Time, stayTime_);
            break;}

        case(9):{ //Mars orbiting
            boundaryState = getCentralGateway_DepState ("Mars", true, specificLPO, orbitFraction_, Time, stayTime_);
            break;}

    default: std::cout << "\nSwitcharg not recognized\n";


    }
    return boundaryState;
}

Eigen::Vector6d B_getArrivalState ( double Time, std::string trajectory, std::vector< ephemerides::EphemerisPointer > ephemerisVector, std::map<double, Eigen::VectorXd> specificLPO, double orbitFraction_, Eigen::Vector3d Vinf_vec){

    Eigen::Vector6d boundaryState = ephemerisVector[1]->getCartesianState(Time);

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char bound = trajectory.back();
    int switcharg{};
    if(bound=='E')switcharg=13;   if(bound=='m')switcharg=10;   if(bound=='M')switcharg=14;  //Earth & Moon & Mars
    if(bound=='1')switcharg=1;    if(bound=='2')switcharg=2;    if(bound=='3')switcharg=3;   //Lagrange points
    if(bound=='4')switcharg=4;    if(bound=='5')switcharg=5;    if(bound=='6')switcharg=6;   //Lagrange points
    if(bound=='7')switcharg=7;    if(bound=='8')switcharg=8;    if(bound=='9')switcharg=9;   //centralOrbits
    if(bound!='E'&& bound!='m' && bound!='M' && bound!='1'&&bound!='2'&&bound!='3'&&bound!='4'&&bound!='5'&&bound!='6'&&bound!='7'&&bound!='8'&&bound!='9')std::cout << "\nNo switcharg was generated!!\n";

    switch (switcharg) {
        case(13):{ //Earth - parking orbit at 250 km
            boundaryState( 3 ) += Vinf_vec[0];
            boundaryState( 4 ) += Vinf_vec[1];
            boundaryState( 5 ) += Vinf_vec[2];
            break;}


         case(10):{ //Moon - parking orbit at 100 km
            boundaryState( 3 ) += Vinf_vec[0];
            boundaryState( 4 ) += Vinf_vec[1];
            boundaryState( 5 ) += Vinf_vec[2];
            break;}

        case(14):{ //Mars - parking orbit at 250 km
            boundaryState( 3 ) += Vinf_vec[0];
            boundaryState( 4 ) += Vinf_vec[1];
            boundaryState( 5 ) += Vinf_vec[2];
            break;}

        case(1):{ //EmL1 - stationary
            boundaryState = getLPGateway_ArrState("Em", specificLPO, orbitFraction_, Time);
            break;}

        case(2):{ //EmL2 - stationary
            boundaryState = getLPGateway_ArrState("Em", specificLPO, orbitFraction_, Time);
            break;}

        case(3):{ //EmL4 - stationary
            boundaryState = getLibrationState ("Em", 4, Time, true );
            break;}

        case(4):{ //EmL5 - stationary
            boundaryState = getLibrationState ("Em", 5, Time, true );
            break;}

        case(5):{ //SEL2 - stationary
            boundaryState = getLPGateway_ArrState("SE", specificLPO, orbitFraction_, Time);
            break;}

        case(6):{ //SML1 - stationary
            boundaryState = getLPGateway_ArrState("SM", specificLPO, orbitFraction_, Time);
            break;}

        case(7):{ //Earth orbiting
            boundaryState = getCentralGateway_ArrState ("Earth", true, specificLPO, orbitFraction_, Time);
            break;}

        case(8):{ //Moon orbiting
            boundaryState = getCentralGateway_ArrState ("Moon", true, specificLPO, orbitFraction_, Time);
            break;}

        case(9):{ //Mars orbiting
            boundaryState = getCentralGateway_ArrState ("Mars", true, specificLPO, orbitFraction_, Time);
            break;}

    default: std::cout << "\nSwitcharg not recognized\n";


    }
    return boundaryState;
}


std::pair< double, double > getTOFbounds(std::string Gateway, std::string Segment){

    std::pair< double, double > TOFbounds = std::make_pair( 1.0, 365.0);

    bool Em_gateway {false};
    bool SE_gateway {false};
    bool SM_gateway {false};

    bool EG_segment {false};
    bool GM_segment {false};


    if(Gateway == "G1" || Gateway == "G2" || Gateway == "G3" || Gateway == "G4" || Gateway == "G7" || Gateway == "G8" ) {Em_gateway = true;}
    if(Gateway == "G5" ) {SE_gateway = true;}
    if(Gateway == "G6" ) {SM_gateway = true;}

    if(Segment == "EG" || Segment == "GE") {EG_segment = true;}
    if(Segment == "GM" || Segment == "MG") {GM_segment = true;}


    if((Em_gateway) && (EG_segment)) //Em gateways- segments EG and GE
    {  TOFbounds = std::make_pair( 1.0 , 25.0  );}

    if((SE_gateway) && (EG_segment)) //SE gateways- segments EG and GE
    { TOFbounds = std::make_pair( 1.0 , 100.0 );}

    if((SM_gateway) && (GM_segment)) //SM gateways- segments GM and MG
    {TOFbounds = std::make_pair( 1.0 , 100.0 );}

    return TOFbounds;

}

std::vector< int > getParameterIndices (std::vector<std::string> EM_trajectory ) {


    //Define segment trajectories
    std::string EG_trajectory = EM_trajectory.at(0);
    std::string GM_trajectory = EM_trajectory.at(1);


    bool EG_isCont {false};
    bool GM_isCont {false};
    if (EG_trajectory.find('c') != std::string::npos)
        EG_isCont = true;
    if (GM_trajectory.find('c') != std::string::npos)
        GM_isCont = true;

    int EG_numDSM = std::count(EG_trajectory.begin(), EG_trajectory.end(), 'd');
    int GM_numDSM = std::count(GM_trajectory.begin(), GM_trajectory.end(), 'd');

    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'd'), EG_trajectory.end());
    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'G'), EG_trajectory.end());
    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'c'), EG_trajectory.end());

    GM_trajectory.erase(std::remove(GM_trajectory.begin(), GM_trajectory.end(), 'd'), GM_trajectory.end());
    GM_trajectory.erase(std::remove(GM_trajectory.begin(), GM_trajectory.end(), 'G'), GM_trajectory.end());
    GM_trajectory.erase(std::remove(GM_trajectory.begin(), GM_trajectory.end(), 'c'), GM_trajectory.end());


    //Determine number of bounds to be created
    int numberBounds{0};
    int EG_Bounds{};
    int GM_Bounds{};
    int contBounds {7}; //number of Cont bounds
    int impLegBounds {1}; //number of Imp bounds per leg
    int impDSMBounds {4}; //number of Imp bounds per leg
    int EG_numLegs = EG_trajectory.size()-1;
    int GM_numLegs = GM_trajectory.size()-1;

    if(EG_isCont){
        EG_Bounds = contBounds;
        numberBounds = numberBounds + EG_Bounds;}
    else{
        EG_Bounds = ((EG_numLegs+1) * impLegBounds)+(EG_numDSM*impDSMBounds);
        numberBounds = numberBounds + EG_Bounds;
    }
    if(GM_isCont){
        GM_Bounds = contBounds;
        numberBounds = numberBounds + GM_Bounds;}
    else{
        GM_Bounds = ((GM_numLegs + 1)* impLegBounds)+(GM_numDSM*impDSMBounds);
        numberBounds = numberBounds + GM_Bounds;
    }

    std::vector<int> ParameterIndices = {numberBounds, EG_Bounds, GM_Bounds};

    return ParameterIndices;
}

Eigen::VectorXd getNsaveTrajectories (std::vector<std::string> EM_trajectory, std::vector<double> Individual, std::vector<int> ParameterIndices, std::string pathSuffix, std::unordered_map<int, std::map<double, Eigen::VectorXd>> orbitLibrary ) {

    Eigen::VectorXd ResultVec (14); //{fitness, deltaV, TOF, DV1, departureDate, arrivalDate, TOF1, DV2, stayTime, DepartureDate, arrivalDate, TOF2, orbitID, orbitFraction}

    //Declare variables to fill
    double fitnessvalue_EG{};
    double DV_EG{};
    Eigen::Vector3d Dates_EG;
    double EG_TOF{};

    double fitnessvalue_GM{};
    double DV_GM{};
    Eigen::Vector4d Dates_GM;
    double GM_TOF{};
    double GM_stayTime{0};


    //Define segment trajectories
    std::string EG_trajectory = EM_trajectory.at(0);
    std::string GM_trajectory = EM_trajectory.at(1);

    int EG_numDSM = std::count(EG_trajectory.begin(), EG_trajectory.end(), 'd');
    int GM_numDSM = std::count(GM_trajectory.begin(), GM_trajectory.end(), 'd');

    //Find out Transfer type
    bool EG_isCont {false};
    bool GM_isCont {false};
    if (EG_trajectory.find('c') != std::string::npos)
        EG_isCont = true;
    if (GM_trajectory.find('c') != std::string::npos)
        GM_isCont = true;

    size_t numberOfParameters_EG = ParameterIndices.at(1);

    //Some defnitions to use later on
    std::vector< ephemerides::EphemerisPointer > EG_ephemerisVector{};
    std::vector< ephemerides::EphemerisPointer > GM_ephemerisVector{};

    std::function< Eigen::Vector6d( const double ) > EG_initialStateFunction = [ = ]( const double currentTime )
    { return EG_ephemerisVector[0]->getCartesianState( currentTime ); };
    std::function< Eigen::Vector6d( const double ) > EG_finalStateFunction = [ = ]( const double currentTime )
    { return EG_ephemerisVector[1]->getCartesianState( currentTime ); };

    std::function< Eigen::Vector6d( const double ) > GM_initialStateFunction = [ = ]( const double currentTime )
    { return GM_ephemerisVector[0]->getCartesianState( currentTime ); };
    std::function< Eigen::Vector6d( const double ) > GM_finalStateFunction = [ = ]( const double currentTime )
    { return GM_ephemerisVector[1]->getCartesianState( currentTime ); };

    int numberOfRevolutions = 0;
    double sunGravitationalParameter = 1.32712428e20;
    double initialMass{1000.0};
    double specificImpulse = 1000;
    int numberOfSteps = 1000;

    int Ind_size = Individual.size();
    int lastElement = Ind_size - 1;
    int secondToLastElement = Ind_size - 2;

    //orbit modelling parameters
    int orbitID = round(Individual.at(secondToLastElement));
    double orbitFraction= Individual.at(lastElement);
    std::map<double, Eigen::VectorXd> specificOrbit = orbitLibrary[orbitID];
    double arrivalTime{};
    double departureTime{};


    ////////////////////////////  Segment 1 //////////////////////////////////////////
    std::string pathSuffix_EG = pathSuffix + "/" + EG_trajectory;
    std::function< std::vector< BaseFunctionVector >( const double ) > basisFunctionsFunction = std::bind( &getShapingBasisFunctions, std::placeholders::_1, numberOfRevolutions );
    std::vector<double> Individual_EG {Individual.begin(), Individual.begin()+numberOfParameters_EG};


    if(EG_isCont) { //Segment 1 is CT

        //create parameters without Individual
        std::vector< int > EG_flybysequence = createFlyBySeq(EG_trajectory);
        EG_ephemerisVector  = createEphemerisVec ( EG_flybysequence);

        //Create parameters using Individual
        double EG_departureTime = Individual_EG.at( 0 ) *  physical_constants::JULIAN_DAY;
        EG_TOF = Individual_EG.at( 1 ) *  physical_constants::JULIAN_DAY;
        double EG_arrivalTime = EG_departureTime + EG_TOF;
        arrivalTime = EG_arrivalTime;

        std::vector< BaseFunctionVector > basisFunctions = basisFunctionsFunction( EG_TOF );

        BaseFunctionVector EG_radialVelocityFunctionComponents = basisFunctions.at( 0 );
        BaseFunctionVector EG_normalVelocityFunctionComponents = basisFunctions.at( 1 );
        BaseFunctionVector EG_axialVelocityFunctionComponents = basisFunctions.at( 2 );

        int numberFreeCoefficientsRadialFunction = EG_radialVelocityFunctionComponents.size( ) - 3;
        int numberFreeCoefficientsNormalFunction = EG_normalVelocityFunctionComponents.size( ) - 3;
        int numberFreeCoefficientsAxialFunction = EG_axialVelocityFunctionComponents.size( ) - 3;

        Eigen::VectorXd EG_freeCoefficientsRadialVelocityFunction( numberFreeCoefficientsRadialFunction );
        Eigen::VectorXd EG_freeCoefficientsNormalVelocityFunction( numberFreeCoefficientsNormalFunction );
        Eigen::VectorXd EG_freeCoefficientsAxialVelocityFunction( numberFreeCoefficientsAxialFunction );

        for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
        {
            EG_freeCoefficientsRadialVelocityFunction[ i ] = Individual_EG[ i + 2 ];
        }
        for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
        {
            EG_freeCoefficientsNormalVelocityFunction[ i ] = Individual_EG[ i + 2 + numberFreeCoefficientsRadialFunction ];
        }
        for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
        {
            EG_freeCoefficientsAxialVelocityFunction[ i ] = Individual_EG[ i + 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
        }

        Eigen::Vector6d EG_departureState =  B_getDepartureState(EG_departureTime, EG_trajectory, EG_ephemerisVector, specificOrbit, orbitFraction, GM_stayTime );
        Eigen::Vector6d EG_arrivalState =  B_getArrivalState(EG_arrivalTime, EG_trajectory, EG_ephemerisVector, specificOrbit, orbitFraction );


        //Compute Trajectory
        tudat::shape_based_methods::HodographicShaping EG_CT = tudat::shape_based_methods::HodographicShaping(
                    EG_departureState, EG_arrivalState,
                    EG_TOF, sunGravitationalParameter,
                    numberOfRevolutions, EG_radialVelocityFunctionComponents,
                    EG_normalVelocityFunctionComponents, EG_axialVelocityFunctionComponents,
                    EG_freeCoefficientsRadialVelocityFunction, EG_freeCoefficientsNormalVelocityFunction,
                    EG_freeCoefficientsAxialVelocityFunction, initialMass );

        //Compute Trajectory parameters
        fitnessvalue_EG = EG_CT.computeDeltaV( );
        DV_EG = EG_CT.computeDeltaV( );
        Dates_EG << EG_departureTime /  physical_constants::JULIAN_DAY, EG_arrivalTime /  physical_constants::JULIAN_DAY, EG_TOF /  physical_constants::JULIAN_DAY;

        //Save files
        double stepSize = EG_TOF / static_cast< double >( numberOfSteps );

        std::vector< double > epochsToSaveResults;
        for ( int i = 0 ; i <= numberOfSteps ; i++ )
        {
            epochsToSaveResults.push_back( i * stepSize );
        }

        std::map< double, Eigen::Vector6d > hodographicShapingTrajectory;
        std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
        std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
        std::map< double, Eigen::VectorXd > hodographicShapingThrustAcceleration;

        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize );

        EG_CT.getTrajectory(
                    epochsToSaveResults, hodographicShapingTrajectory );

        EG_CT.getMassProfile(
                    epochsToSaveResults, hodographicShapingMassProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
        EG_CT.getThrustForceProfile(
                    epochsToSaveResults, hodographicShapingThrustProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
        EG_CT.getCylindricalThrustAccelerationProfile(
                    epochsToSaveResults, hodographicShapingThrustAcceleration );

        input_output::writeDataMapToTextFile(
                    hodographicShapingTrajectory, "hodographicShapingOptimalTrajectory.dat", tudat_applications::getOutputPath( )+pathSuffix_EG);

        input_output::writeDataMapToTextFile(
                    hodographicShapingMassProfile, "hodographicShapingOptimalMassProfile.dat", tudat_applications::getOutputPath( )+pathSuffix_EG);

        input_output::writeDataMapToTextFile(
                    hodographicShapingThrustProfile, "hodographicShapingOptimalThrustProfile.dat", tudat_applications::getOutputPath( )+pathSuffix_EG);

        input_output::writeDataMapToTextFile(
                    hodographicShapingThrustAcceleration, "hodographicShapingOptimalThrustAcceleration.dat", tudat_applications::getOutputPath( )+pathSuffix_EG);


    }
    else { //Segment 1 is IT

        //Create parameters without Individual
        std::vector< B_LegTypeTransfer > EG_LegTypeVector = createLegType_B(EG_trajectory);
        std::vector< int > EG_flybysequence = createFlyBySeq(EG_trajectory);
        int EG_numLegs = EG_flybysequence.size();
        std::vector< ephemerides::EphemerisPointer > EG_ephemerisVector ( EG_numLegs );
        Eigen::VectorXd EG_gravitationalParameterVector( EG_numLegs );
        Eigen::VectorXd EG_minimumPericenterRadii( EG_numLegs );
        EG_ephemerisVector = createEphemerisVec (EG_flybysequence);
        EG_gravitationalParameterVector=createGravParamVec(EG_flybysequence);
        if(EG_numDSM==0) {EG_minimumPericenterRadii=createMinPerRad(EG_flybysequence);}
        else {
            for (int i{0}; i < EG_numLegs; i++){
                EG_minimumPericenterRadii [i] = TUDAT_NAN;
            }
        }
        Eigen::Vector2d EG_semiMajorAxes = getSMA(EG_trajectory);
        Eigen::Vector2d EG_eccentricities = getEccentricities(EG_trajectory);

        //Create parameters using Individual
        Eigen::VectorXd EG_variableVector (EG_numLegs +1+(4*EG_numDSM));

        EG_TOF = 0;
        for(int i = 0; i < EG_numLegs ; i++){
            EG_variableVector[ i ] = Individual_EG [ i ];
            if( i > 0 ){
                EG_TOF += Individual_EG[i];
            }
        }

        double EG_departureTime = Individual_EG [0];
        double EG_arrivalTime = EG_departureTime  + EG_TOF;
        arrivalTime = EG_arrivalTime * physical_constants::JULIAN_DAY;
        EG_variableVector *= physical_constants::JULIAN_DAY;
        EG_variableVector[ EG_numLegs ] = 1;//dummy

        if(EG_numDSM>0){
            int count {EG_numLegs+1};
            for(size_t i{0}; i<EG_numDSM; i++){
                EG_variableVector[ count ] = Individual_EG[ count-1 ];
                EG_variableVector[ count + 1 ] = Individual_EG[ count ];
                EG_variableVector[ count + 2 ] = Individual_EG[ count + 1 ];
                EG_variableVector[ count + 3 ] = Individual_EG[ count + 2 ];
                count+=4;
            }
        }


        //Compute Trajectory
        B_makeTraj EG_IT ( EG_numLegs, EG_LegTypeVector, EG_trajectory, EG_ephemerisVector,
                              EG_gravitationalParameterVector, EG_variableVector, sunGravitationalParameter,
                              EG_minimumPericenterRadii, EG_semiMajorAxes, EG_eccentricities, specificOrbit, orbitFraction );

        //Compute Trajectory Parameters
        EG_IT.calculateTrajectory( fitnessvalue_EG );
        EG_IT.calculateTrajectory( DV_EG );
        Dates_EG << EG_departureTime, EG_arrivalTime, EG_TOF;


        //Save files
        // Define vectors to calculate intermediate points
        std::vector< Eigen::Vector3d > interPositionVectorTop;
        std::vector< double > interTimeVectorTop;

        // Calculate intermediate points
        EG_IT.intermediatePoints( numberOfSteps , interPositionVectorTop, interTimeVectorTop );

        // Define vectors to calculate manoeuvres
        std::vector< Eigen::Vector3d > manPositionVectorTop;
        std::vector< double > manTimeVectorTop;
        std::vector< double > manDeltaVVectorTop;

        // Calculate maneuvers
        EG_IT.maneuvers( manPositionVectorTop, manTimeVectorTop, manDeltaVVectorTop );

        //Writing intermediate points to file
         std::string outputFileTraj = tudat_applications::getOutputPath( ) + pathSuffix_EG + "/";
         tudat::transfer_trajectories::writeTrajectoryToFile_thesis( interPositionVectorTop, interTimeVectorTop, "TopTrajectory.dat", outputFileTraj );

        // Writing Maneuvres to file
         std::string outputFileMan = tudat_applications::getOutputPath( ) + pathSuffix_EG + "/";
         tudat::transfer_trajectories::writeTrajectoryToFile_thesis( manPositionVectorTop, manTimeVectorTop, "TopManeuvers.dat", outputFileMan );

        // Writing DeltaV's to file
        std::ofstream file(tudat_applications::getOutputPath( ) + pathSuffix_EG +  "/TopDeltaV.dat");
        for (size_t i {0}; i < manDeltaVVectorTop.size(); i++){
            file << i << " " <<manDeltaVVectorTop.at(i) << std::endl;
        }
        file.close();

    }




    ////////////////////////////  Segment 2 //////////////////////////////////////////
    std::string pathSuffix_GM = pathSuffix + "/" + GM_trajectory;

    std::vector<double> Individual_GM {Individual.begin() + numberOfParameters_EG, Individual.end()};


    if(GM_isCont) { //Segment 2 is CT

        //create parameters without Individual
        std::vector< int > GM_flybysequence = createFlyBySeq(GM_trajectory);
        GM_ephemerisVector  = createEphemerisVec ( GM_flybysequence);

        //Create parameters using Individual
        GM_stayTime = Individual_GM.at( 0 );
        double GM_departureTime = (Dates_EG [1] + GM_stayTime)*physical_constants::JULIAN_DAY;
        GM_TOF = Individual_GM.at( 1 ) * physical_constants::JULIAN_DAY;
        double GM_arrivalTime = GM_departureTime + GM_TOF;
        departureTime = GM_departureTime;

        std::vector< BaseFunctionVector > basisFunctions = basisFunctionsFunction( GM_TOF );


        BaseFunctionVector GM_radialVelocityFunctionComponents = basisFunctions.at( 0 );
        BaseFunctionVector GM_normalVelocityFunctionComponents = basisFunctions.at( 1 );
        BaseFunctionVector GM_axialVelocityFunctionComponents = basisFunctions.at( 2 );

        int numberFreeCoefficientsRadialFunction = GM_radialVelocityFunctionComponents.size( ) - 3;
        int numberFreeCoefficientsNormalFunction = GM_normalVelocityFunctionComponents.size( ) - 3;
        int numberFreeCoefficientsAxialFunction = GM_axialVelocityFunctionComponents.size( ) - 3;

        Eigen::VectorXd GM_freeCoefficientsRadialVelocityFunction( numberFreeCoefficientsRadialFunction );
        Eigen::VectorXd GM_freeCoefficientsNormalVelocityFunction( numberFreeCoefficientsNormalFunction );
        Eigen::VectorXd GM_freeCoefficientsAxialVelocityFunction( numberFreeCoefficientsAxialFunction );


        for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
        {
            GM_freeCoefficientsRadialVelocityFunction[ i ] = Individual_GM[ i + 2 ];
        }
        for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
        {
            GM_freeCoefficientsNormalVelocityFunction[ i ] = Individual_GM[ i + 2 + numberFreeCoefficientsRadialFunction ];
        }
        for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
        {
            GM_freeCoefficientsAxialVelocityFunction[ i ] = Individual_GM[ i + 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
        }

        Eigen::Vector6d GM_departureState =  B_getDepartureState(GM_departureTime, GM_trajectory, GM_ephemerisVector, specificOrbit, orbitFraction, GM_stayTime );
        Eigen::Vector6d GM_arrivalState =  B_getArrivalState(GM_arrivalTime, GM_trajectory, GM_ephemerisVector, specificOrbit, orbitFraction );


        //Compute Trajectory
        tudat::shape_based_methods::HodographicShaping GM_CT = tudat::shape_based_methods::HodographicShaping(
                    GM_departureState, GM_arrivalState,
                    GM_TOF, sunGravitationalParameter,
                    numberOfRevolutions, GM_radialVelocityFunctionComponents,
                    GM_normalVelocityFunctionComponents, GM_axialVelocityFunctionComponents,
                    GM_freeCoefficientsRadialVelocityFunction, GM_freeCoefficientsNormalVelocityFunction,
                    GM_freeCoefficientsAxialVelocityFunction, initialMass );


        //Compute Trajectory parameters
                fitnessvalue_GM = GM_CT.computeDeltaV( );
                DV_GM = GM_CT.computeDeltaV( );
                Dates_GM << GM_stayTime, GM_departureTime/ physical_constants::JULIAN_DAY, GM_arrivalTime/ physical_constants::JULIAN_DAY, GM_TOF / physical_constants::JULIAN_DAY;    //stayTime, DepartureDate, arrivalDate, TOF2

                //Save files
                double stepSize = GM_TOF / static_cast< double >( numberOfSteps );

                std::vector< double > epochsToSaveResults;
                for ( int i = 0 ; i <= numberOfSteps ; i++ )
                {
                    epochsToSaveResults.push_back( i * stepSize );
                }

                std::map< double, Eigen::Vector6d > hodographicShapingTrajectory;
                std::map< double, Eigen::VectorXd > hodographicShapingMassProfile;
                std::map< double, Eigen::VectorXd > hodographicShapingThrustProfile;
                std::map< double, Eigen::VectorXd > hodographicShapingThrustAcceleration;

                std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
                        std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize );

                GM_CT.getTrajectory(
                            epochsToSaveResults, hodographicShapingTrajectory );

                GM_CT.getMassProfile(
                            epochsToSaveResults, hodographicShapingMassProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
                GM_CT.getThrustForceProfile(
                            epochsToSaveResults, hodographicShapingThrustProfile, [ = ]( const double ){ return specificImpulse; }, integratorSettings );
                GM_CT.getCylindricalThrustAccelerationProfile(
                            epochsToSaveResults, hodographicShapingThrustAcceleration );

                input_output::writeDataMapToTextFile(
                            hodographicShapingTrajectory, "hodographicShapingOptimalTrajectory.dat", tudat_applications::getOutputPath( )+pathSuffix_GM);

                input_output::writeDataMapToTextFile(
                            hodographicShapingMassProfile, "hodographicShapingOptimalMassProfile.dat", tudat_applications::getOutputPath( )+pathSuffix_GM);

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustProfile, "hodographicShapingOptimalThrustProfile.dat", tudat_applications::getOutputPath( )+pathSuffix_GM);

                input_output::writeDataMapToTextFile(
                            hodographicShapingThrustAcceleration, "hodographicShapingOptimalThrustAcceleration.dat", tudat_applications::getOutputPath( )+pathSuffix_GM);

    }
    else { //Segment 2 is IT

        //Create parameters without Individual
        std::vector< B_LegTypeTransfer > GM_LegTypeVector = createLegType_B(GM_trajectory);
        std::vector< int > GM_flybysequence = createFlyBySeq(GM_trajectory);
        int GM_numLegs = GM_flybysequence.size();
        std::vector< ephemerides::EphemerisPointer > GM_ephemerisVector ( GM_numLegs );
        Eigen::VectorXd GM_gravitationalParameterVector( GM_numLegs );
        Eigen::VectorXd GM_minimumPericenterRadii( GM_numLegs );
        GM_ephemerisVector = createEphemerisVec (GM_flybysequence);
        GM_gravitationalParameterVector=createGravParamVec(GM_flybysequence);
        if(EG_numDSM==0) {GM_minimumPericenterRadii=createMinPerRad(GM_flybysequence);}
        else {
            for (int i{0}; i < GM_numLegs; i++){
                GM_minimumPericenterRadii [i] = TUDAT_NAN;
            }
        }
        Eigen::Vector2d GM_semiMajorAxes = getSMA(GM_trajectory);
        Eigen::Vector2d GM_eccentricities = getEccentricities(GM_trajectory);


        //Create parameters using Individual
        Eigen::VectorXd GM_variableVector (GM_numLegs +1+(4*GM_numDSM));

        GM_stayTime = Individual_GM.at( 0 );
        double GM_departureTime = Dates_EG [1] + GM_stayTime;

        GM_variableVector[ 0 ] = GM_departureTime; //Departure time = arrival time EG PLUS stay time (xv[numEGvars_])
        departureTime = GM_departureTime * physical_constants::JULIAN_DAY;

        GM_TOF = 0;
        for(int i = 1; i < GM_numLegs ; i++){
            GM_variableVector[ i ] = Individual_GM [ i ];
            if( i > 0 ){
                GM_TOF += Individual_GM[i];
            }
        }

        double GM_arrivalTime = GM_departureTime  + GM_TOF;
        GM_variableVector *= physical_constants::JULIAN_DAY;
        GM_variableVector[ GM_numLegs ] = 1;//dummy

        if(GM_numDSM>0){
            int count {GM_numLegs+1};
            for(size_t i{0}; i<GM_numDSM; i++){
                GM_variableVector[ count ] = Individual_GM[ count-1 ];
                GM_variableVector[ count + 1 ] = Individual_GM[ count ];
                GM_variableVector[ count + 2 ] = Individual_GM[ count + 1 ];
                GM_variableVector[ count + 3 ] = Individual_GM[ count + 2 ];
                count+=4;
            }
        }

        //Compute Trajectory
        B_makeTraj GM_IT ( GM_numLegs, GM_LegTypeVector, GM_trajectory, GM_ephemerisVector,
                              GM_gravitationalParameterVector, GM_variableVector, sunGravitationalParameter,
                              GM_minimumPericenterRadii, GM_semiMajorAxes, GM_eccentricities, specificOrbit, orbitFraction, GM_stayTime );

        //Compute Trajectory Parameters
        GM_IT.calculateTrajectory( fitnessvalue_GM );
        GM_IT.calculateTrajectory( DV_GM );
        Dates_GM << GM_stayTime, GM_departureTime, GM_arrivalTime, GM_TOF;


                //Save files
                // Define vectors to calculate intermediate points
                std::vector< Eigen::Vector3d > interPositionVectorTop;
                std::vector< double > interTimeVectorTop;

                // Calculate intermediate points
                GM_IT.intermediatePoints( numberOfSteps , interPositionVectorTop, interTimeVectorTop );

                // Define vectors to calculate manoeuvres
                std::vector< Eigen::Vector3d > manPositionVectorTop;
                std::vector< double > manTimeVectorTop;
                std::vector< double > manDeltaVVectorTop;

                // Calculate maneuvers
                GM_IT.maneuvers( manPositionVectorTop, manTimeVectorTop, manDeltaVVectorTop );

                //Writing intermediate points to file
                 std::string outputFileTraj = tudat_applications::getOutputPath( ) + pathSuffix_GM + "/";
                 tudat::transfer_trajectories::writeTrajectoryToFile_thesis( interPositionVectorTop, interTimeVectorTop, "TopTrajectory.dat", outputFileTraj );

                // Writing Maneuvres to file
                 std::string outputFileMan = tudat_applications::getOutputPath( ) + pathSuffix_GM + "/";
                 tudat::transfer_trajectories::writeTrajectoryToFile_thesis( manPositionVectorTop, manTimeVectorTop, "TopManeuvers.dat", outputFileMan );

                // Writing DeltaV's to file
                std::ofstream file(tudat_applications::getOutputPath( ) + pathSuffix_GM +  "/TopDeltaV.dat");
                for (size_t i {0}; i < manDeltaVVectorTop.size(); i++){
                    file << i << " " <<manDeltaVVectorTop.at(i) << std::endl;
                }
                file.close();
            }


    ////////////////////////////  Gateway //////////////////////////////////////////

    //Initiating parameters
    bool isLPO {false};
    std::string system {};
    std::string planet {};

    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'c'), EG_trajectory.end());
    char GW_char = EG_trajectory.back();
    int GW_int = GW_char - '0';

    switch(GW_int)
    {
    case( 1 ):{ //EmL1
        isLPO = true;
        system = "Em";
        break;}

    case( 2 ):{ //EmL2
        isLPO = true;
        system = "Em";
        break;}

    case( 3 ):{ //EmL4
        isLPO = true;
        system = "Em";
        break;}

    case( 4 ):{ //EmL5
        isLPO = true;
        system = "Em";
        break;}

    case( 5 ):{ //SEL2
        isLPO = true;
        system = "SE";
        break;}

    case( 6 ):{ //SML1
        isLPO = true;
        system = "SM";
        break;}

    case( 7 ):{ //Earth orbitting
        isLPO = false;
        planet = "Earth";
        break;}

    case( 8 ):{ //Moon orbitting
        isLPO = false;
        planet = "Moon";
        break;}

    case( 9 ):{ //Mars orbitting
        isLPO = false;
        planet = "Mars";
        break;}
    }

    std::string outputFileGateway = tudat_applications::getOutputPath( ) + pathSuffix + "/Gateway/";

    if(isLPO) {
        std::map<double,Eigen::Vector6d> Gateway_orbit_cart =  getLPGateway_orbit (system,  specificOrbit,  orbitFraction, arrivalTime, GM_stayTime, true );
        input_output::writeDataMapToTextFile(Gateway_orbit_cart, "Gateway_orbit_cart.dat", outputFileGateway );
        std::map<double,Eigen::Vector6d> Gateway_orbit_cr3b =  getLPGateway_orbit (system,  specificOrbit,  orbitFraction, arrivalTime, GM_stayTime, false );
        input_output::writeDataMapToTextFile(Gateway_orbit_cr3b, "Gateway_orbit_cr3b.dat", outputFileGateway );

        Eigen::Vector6d Gateway_ArrState_cart = getLPGateway_ArrState(system, specificOrbit, orbitFraction, arrivalTime, true );
        Eigen::Vector6d Gateway_DepState_cart = getLPGateway_DepState(system, specificOrbit, orbitFraction, departureTime, GM_stayTime, true);
        std::map<int, Eigen::Vector6d> Gateway_ArrDepStates_cart;
        Gateway_ArrDepStates_cart[1] << Gateway_ArrState_cart;
        Gateway_ArrDepStates_cart[2] << Gateway_DepState_cart;
        input_output::writeDataMapToTextFile(Gateway_ArrDepStates_cart, "Gateway_ArrDepStates_cart.dat", outputFileGateway );

        Eigen::Vector6d Gateway_ArrState_cr3b = getLPGateway_ArrState(system, specificOrbit, orbitFraction, arrivalTime, false );
        Eigen::Vector6d Gateway_DepState_cr3b = getLPGateway_DepState(system, specificOrbit, orbitFraction, departureTime, GM_stayTime, false);
        std::map<int, Eigen::Vector6d> Gateway_ArrDepStates_cr3b;
        Gateway_ArrDepStates_cr3b[1] << Gateway_ArrState_cr3b;
        Gateway_ArrDepStates_cr3b[2] << Gateway_DepState_cr3b;
        input_output::writeDataMapToTextFile(Gateway_ArrDepStates_cr3b, "Gateway_ArrDepStates_cr3b.dat", outputFileGateway );
    }
    else {

        //std::map<double,Eigen::Vector6d> Gateway_orbit_Sun =  getCentralGateway_orbit (planet,  specificOrbit,  orbitFraction, arrivalTime, GM_stayTime, true );
        //input_output::writeDataMapToTextFile(Gateway_orbit_Sun, "Gateway_orbit_Sun.dat", outputFileGateway );
        std::map<double,Eigen::Vector6d> Gateway_orbit_Planet =  getCentralGateway_orbit (planet,  specificOrbit,  orbitFraction, arrivalTime, GM_stayTime, false );
        input_output::writeDataMapToTextFile(Gateway_orbit_Planet, "Gateway_orbit_Planet.dat", outputFileGateway );

        Eigen::Vector6d Gateway_ArrState_Sun = getCentralGateway_ArrState(planet, true, specificOrbit, orbitFraction, arrivalTime );
        Eigen::Vector6d Gateway_DepState_Sun = getCentralGateway_DepState(planet, true, specificOrbit, orbitFraction, departureTime, GM_stayTime);
        std::map<int, Eigen::Vector6d> Gateway_ArrDepStates_Sun;
        Gateway_ArrDepStates_Sun[1] << Gateway_ArrState_Sun;
        Gateway_ArrDepStates_Sun[2] << Gateway_DepState_Sun;
        input_output::writeDataMapToTextFile(Gateway_ArrDepStates_Sun, "Gateway_ArrDepStates_Sun.dat", outputFileGateway );

        Eigen::Vector6d Gateway_ArrState_Planet = getCentralGateway_ArrState(planet, false, specificOrbit, orbitFraction, arrivalTime );
        Eigen::Vector6d Gateway_DepState_Planet = getCentralGateway_DepState(planet, false, specificOrbit, orbitFraction, departureTime, GM_stayTime);
        std::map<int, Eigen::Vector6d> Gateway_ArrDepStates_Planet;
        Gateway_ArrDepStates_Sun[1] << Gateway_ArrState_Planet;
        Gateway_ArrDepStates_Sun[2] << Gateway_DepState_Planet;
        input_output::writeDataMapToTextFile(Gateway_ArrDepStates_Planet, "Gateway_ArrDepStates_Planet.dat", outputFileGateway );

    }


    double TotalFitness = fitnessvalue_EG + fitnessvalue_GM;
    double TotalDV = DV_EG + DV_GM;
    double TotalTOF = Dates_EG[2] + Dates_GM[0] + Dates_GM[3];

    //{fitness, deltaV, TOF, DV1, departureDate, arrivalDate, TOF1, DV2, stayTime, DepartureDate, arrivalDate, TOF2}
    ResultVec << TotalFitness, TotalDV, TotalTOF, DV_EG, Dates_EG, DV_GM, Dates_GM, orbitID, orbitFraction;


    return ResultVec;

}


/*
Eigen::VectorXd getVVsfromInd (std::vector<std::string> EM_trajectory, std::vector<double> Individual, std::vector<int> ParameterIndices) {

    //Define segment trajectories
    std::string EG_trajectory = EM_trajectory.at(0);
    std::string GM_trajectory = EM_trajectory.at(1);


    bool EG_isCont {false};
    bool GM_isCont {false};
    if (EG_trajectory.find('c') != std::string::npos)
        EG_isCont = true;
    if (GM_trajectory.find('c') != std::string::npos)
        GM_isCont = true;



    //Get VV of segment1

    //get VV of segment2

    //if segment is ct, than return

    IDEE IS NU:
    in main: splits individual naar segment 1 en segment 2
    if (ct) -> functie SaveCTresults in map Segment
    if (it) -> functie SaveITresults in map Segment
    beide functies returnen ook vector met (fitness, departDate, tof, arrivalDate), deze kan gebruikt worden in final map.

}
*/





std::vector< std::vector< double > > B_createBounds (std::string Gateway, std::vector<std::string> EM_trajectory, double start_date, double launch_window, std::pair <double, double> Stay, std::pair<int, int> orbitID, bool returnFlight ) {

    double maxVinf = 5e3;

    std::vector<std::string> Segments = {"EG", "GM"};
    if(returnFlight){Segments = {"MG", "GE"};}

    //Define segment trajectories
    std::string EG_trajectory = EM_trajectory.at(0);
    std::string GM_trajectory = EM_trajectory.at(1);


    bool EG_isCont {false};
    bool GM_isCont {false};
    if (EG_trajectory.find('c') != std::string::npos)
        EG_isCont = true;
    if (GM_trajectory.find('c') != std::string::npos)
        GM_isCont = true;

    int EG_numDSM = std::count(EG_trajectory.begin(), EG_trajectory.end(), 'd');
    int GM_numDSM = std::count(GM_trajectory.begin(), GM_trajectory.end(), 'd');
     std::vector< LegTypeTransfer > EG_LegTypeVector;
     std::vector< int > EG_flybysequence;
     std::vector<int> EG_sequence;

     std::vector< LegTypeTransfer > GM_LegTypeVector;
     std::vector< int > GM_flybysequence;
     std::vector<int> GM_sequence;


    if(!EG_isCont) {
        EG_LegTypeVector = createLegType(EG_trajectory);
        EG_flybysequence = createFlyBySeq(EG_trajectory);
        EG_sequence = createFlyBySeq(EG_trajectory);

    }

    if(!GM_isCont) {
        GM_LegTypeVector = createLegType(GM_trajectory);
        GM_flybysequence = createFlyBySeq(GM_trajectory);
        GM_sequence = createFlyBySeq(GM_trajectory);
    }

    //Erase d, G and c
    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'd'), EG_trajectory.end());
    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'G'), EG_trajectory.end());
    EG_trajectory.erase(std::remove(EG_trajectory.begin(), EG_trajectory.end(), 'c'), EG_trajectory.end());


    GM_trajectory.erase(std::remove(GM_trajectory.begin(), GM_trajectory.end(), 'd'), GM_trajectory.end());
    GM_trajectory.erase(std::remove(GM_trajectory.begin(), GM_trajectory.end(), 'G'), GM_trajectory.end());
    GM_trajectory.erase(std::remove(GM_trajectory.begin(), GM_trajectory.end(), 'c'), GM_trajectory.end());


    //Determine number of bounds to be created
    int numberTrajectoryBounds{0};
    int numberOrbitBounds{2};
    int totalNumberBounds{0};

    int EG_Bounds{};
    int GM_Bounds{};
    int contBounds {10}; //number of Cont bounds
    int impLegBounds {1}; //number of Imp bounds per leg
    int impDSMBounds {4}; //number of Imp bounds per leg
    int EG_numLegs = EG_trajectory.size()-1;
    int GM_numLegs = GM_trajectory.size()-1;

    if(EG_isCont){
        EG_Bounds = contBounds;
        numberTrajectoryBounds = numberTrajectoryBounds + EG_Bounds;}
    else{
        EG_Bounds = ((EG_numLegs+1) * impLegBounds)+(EG_numDSM*impDSMBounds);
        numberTrajectoryBounds = numberTrajectoryBounds + EG_Bounds;
    }
    if(GM_isCont){
        GM_Bounds = contBounds;
        numberTrajectoryBounds = numberTrajectoryBounds + GM_Bounds;}
    else{
        GM_Bounds = ((GM_numLegs + 1)* impLegBounds)+(GM_numDSM*impDSMBounds);
        numberTrajectoryBounds = numberTrajectoryBounds + GM_Bounds;
    }

    totalNumberBounds = numberTrajectoryBounds + numberOrbitBounds;


    //create bounds variable
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( totalNumberBounds, 0.0 ) );

    if(EG_isCont){ //create Continuous Thrust bounds for first segment
     std::pair< double, double > TOFbound = getTOFbounds(Gateway, Segments.at(0) );

        // Only Dep orbit
        bounds[ 0 ][ 0 ] =  start_date;                 //1a
        bounds[ 1 ][ 0 ] = start_date + launch_window;  //1b
        bounds[ 0 ][ 1 ] = TOFbound.first;              //2a
        bounds[ 1 ][ 1 ] = TOFbound.second;             //2b
        bounds[ 0 ][ 2 ] = - 10000.0;                   //3a
        bounds[ 1 ][ 2 ] = 10000.0;                     //3b
        bounds[ 0 ][ 3 ] = -10000;                      //4a
        bounds[ 1 ][ 3 ] = 10000.0;                     //4b
        bounds[ 0 ][ 4 ] = - 10000.0;                   //5a
        bounds[ 1 ][ 4 ] = 10000.0;                     //5b
        bounds[ 0 ][ 5 ] = -10000;                      //6a
        bounds[ 1 ][ 5 ] = 10000.0;                     //6b
        bounds[ 0 ][ 6 ] = -10000;                      //7a
        bounds[ 1 ][ 6 ] = 10000.0;                     //7b
        bounds[ 0 ][ 7 ] = -1 * maxVinf;                //8a Dep Vinf_x
        bounds[ 1 ][ 7 ] = 1 * maxVinf;                 //8b
        bounds[ 0 ][ 8 ] = -1 * maxVinf;                //9a Dep Vinf_y
        bounds[ 1 ][ 8 ] = 1 * maxVinf;                 //9b
        bounds[ 0 ][ 9 ] = -1 * maxVinf;                //10a Dep Vinf_z
        bounds[ 1 ][ 9 ] = 1 * maxVinf;                 //10b

    }

    else {
        //create Impulsive Thrust bounds for first segment
        bounds[ 0 ][ 0 ] = start_date; //MJD2000
        bounds[ 1 ][ 0 ] = start_date + launch_window; //MJD2000

        // create traditional ToF Bounds
        enum LegOptions{ EE, Em, EM, mE, mm, mM, ME, Mm, MM};
        int casepointer {};

        std::stringstream EG_sequenceSS;
        std::copy(EG_sequence.begin(), EG_sequence.end(), std::ostream_iterator<int>(EG_sequenceSS));
        std::string EG_sequence_String = EG_sequenceSS.str();



        for (int i {0}; i < EG_numLegs; i++){
            int countLeg{i+1};
            std::string legString {EG_sequence_String.substr(i,2)};


            if( legString == "33" ) casepointer = EE;
            if( legString == "30" ) casepointer = Em;
            if( legString == "34" ) casepointer = EM;

            if( legString == "03" ) casepointer = mE;
            if( legString == "00" ) casepointer = mm;
            if( legString == "04" ) casepointer = mM;

            if( legString == "43" ) casepointer = ME;
            if( legString == "40" ) casepointer = Mm;
            if( legString == "44" ) casepointer = MM;

            switch(casepointer){
            case (EE):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (Em):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (EM):
                bounds[ 0 ][ countLeg ] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;

            case (mE):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (mm):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (mM):
                bounds[ 0 ][ countLeg ] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;

            case (ME):
                bounds[ 0 ][ countLeg] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;
            case (Mm):
                bounds[ 0 ][ countLeg] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;
            case (MM):
                bounds[ 0 ][ countLeg] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            default: std::cout << "Leg input could not be recognized\n";
            }

           }


        // create additional DSM Bounds
        int counter{EG_numLegs+1};
        std::vector <int> DSMpos;
         for(int i{0}; i < EG_LegTypeVector.size(); i++){
             if(EG_LegTypeVector.at(i)==6) {DSMpos.push_back(i);}
              }

        if(EG_LegTypeVector[0]==5){
            bounds[ 0 ][ counter ] = 0.01;                                                                       //1st additional var - fraction of DSM
            bounds[ 1 ][ counter ] = 0.99;
            bounds[ 0 ][ counter+1 ] = 0;                                                                     //2nd additional var - V_inf at departure
            bounds[ 1 ][ counter+1 ] = 5000;
            bounds[ 0 ][ counter+2 ] = 0 * 2 * mathematical_constants::getPi< double >( );                       //3rd additional var - in-plane angle for the hyperbolic excess velocity
            bounds[ 1 ][ counter+2 ] = 1 * 2 * mathematical_constants::getPi< double >( );
            bounds[ 0 ][ counter+3 ] = std::acos(  2 * 1 - 1. ) - mathematical_constants::getPi< double >( )/2; //4th additional var - out-of-plane angle for the hyperbolic excess velocity
            bounds[ 1 ][ counter+3 ] = std::acos(  2 * 0 - 1. ) - mathematical_constants::getPi< double >( )/2;
            counter+=4;

            for (int i{1}; i<EG_numDSM;i++){
                int pos = DSMpos.at(i-1);
                std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(EG_flybysequence.at(pos));
                Eigen::VectorXd minPerRad = createMinPerRad(SwingbyPLanet);
                Eigen::VectorXd maxPerRad = createMaxPerRad(SwingbyPLanet);

                bounds[ 0 ][ counter ] = 0.01;                                                 //1st additional var - fraction of DSM
                bounds[ 1 ][ counter ] = 0.99;
                bounds[ 0 ][ counter+1 ] = -1 * mathematical_constants::getPi< double >( ); //2nd additional var - rotation angle of the GA
                bounds[ 1 ][ counter+1 ] = 1 * mathematical_constants::getPi< double >();
                bounds[ 0 ][ counter+2 ] = minPerRad[0];                                    //3rd additional var - pericenter radius of the GA
                bounds[ 1 ][ counter+2 ] = maxPerRad[0];
                bounds[ 0 ][ counter+3 ] = 0; //4th additional var - V added for the powered GA
                bounds[ 1 ][ counter+3 ] = 2e3;
                counter+=4;
            }
        }

        if(EG_LegTypeVector[0]!=5){

        for (int i{0}; i<EG_numDSM;i++){
            int pos = DSMpos.at(i-1);
            std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(EG_flybysequence.at(pos));
            Eigen::VectorXd minPerRad = createMinPerRad(SwingbyPLanet);
            Eigen::VectorXd maxPerRad = createMaxPerRad(SwingbyPLanet);

            bounds[ 0 ][ counter ] = 0.01;                                                 //1st additional var - fraction of DSM
            bounds[ 1 ][ counter ] = 0.99;
            bounds[ 0 ][ counter+1 ] = -1 * mathematical_constants::getPi< double >( ); //2nd additional var - rotation angle of the GA
            bounds[ 1 ][ counter+1 ] = 1 * mathematical_constants::getPi< double >();
            bounds[ 0 ][ counter+2 ] = minPerRad[0];                                    //3rd additional var - pericenter radius of the GA
            bounds[ 1 ][ counter+2 ] = maxPerRad[0];
            bounds[ 0 ][ counter+3 ] = 0; //4th additional var - V added for the powered GA
            bounds[ 1 ][ counter+3 ] = 2e3;
            counter+=4;
        }
        }

    }

    // Second segment

    if(GM_isCont){ //create Continuous Thrust bounds for first segment
     std::pair< double, double > TOFbound = getTOFbounds(Gateway, Segments.at(1) );

        // Only Dep orbit
        bounds[ 0 ][ EG_Bounds ] =  Stay.first;                     //1a
        bounds[ 1 ][ EG_Bounds ] = Stay.second;                     //1b
        bounds[ 0 ][ EG_Bounds + 1 ] = TOFbound.first;              //2a
        bounds[ 1 ][ EG_Bounds + 1 ] = TOFbound.second;             //2b
        bounds[ 0 ][ EG_Bounds + 2 ] = - 10000.0;                   //3a
        bounds[ 1 ][ EG_Bounds + 2 ] = 10000.0;                     //3b
        bounds[ 0 ][ EG_Bounds + 3 ] = -10000;                      //4a
        bounds[ 1 ][ EG_Bounds + 3 ] = 10000.0;                     //4b
        bounds[ 0 ][ EG_Bounds + 4 ] = - 10000.0;                   //5a
        bounds[ 1 ][ EG_Bounds + 4 ] = 10000.0;                     //5b
        bounds[ 0 ][ EG_Bounds + 5 ] = -10000;                      //6a
        bounds[ 1 ][ EG_Bounds + 5 ] = 10000.0;                     //6b
        bounds[ 0 ][ EG_Bounds + 6 ] = -10000;                      //7a
        bounds[ 1 ][ EG_Bounds + 6 ] = 10000.0;                     //7b
        bounds[ 0 ][ EG_Bounds + 7 ] = -1 * maxVinf;                //8a Dep Vinf_x
        bounds[ 1 ][ EG_Bounds + 7 ] = 1 * maxVinf;                 //8b
        bounds[ 0 ][ EG_Bounds + 8 ] = -1 * maxVinf;                //9a Dep Vinf_y
        bounds[ 1 ][ EG_Bounds + 8 ] = 1 * maxVinf;                 //9b
        bounds[ 0 ][ EG_Bounds + 9 ] = -1 * maxVinf;                //10a Dep Vinf_z
        bounds[ 1 ][ EG_Bounds + 9 ] = 1 * maxVinf;                 //10b

    }
    else { //create Impulsive Thrust bounds for first segment

        bounds[ 0 ][ EG_Bounds  ] = Stay.first; //MJD2000
        bounds[ 1 ][ EG_Bounds ] = Stay.second; //MJD2000

        // create traditional ToF Bounds
        enum LegOptions{ EE, Em, EM, mE, mm, mM, ME, Mm, MM};
        int casepointer {};

        std::stringstream GM_sequenceSS;
        std::copy(GM_sequence.begin(), GM_sequence.end(), std::ostream_iterator<int>(GM_sequenceSS));
        std::string GM_sequence_String = GM_sequenceSS.str();


        for (int i {0}; i < GM_numLegs; i++){
            int countLeg{EG_Bounds + i+1};

            std::string legString {GM_sequence_String.substr(i,2)};


            if( legString == "33" ) casepointer = EE;
            if( legString == "30" ) casepointer = Em;
            if( legString == "34" ) casepointer = EM;

            if( legString == "03" ) casepointer = mE;
            if( legString == "00" ) casepointer = mm;
            if( legString == "04" ) casepointer = mM;

            if( legString == "43" ) casepointer = ME;
            if( legString == "40" ) casepointer = Mm;
            if( legString == "44" ) casepointer = MM;

            switch(casepointer){
            case (EE):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (Em):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (EM):
                bounds[ 0 ][ countLeg ] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;

            case (mE):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (mm):
                bounds[ 0 ][ countLeg ] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            case (mM):
                bounds[ 0 ][ countLeg ] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;

            case (ME):
                bounds[ 0 ][ countLeg] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;
            case (Mm):
                bounds[ 0 ][ countLeg] = 100;
                bounds[ 1 ][ countLeg ] = 365;
                break;
            case (MM):
                bounds[ 0 ][ countLeg] = 1;
                bounds[ 1 ][ countLeg ] = 100;
                break;
            default: std::cout << "Leg input could not be recognized\n";
            }

           }

        // create additional DSM Bounds
        int counter{EG_Bounds + GM_numLegs+1};
        std::vector <int> DSMpos;
         for(int i{0}; i < GM_LegTypeVector.size(); i++){
             if(GM_LegTypeVector.at(i)==6) {DSMpos.push_back(i);}
              }

        if(GM_LegTypeVector[0]==5){
            bounds[ 0 ][ counter ] = 0.01;                                                                       //1st additional var - fraction of DSM
            bounds[ 1 ][ counter ] = 0.99;
            bounds[ 0 ][ counter+1 ] = 0;                                                                     //2nd additional var - V_inf at departure
            bounds[ 1 ][ counter+1 ] = 5000;
            bounds[ 0 ][ counter+2 ] = 0 * 2 * mathematical_constants::getPi< double >( );                       //3rd additional var - in-plane angle for the hyperbolic excess velocity
            bounds[ 1 ][ counter+2 ] = 1 * 2 * mathematical_constants::getPi< double >( );
            bounds[ 0 ][ counter+3 ] = std::acos(  2 * 1 - 1. ) - mathematical_constants::getPi< double >( )/2; //4th additional var - out-of-plane angle for the hyperbolic excess velocity
            bounds[ 1 ][ counter+3 ] = std::acos(  2 * 0 - 1. ) - mathematical_constants::getPi< double >( )/2;
            counter+=4;

            for (int i{1}; i<GM_numDSM;i++){
                int pos = DSMpos.at(i-1);
                std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(GM_flybysequence.at(pos));
                Eigen::VectorXd minPerRad = createMinPerRad(SwingbyPLanet);
                Eigen::VectorXd maxPerRad = createMaxPerRad(SwingbyPLanet);

                bounds[ 0 ][ counter ] = 0.01;                                                 //1st additional var - fraction of DSM
                bounds[ 1 ][ counter ] = 0.99;
                bounds[ 0 ][ counter+1 ] = -1 * mathematical_constants::getPi< double >( ); //2nd additional var - rotation angle of the GA
                bounds[ 1 ][ counter+1 ] = 1 * mathematical_constants::getPi< double >();
                bounds[ 0 ][ counter+2 ] = minPerRad[0];                                    //3rd additional var - pericenter radius of the GA
                bounds[ 1 ][ counter+2 ] = maxPerRad[0];
                bounds[ 0 ][ counter+3 ] = 0; //4th additional var - V added for the powered GA
                bounds[ 1 ][ counter+3 ] = 2e3;
                counter+=4;
            }


        }
        if(GM_LegTypeVector[0]!=5){

        for (int i{0}; i<GM_numDSM;i++){
            int pos = DSMpos.at(i-1);
            std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(GM_flybysequence.at(pos));
            Eigen::VectorXd minPerRad = createMinPerRad(SwingbyPLanet);
            Eigen::VectorXd maxPerRad = createMaxPerRad(SwingbyPLanet);

            bounds[ 0 ][ counter ] = 0.01;                                                 //1st additional var - fraction of DSM
            bounds[ 1 ][ counter ] = 0.99;
            bounds[ 0 ][ counter+1 ] = -1 * mathematical_constants::getPi< double >( ); //2nd additional var - rotation angle of the GA
            bounds[ 1 ][ counter+1 ] = 1 * mathematical_constants::getPi< double >();
            bounds[ 0 ][ counter+2 ] = minPerRad[0];                                    //3rd additional var - pericenter radius of the GA
            bounds[ 1 ][ counter+2 ] = maxPerRad[0];
            bounds[ 0 ][ counter+3 ] = 0; //4th additional var - V added for the powered GA
            bounds[ 1 ][ counter+3 ] = 2e3;
            counter+=4;
        }
        }



    }

    // Add bounds needed to create orbit modelling parameters
    bounds[ 0 ][ numberTrajectoryBounds ] = orbitID.first; //orbitID
    bounds[ 1 ][ numberTrajectoryBounds ] = orbitID.second;
    bounds[ 0 ][ numberTrajectoryBounds + 1] = 0; //orbitFraction
    bounds[ 1 ][ numberTrajectoryBounds + 1 ] = 1.0;


    return bounds;

}


std::map<double, Eigen::VectorXd> getLPO_cr3b (std::string system, std::string LPO_type, std::string LP,  int orbitID, bool save) {

    std::string filePath_prefix = tudat_applications::getOutputPath( ) + "LPO-Library/";
    std::string fileName = LP + "_" + LPO_type + "_" + std::to_string(orbitID) + ".txt";
    std::string filePath = filePath_prefix + system + "/" + LPO_type + "/" + LP + "/";

    std::vector<double> raw_cr3b{};
    double element{};

    std::ifstream myfile (filePath + fileName);
    if (myfile.is_open())
    {
      while ( !myfile.eof() )
          {
          myfile >> element;
          raw_cr3b.push_back(element);
          }

      myfile.close();
    }
    else {
        std::cout << "\nError in getLPO_cr3b  -  filename not found:   " << filePath + fileName << std::endl;
    }

    std::map<double, Eigen::VectorXd> cr3bMap;
    int loopSize = (raw_cr3b.size()-7);
    for (int i{0}; i < loopSize; i+=7 ) {
        double key = raw_cr3b.at(i);
        Eigen::Vector6d value; value <<  raw_cr3b.at(i+1), raw_cr3b.at(i+2), raw_cr3b.at(i+3), raw_cr3b.at(i+4), raw_cr3b.at(i+5), raw_cr3b.at(i+6);
        cr3bMap[key] = value;
    }

    if(save) {

        std::string outputPath = filePath + "CR3B-Maps/";
        std::string outputFileName = std::to_string(orbitID) + ".dat";
        input_output::writeDataMapToTextFile(cr3bMap, outputFileName, outputPath );

    }


    return cr3bMap;

}
/*
std::map<double, Eigen::Vector6d> getLPO_cart (std::string system, std::map<double, Eigen::Vector6d> specificLPO, double initialTime, int numberRevs) {

    //Initializing parameters
    double primGrav{};
    double secGrav{};
    double initialDistance{};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> P1P2ephemeris {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> SunEarthephemeris {};


    if(system == "SE") {
        //distance
         std::vector<int> SunEB = {3};
         P1P2ephemeris = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true);
         initialDistance = P1P2ephemeris[0]->getCartesianState( initialTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
    }
    else if (system == "Em") {
        //distance
        std::vector<int> EarthMoon = {0};
        std::vector<int> SunEarth = {3};

        P1P2ephemeris = createEphemerisVec(EarthMoon, planet, "ECLIPJ2000", false);
        SunEarthephemeris = createEphemerisVec(SunEarth, "Sun", "ECLIPJ2000", false);

        initialDistance = P1P2ephemeris[0]->getCartesianState( initialTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties(planet, "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        }

    else if (system == "SM") {
        //distance
        std::vector<int> SunMars = {4};
        P1P2ephemeris = createEphemerisVec(SunMars, "Sun", "ECLIPJ2000", false);
        initialDistance = P1P2ephemeris[0]->getCartesianState( initialTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;
        }
    else {
        std::cout << "\nError - system not recognized in FUNCTION: getLPO_cart. \n";
    }


    std::map<double, Eigen::Vector6d> cr3bMap = specificLPO;

    double LPO_period = getLPO_period (system, cr3bMap, initialTime) * physical_constants::JULIAN_DAY;

    //create iterator
    std::map<double, Eigen::Vector6d>::iterator it = cr3bMap.begin();

    //Initiate map
    std::map<double, Eigen::Vector6d> cartMap{};
    Eigen::Vector6d dimensionalCartState {};

    for (int i {0}; i < (numberRevs + 1); i++ ) {

        while (it != cr3bMap.end())
        {
            //converting dimensionless time to dimensional time using initialDistance between primary and secondary
            double dimensionalTime = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(it->first, primGrav, secGrav, initialDistance );

            double adjustedDimensionalTime = dimensionalTime + initialTime + (i * LPO_period);

            Eigen::Vector6d P1P2state = P1P2ephemeris[0]->getCartesianState( adjustedDimensionalTime );


            if(system == "Em") {
                Eigen::Vector6d dimensionalCartState_EarthFrame = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, it->second, P1P2state);

                dimensionalCartState = dimensionalCartState_EarthFrame + SunEarthephemeris[0]->getCartesianState( adjustedDimensionalTime );

            }
            else {
                dimensionalCartState = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, it->second, P1P2state);
            }

            cartMap[adjustedDimensionalTime] = dimensionalCartState;

            it++;
        }

        it = cr3bMap.begin();

    }

    return cartMap;

}


double getLPO_period (std::string system, std::map<double, Eigen::Vector6d> specificLPO, double arrivalTime) {

    //Initializing parameters
    double primGrav{};
    double secGrav{};
    double initialDistance{};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> P1P2ephemeris {};


    if(system == "SE") {
        //distance
         std::vector<int> SunEB = {3};
         P1P2ephemeris = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true);
         initialDistance = P1P2ephemeris[0]->getCartesianState( arrivalTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
    }
    else if (system == "Em") {
        //distance
        std::vector<int> EarthMoon = {0};
        P1P2ephemeris = createEphemerisVec(EarthMoon, planet, "ECLIPJ2000", false);
        initialDistance = P1P2ephemeris[0]->getCartesianState( arrivalTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties(planet, "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        }

    else if (system == "SM") {
        //distance
        std::vector<int> SunMars = {4};
        P1P2ephemeris = createEphemerisVec(SunMars, "Sun", "ECLIPJ2000", false);
        initialDistance = P1P2ephemeris[0]->getCartesianState( arrivalTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;
        }
    else {
        std::cout << "\nError - system not recognized in FUNCTION: getLPO_cart. \n";
    }

    std::map<double, Eigen::Vector6d> cr3bMap = specificLPO;

    double dimensionalTime = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(cr3bMap.rbegin()->first, primGrav, secGrav, initialDistance );

    return dimensionalTime / physical_constants::JULIAN_DAY; //return period in days

}


Eigen::Vector6d getLPO_ArrState (std::string system, std::map<double, Eigen::Vector6d> specificLPO,  double fraction, double arrivalTime) {


    int numRev = 1;

    std::map<double, Eigen::Vector6d> LPOcart_map =  getLPO_cart (system, specificLPO, arrivalTime, numRev);

    int numberOfTimeSteps = LPOcart_map.size();

    int arrivalTimeStep = round(fraction * numberOfTimeSteps );

    if(arrivalTimeStep == numberOfTimeSteps) arrivalTimeStep--;

    std::map<double, Eigen::Vector6d>::iterator it = LPOcart_map.begin();

    std::advance(it, arrivalTimeStep);

    Eigen::Vector6d LPO_ArrState = it->second;


     //Eigen::Vector6d LPO_ArrState = getLibrationState ("Em", 1, arrivalTime, true);

    return LPO_ArrState;


}

Eigen::Vector6d getLPO_DepState (std::string system, std::map<double, Eigen::Vector6d> specificLPO, double fraction, double departureTime, double stayTime) {


    double arrivalTime = departureTime - (stayTime * physical_constants::JULIAN_DAY);

    double LPO_period = getLPO_period (system, specificLPO, arrivalTime);

    double stayPeriod = stayTime / LPO_period;

    int numRev = round(stayPeriod) + 1;

    std::map<double, Eigen::Vector6d> LPOcart_map =  getLPO_cart (system, specificLPO, arrivalTime, numRev);

    int numberOfTimeSteps_perOrbit = LPOcart_map.size()/numRev;

    int numberOfTimeSteps = LPOcart_map.size();

    int arrivalTimeStep = round(fraction * numberOfTimeSteps_perOrbit );

    int departureTimeStep = arrivalTimeStep + round( stayPeriod * numberOfTimeSteps_perOrbit);

    if(departureTimeStep == numberOfTimeSteps) departureTimeStep--;


    std::map<double, Eigen::Vector6d>::iterator it = LPOcart_map.begin();

    std::advance(it, departureTimeStep);

    Eigen::Vector6d LPO_DepState = it->second;



    //Eigen::Vector6d LPO_DepState = getLibrationState ("Em", 1, departureTime, true);


    return LPO_DepState;

}
*/


std::vector< B_LegTypeTransfer > createLegType_B (std::string flyby_seq){
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'G'), flyby_seq.end());
    size_t numLegs = flyby_seq.length( )-1;
    std::vector< B_LegTypeTransfer > legTypeVector;
    size_t count{0};
    while (count <= numLegs ){
        char next{};
        if((count)<numLegs){next=flyby_seq.at(count+1);}
        if(count==0 && next=='d') {legTypeVector.push_back(B_MGA1DsmVelocity_Departure); count=count+2; continue;}
        if(count==0 && next!='d'){legTypeVector.push_back(B_MGA_Departure);count++;continue;}
        if(count!=0 && next=='d'){legTypeVector.push_back(B_MGA1DsmVelocity_Swingby); count=count+2; continue;}
        if(count!=0 && count!=numLegs && next!='d'){legTypeVector.push_back(B_MGA_Swingby);count++;continue;}
        if(count==numLegs){legTypeVector.push_back(B_Capture);count++;}
    }
return legTypeVector;
}

std::unordered_map<int, std::map<double, Eigen::VectorXd>> importLPO_cr3b (std::string system, std::string LPO_type, std::string LP, std::pair<int, int> orbitID_bounds) {
    std::unordered_map<int, std::map<double, Eigen::VectorXd>> allLPOs{};

    for(int i{orbitID_bounds.first}; i < (orbitID_bounds.second + 1); i++ ) {
        std::map<double, Eigen::VectorXd> singleLPO{};
        singleLPO = getLPO_cr3b(system, LPO_type, LP, i, false);
        allLPOs[i] = singleLPO;
    }

    return allLPOs;

}


Eigen::Vector6d getLPGateway_ArrState (std::string system, std::map<double, Eigen::VectorXd> Gateway_LPO, double arrivalFraction, double arrivalTime, bool useRefFrame) {
    Eigen::VectorXd Gateway_ArrState_Xd{};
    Eigen::Vector6d Gateway_ArrState{};


    //Initializing parameters
    double primGrav{};
    double secGrav{};
    double averageDistance{};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> P1P2ephemeris {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> SunEarthephemeris {};


    if(system == "SE") {
        //distance
         std::vector<int> SunEB = {3};
         P1P2ephemeris = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true);
         averageDistance = physical_constants::ASTRONOMICAL_UNIT; // average distance Sun - Earth

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
    }
    else if (system == "Em") {
        //distance
        std::vector<int> EarthMoon = {0};
        std::vector<int> SunEarth = {3};

        P1P2ephemeris = createEphemerisVec(EarthMoon, "Earth", "ECLIPJ2000", false);
        SunEarthephemeris = createEphemerisVec(SunEarth, "Sun", "ECLIPJ2000", false);

        averageDistance = 384400000; // average distance Earth - Moon

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        }

    else if (system == "SM") {
        //distance
        std::vector<int> SunMars = {4};
        P1P2ephemeris = createEphemerisVec(SunMars, "Sun", "ECLIPJ2000", false);
        averageDistance = 1.524 * physical_constants::ASTRONOMICAL_UNIT; // average distance Sun - Mars (taken from NASA factsheet)

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;
        }
    else {
        std::cout << "\nError - system not recognized in FUNCTION: getLPO_cart. \n";
    }

    //Get Gateway period in days
    std::map<double, Eigen::VectorXd>::reverse_iterator  it_last = Gateway_LPO.rbegin();
    double LP_period = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(it_last->first, primGrav, secGrav, averageDistance ) / physical_constants::JULIAN_DAY;

    std::map<double, Eigen::VectorXd>::iterator it_arr = Gateway_LPO.begin();
    int timeStepofInterest = round(arrivalFraction * Gateway_LPO.size());
    int advancement = (timeStepofInterest!=0) ? (timeStepofInterest-1) : timeStepofInterest;
    std::advance(it_arr, advancement);
    Eigen::VectorXd state_cr3b = it_arr->second;
    if(!useRefFrame) {
        return convertXdtoSixd(state_cr3b);
    }

    if (system == "Em") {
        Eigen::VectorXd state_cart_E = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, state_cr3b, P1P2ephemeris[0]->getCartesianState(arrivalTime));
        Gateway_ArrState_Xd = state_cart_E + SunEarthephemeris[0]->getCartesianState(arrivalTime);
    }
    else {
        Gateway_ArrState_Xd = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, state_cr3b, P1P2ephemeris[0]->getCartesianState(arrivalTime));

    }

    Gateway_ArrState = convertXdtoSixd(Gateway_ArrState_Xd);

    return Gateway_ArrState;
}


Eigen::Vector6d getLPGateway_DepState (std::string system, std::map<double, Eigen::VectorXd> Gateway_LPO, double arrivalFraction, double departureTime, double stayTime, bool useRefFrame ) {


    Eigen::VectorXd Gateway_DepState_Xd{};
    Eigen::Vector6d Gateway_DepState{};


    //Initializing parameters
    double primGrav{};
    double secGrav{};
    double averageDistance{};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> P1P2ephemeris {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> SunEarthephemeris {};


    if(system == "SE") {
        //distance
         std::vector<int> SunEB = {3};
         P1P2ephemeris = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true);
         averageDistance = physical_constants::ASTRONOMICAL_UNIT; // average distance Sun - Earth

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
    }
    else if (system == "Em") {
        //distance
        std::vector<int> EarthMoon = {0};
        std::vector<int> SunEarth = {3};

        P1P2ephemeris = createEphemerisVec(EarthMoon, "Earth", "ECLIPJ2000", false);
        SunEarthephemeris = createEphemerisVec(SunEarth, "Sun", "ECLIPJ2000", false);

        averageDistance = 384400000; // average distance Earth - Moon

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        }

    else if (system == "SM") {
        //distance
        std::vector<int> SunMars = {4};
        P1P2ephemeris = createEphemerisVec(SunMars, "Sun", "ECLIPJ2000", false);
        averageDistance = 1.524 * physical_constants::ASTRONOMICAL_UNIT; // average distance Sun - Mars (taken from NASA factsheet)

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;
        }
    else {
        std::cout << "\nError - system not recognized in FUNCTION: getLPO_cart. \n";
    }


    //Get Gateway period in days
    std::map<double, Eigen::VectorXd>::reverse_iterator  it_last = Gateway_LPO.rbegin();
    double LP_period = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(it_last->first, primGrav, secGrav, averageDistance ) / physical_constants::JULIAN_DAY;


    double DepFraction = (arrivalFraction + (stayTime/LP_period));
            while(DepFraction>1.0) { DepFraction = DepFraction-1; }

    int timeStepofInterest_Dep = round(DepFraction * Gateway_LPO.size());

    int advancement_Dep = (timeStepofInterest_Dep!=0) ? (timeStepofInterest_Dep-1) : timeStepofInterest_Dep;
    std::map<double, Eigen::VectorXd>::iterator it_dep = Gateway_LPO.begin();
    std::advance(it_dep, advancement_Dep);
    Eigen::VectorXd stateDep_cr3b = it_dep->second;
    if(!useRefFrame) {
        return convertXdtoSixd(stateDep_cr3b);
    }

    if (system == "Em") {
        Eigen::VectorXd state_cart_E = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, stateDep_cr3b, P1P2ephemeris[0]->getCartesianState(departureTime));
        Gateway_DepState_Xd = state_cart_E + SunEarthephemeris[0]->getCartesianState(departureTime);
    }
    else {
        Gateway_DepState_Xd = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, stateDep_cr3b, P1P2ephemeris[0]->getCartesianState(departureTime));

    }

    Gateway_DepState = convertXdtoSixd(Gateway_DepState_Xd);

    return Gateway_DepState;

}


std::map<double,Eigen::Vector6d> getLPGateway_orbit (std::string system, std::map<double, Eigen::VectorXd> Gateway_LPO, double arrivalFraction, double arrivalTime, double stayTime, bool useRefFrame) {
    std::map<double,Eigen::Vector6d> Gateway_orbit{};


    //std::cout << "\n*********  Researching functioning of getLPGateway_orbit  *************\n\n ";

    //Initializing parameters
    double primGrav{};
    double secGrav{};
    double averageDistance{};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> P1P2ephemeris {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> SunEarthephemeris {};


    if(system == "SE") {
        //distance
         std::vector<int> SunEB = {3};
         P1P2ephemeris = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true);
         averageDistance = physical_constants::ASTRONOMICAL_UNIT; // average distance Sun - Earth

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
    }
    else if (system == "Em") {
        //distance
        std::vector<int> EarthMoon = {0};
        std::vector<int> SunEarth = {3};

        P1P2ephemeris = createEphemerisVec(EarthMoon, "Earth", "ECLIPJ2000", false);
        SunEarthephemeris = createEphemerisVec(SunEarth, "Sun", "ECLIPJ2000", false);

        averageDistance = 384400000; // average distance Earth - Moon

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        }

    else if (system == "SM") {
        //distance
        std::vector<int> SunMars = {4};
        P1P2ephemeris = createEphemerisVec(SunMars, "Sun", "ECLIPJ2000", false);
        averageDistance = 1.524 * physical_constants::ASTRONOMICAL_UNIT; // average distance Sun - Mars (taken from NASA factsheet)

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;
        }
    else {
        std::cout << "\nError - system not recognized in FUNCTION: getLPO_cart. \n";
    }

    //Get Gateway period in days
    std::map<double, Eigen::VectorXd>::reverse_iterator  it_last = Gateway_LPO.rbegin();
    double LP_period = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(it_last->first, primGrav, secGrav, averageDistance );

    int numberRevs{};
    if (useRefFrame) {numberRevs = stayTime/(LP_period/physical_constants::JULIAN_DAY) + 1;}
    else {numberRevs = 0;}

    //std::cout << "\n---LP_period: " << LP_period / physical_constants::JULIAN_DAY;

    std::map<double, Eigen::VectorXd>::iterator it_arr = Gateway_LPO.begin();
    int timeStepofInterest = round(arrivalFraction * Gateway_LPO.size());
    int advancement = (timeStepofInterest!=0) ? (timeStepofInterest-1) : timeStepofInterest;
    std::advance(it_arr, advancement);
    double timeStepOfArr_cr3b = it_arr->first;
    double timeStepOfArr_cart = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(timeStepOfArr_cr3b, primGrav, secGrav, averageDistance );
    double startTimeOrbit = arrivalTime - timeStepOfArr_cart;

    double key {};
    Eigen::Vector6d value{};

    for (int i {0}; i < (numberRevs + 1); i++ ) {


        for ( std::map<double, Eigen::VectorXd>::iterator it { Gateway_LPO.begin() }; it != Gateway_LPO.end(); it++) {
            if(useRefFrame){
                key = startTimeOrbit + tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(it->first, primGrav, secGrav, averageDistance ) + (i * LP_period);
                value = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, it->second, P1P2ephemeris[0]->getCartesianState(key));
                if(system=="Em"){
                    value = value + SunEarthephemeris[0]->getCartesianState(key);
                }
                Gateway_orbit[key] = value;
              }
            else {
                key = it->first;
                value = it->second;
                Gateway_orbit[key] = value;
            }

        }
    }

    return Gateway_orbit;


}

std::map<double,Eigen::Vector6d> getBody_orbit (std::string Body, double startTime, std::string centralBody, int numRev) {

    std::map<double,Eigen::Vector6d> BodyOrbit{};

    std::vector<int> BodyVec {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> BodyEphemeris {};
    double bodyPeriod{};


    if(Body == "Earth" || Body == "earth") {
        BodyVec = {3};
        BodyEphemeris = createEphemerisVec(BodyVec, centralBody, "ECLIPJ2000", true);
        bodyPeriod = 365.25 * physical_constants::JULIAN_DAY;
    }
    if(Body == "Moon" || Body == "moon") {
        BodyVec = {0};
        BodyEphemeris = createEphemerisVec(BodyVec, centralBody, "ECLIPJ2000", false);
        bodyPeriod = 27 * physical_constants::JULIAN_DAY;

    }
    if(Body == "Mars" || Body == "mars") {
        BodyVec = {4};
        BodyEphemeris = createEphemerisVec(BodyVec, centralBody, "ECLIPJ2000", true);
        bodyPeriod = 687 * physical_constants::JULIAN_DAY;
    }

    double numberStepsPerPeriod = 10000;
    double stepSize = bodyPeriod / numberStepsPerPeriod;

    for (double time {startTime}; time < (startTime + (bodyPeriod * numRev)); time = time + stepSize) {
        BodyOrbit[time] = BodyEphemeris[0]->getCartesianState(time);
    }

    return BodyOrbit;

}

std::unordered_map<int, std::map<double, Eigen::VectorXd>> generateCentralOrbits (std::string planet, std::pair<double, double> SMA_bounds, std::pair<double, double> eccentricity_bounds, std::pair<double, double> AOP_bounds, int Nsteps, std::string pathDirectory ){
    std::unordered_map<int, std::map<double, Eigen::VectorXd>> AllCentralOrbits {};
    std::map<int, Eigen::Vector3d> orbitID_map{};
    int orbitID = 0;
    const double SMA_min = SMA_bounds.first;
    const double SMA_max = SMA_bounds.second;
    const double Ecc_min = eccentricity_bounds.first;
    const double Ecc_max = eccentricity_bounds.second;
    const double AOP_min = AOP_bounds.first;
    const double AOP_max = AOP_bounds.second;

    std::string outputpreffix = tudat_applications::getOutputPath( ) + "/Orbits/CentralOrbits/" + planet + "/";

    ////////////////////////////////////////////////////////// Set variables //////////////////////////////////////////////////////////////////////////////////

    //Integration variables
    const double startTime = 0;
    const double N_timeSteps = 1000;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const double Ecc_stepSize = (Ecc_max - Ecc_min)/Nsteps;
    const double SMA_stepSize = (SMA_max - SMA_min)/Nsteps;
    const double AOP_stepSize = (AOP_max - AOP_min)/Nsteps;


    for (int Ecc_iter{0}; Ecc_iter < Nsteps; Ecc_iter++) {

        double Ecc = Ecc_min + Ecc_iter * Ecc_stepSize;

        for (int SMA_iter{0}; SMA_iter < Nsteps; SMA_iter++){

            double SMA = SMA_min + SMA_iter * SMA_stepSize;

            for (int AOP_iter{0}; AOP_iter < Nsteps; AOP_iter++) {

                double AOP = AOP_min + AOP_iter * AOP_stepSize;

                orbitID_map[orbitID] << Ecc, SMA, AOP;

            ////////////////////////////////////////////////////////// Set initial state of Gateway //////////////////////////////////////////////////////////////////

            // Initial gateway state in Kepler elements.
            Eigen::Vector6d initialGatewayState_kepler;
            initialGatewayState_kepler( semiMajorAxisIndex ) = SMA;
            initialGatewayState_kepler( eccentricityIndex ) = Ecc;
            initialGatewayState_kepler( inclinationIndex ) = tudat::unit_conversions::convertDegreesToRadians( 0 );
            initialGatewayState_kepler( argumentOfPeriapsisIndex ) = tudat::unit_conversions::convertDegreesToRadians( AOP );
            initialGatewayState_kepler( longitudeOfAscendingNodeIndex ) = tudat::unit_conversions::convertDegreesToRadians( 0 );
            initialGatewayState_kepler( trueAnomalyIndex ) = tudat::unit_conversions::convertDegreesToRadians( 0 );

            /////////////////////////////  Create environment  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Create body objects.
            std::vector< std::string > bodiesToCreate;
            bodiesToCreate.push_back( planet );
            std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                    getDefaultBodySettings( bodiesToCreate );
            bodySettings[ planet ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                        Eigen::Vector6d::Zero( ) );

            // Create Earth object
            NamedBodyMap bodyMap = createBodies( bodySettings );

            //Create gateway
            bodyMap[ "Gateway" ] = std::make_shared< simulation_setup::Body >( );

            //Set frame
            setGlobalFrameBodyEphemerides( bodyMap, planet, "ECLIPJ2000" );

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            bodiesToPropagate.push_back( "Gateway" );
            centralBodies.push_back( planet );

            // Define propagation settings.
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > Gateway_acceleration;
            Gateway_acceleration[ planet ].push_back( std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
            accelerationMap[ "Gateway" ] = Gateway_acceleration;

            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Convert Gateway state from Keplerian elements to Cartesian elements.
            double centralGravitationalParameter = bodyMap.at( planet )->getGravityFieldModel( )->getGravitationalParameter( );
            Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                        initialGatewayState_kepler,
                        centralGravitationalParameter );

            // Set simulation parameters.
            double simulationEndEpoch = 2 * mathematical_constants::PI * sqrt(pow(SMA,3)/centralGravitationalParameter);

            // Create propagator settings.
            std::shared_ptr< tudat::propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< tudat::propagators::TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );

            // Create numerical integrator settings.
            double simulationStartEpoch = startTime;
            const double fixedStepSize = simulationEndEpoch / N_timeSteps;
            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, fixedStepSize );


            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Create simulation object and propagate dynamics.
            tudat::propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            AllCentralOrbits[orbitID] = integrationResult;

            ////////////////////////////////////////////////////////// Save result ///////////////////////////////////////////////////////////////////////////////////////
            input_output::writeDataMapToTextFile( integrationResult,
                                                  std::to_string(orbitID) + ".dat",
                                                  outputpreffix,
                                                  "",
                                                  std::numeric_limits< double >::digits10,
                                                  std::numeric_limits< double >::digits10,
                                                  "," );


            //std::cout << "\n-- saved orbitID - " << orbitID;
            orbitID++;
            input_output::writeDataMapToTextFile(orbitID_map, "orbitID_map.dat", outputpreffix );
            input_output::writeDataMapToTextFile(orbitID_map, "orbitID_map.dat", pathDirectory );


            }
        }
    }

return AllCentralOrbits;

}

Eigen::Vector6d getCentralGateway_ArrState (std::string planet, bool useRefFrame, std::map<double, Eigen::VectorXd> centralOrbit, double arrivalFraction, double arrivalTime) {
    Eigen::VectorXd Gateway_ArrState_Xd{};

    std::vector<int> BodyVec {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> BodyEphemeris {};


    if(planet == "Earth" || planet == "earth") {
        BodyVec = {3};
        BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", true);
    }
    if(planet == "Moon" || planet == "moon") {
        BodyVec = {0};
        BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", false);

    }
    if(planet == "Mars" || planet == "mars") {
        BodyVec = {4};
        BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", true);
    }

    std::map<double, Eigen::VectorXd>::iterator it_arr = centralOrbit.begin();
    int timeStepofInterest = round(arrivalFraction * centralOrbit.size());
    int advancement = (timeStepofInterest!=0) ? (timeStepofInterest-1) : timeStepofInterest;
    std::advance(it_arr, advancement);

    Eigen::VectorXd ArrState_PlanetFrame = it_arr->second;
    Eigen::Vector6d PlanetState_SunFrame = BodyEphemeris[0]->getCartesianState(arrivalTime);

    if(useRefFrame){
    Gateway_ArrState_Xd = ArrState_PlanetFrame + PlanetState_SunFrame;
    }
    else {
        Gateway_ArrState_Xd = ArrState_PlanetFrame;
    }

    return convertXdtoSixd(Gateway_ArrState_Xd);

}

Eigen::Vector6d getCentralGateway_DepState (std::string planet, bool useRefFrame, std::map<double, Eigen::VectorXd> centralOrbit, double arrivalFraction, double departureTime, double stayTime) {
    Eigen::VectorXd Gateway_DepState_Xd{};

    std::vector<int> BodyVec {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> BodyEphemeris {};


    if(planet == "Earth" || planet == "earth") {
        BodyVec = {3};
        BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", true);
    }
    if(planet == "Moon" || planet == "moon") {
        BodyVec = {0};
        BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", false);

    }
    if(planet == "Mars" || planet == "mars") {
        BodyVec = {4};
        BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", true);
    }

    std::map<double, Eigen::VectorXd>::reverse_iterator  it_last = centralOrbit.rbegin();
    double centralOrbit_period = (it_last->first)/physical_constants::JULIAN_DAY;


    std::map<double, Eigen::VectorXd>::iterator it_dep = centralOrbit.begin();
    double departureFraction = (arrivalFraction + (stayTime/centralOrbit_period));
    while(departureFraction>1.0) { departureFraction = departureFraction-1; }

    int timeStepofInterest = round(departureFraction * centralOrbit.size());
    int advancement = (timeStepofInterest!=0) ? (timeStepofInterest-1) : timeStepofInterest;
    std::advance(it_dep, advancement);

    Eigen::Vector6d DepState_PlanetFrame = it_dep->second;
    Eigen::Vector6d PlanetState_SunFrame = BodyEphemeris[0]->getCartesianState(departureTime);

    if(useRefFrame){
        Gateway_DepState_Xd = DepState_PlanetFrame + PlanetState_SunFrame;
    }
    else {
        Gateway_DepState_Xd = DepState_PlanetFrame;
    }

    return convertXdtoSixd(Gateway_DepState_Xd);


}

std::map<double,Eigen::Vector6d> getCentralGateway_orbit (std::string planet, std::map<double, Eigen::VectorXd> centralOrbit, double arrivalFraction, double arrivalTime, double stayTime, bool useRefFrame ) {
    std::map<double,Eigen::Vector6d> centralGateway_orbit{};


    std::map<double, Eigen::VectorXd>::reverse_iterator  it_last = centralOrbit.rbegin();
    double centralOrbit_period = it_last->first; //in seconds

    int numRev{};
    if(!useRefFrame) {
        numRev = 1;
    }
    else {
        numRev = stayTime/(centralOrbit_period/physical_constants::JULIAN_DAY) + 1;
    }

    std::vector<int> BodyVec {};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> BodyEphemeris {};

    if(useRefFrame){

        if(planet == "Earth" || planet == "earth") {
            BodyVec = {3};
            BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", true);
        }
        if(planet == "Moon" || planet == "moon") {
            BodyVec = {0};
            BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", false);

        }
        if(planet == "Mars" || planet == "mars") {
            BodyVec = {4};
            BodyEphemeris = createEphemerisVec(BodyVec, "Sun", "ECLIPJ2000", true);
        }
    }

        std::map<double, Eigen::VectorXd>::iterator it_arr = centralOrbit.begin();
        int timeStepofInterest = round(arrivalFraction * centralOrbit.size());
        int advancement = (timeStepofInterest!=0) ? (timeStepofInterest-1) : timeStepofInterest;
        std::advance(it_arr, advancement);
        double arrivalTime_orbit = it_arr->first;
        double startTime_orbit = arrivalTime - arrivalTime_orbit;

        for(int Rev{0}; Rev < numRev; Rev++) {
            for(std::map<double, Eigen::VectorXd>::iterator iterator {centralOrbit.begin()}; iterator != centralOrbit.end(); iterator++ )
            {
                double time = startTime_orbit + iterator->first + (Rev * centralOrbit_period);
                Eigen::VectorXd gwState_planetframe = iterator->second;

                if(useRefFrame){
                    Eigen::VectorXd planetState_Sunframe = BodyEphemeris[0]->getCartesianState(time);
                    centralGateway_orbit[time] = convertXdtoSixd(gwState_planetframe + planetState_Sunframe);
                }
                else {
                    centralGateway_orbit[time] = convertXdtoSixd(gwState_planetframe);
                }

            }
        }


return centralGateway_orbit;

}


Eigen::Vector6d convertXdtoSixd (Eigen::VectorXd Xd_vector) {

    Eigen::Vector6d Sixd_vector {};
    Sixd_vector << Xd_vector[0], Xd_vector[1], Xd_vector[2], Xd_vector[3], Xd_vector[4], Xd_vector[5];

    return Sixd_vector;
}

std::unordered_map<int, std::map<double, Eigen::VectorXd>> createOrbitLibrary (std::string Gateway, std::map<std::string, std::pair<double, double>> boundsMap, std::string LPO_type, int Nsteps, std::string pathDirectory) {
    std::unordered_map<int, std::map<double, Eigen::VectorXd>> orbitLibrary {};

    //Getting LP number from Gateway string
    Gateway.erase(std::remove(Gateway.begin(), Gateway.end(), 'G'), Gateway.end());
    int LP_int = std::stoi( Gateway );
    std::string LP = "L" + std::to_string(LP_int);

    //Initiating parameters
    bool isLPO {false};
    std::string system {};
    std::string planet {};

    switch(LP_int)
    {
    case( 1 ):{ //EmL1
        isLPO = true;
        system = "Em";
        break;}

    case( 2 ):{ //EmL2
        isLPO = true;
        system = "Em";
        break;}

    case( 3 ):{ //EmL4
        isLPO = true;
        system = "Em";
        break;}

    case( 4 ):{ //EmL5
        isLPO = true;
        system = "Em";
        break;}

    case( 5 ):{ //SEL2
        isLPO = true;
        system = "SE";
        break;}

    case( 6 ):{ //SML1
        isLPO = true;
        system = "SM";
        break;}

    case( 7 ):{ //Earth orbitting
        isLPO = false;
        planet = "Earth";
        break;}

    case( 8 ):{ //Moon orbitting
        isLPO = false;
        planet = "Moon";
        break;}

    case( 9 ):{ //Mars orbitting
        isLPO = false;
        planet = "Mars";
        break;}
    }

    if(isLPO){
                std::pair<double, double> orbitID_bounds = boundsMap["orbitID"];
        orbitLibrary = importLPO_cr3b (system, LPO_type, LP, orbitID_bounds);
    }
    else {
        std::pair<double, double> SMA_bounds = boundsMap["SMA"];
        std::pair<double, double> eccentricity_bounds = boundsMap["Ecc"];
        std::pair<double, double> AOP_bounds = boundsMap["AOP"];
        orbitLibrary = generateCentralOrbits (planet, SMA_bounds, eccentricity_bounds, AOP_bounds, Nsteps, pathDirectory );
    }

    std::cout << "\n----\nOrbit Library created\n---\n ";

    return orbitLibrary;
}

std::pair<double, double> getOrbitIDbounds(int Nsteps, std::string Gateway, std::string LPO_type, int lower, int upper) {

    //parameters
    std::pair<double, double> haloG1 = std::make_pair(0, 2561);
    std::pair<double, double> haloG2 = std::make_pair(0, 1615);
    std::pair<double, double> horizontalG1 = std::make_pair(0, 1176);
    std::pair<double, double> horizontalG2 = std::make_pair(0, 1514);
    std::pair<double, double> verticalG1 = std::make_pair(0, 6590);
    std::pair<double, double> verticalG2 = std::make_pair(0, 306);

    //if set manually
    if(lower!=0 && upper!=0) {
        return std::make_pair(lower, upper);
    }

    //if set automatically
    Gateway.erase(std::remove(Gateway.begin(), Gateway.end(), 'G'), Gateway.end());
    int GW_int = std::stoi( Gateway );

    switch(GW_int)
    {
    case( 1 ):
    { //EmL1 gateway
        if(LPO_type=="halo") return haloG1;
        if(LPO_type=="horizontal") return horizontalG1;
        if(LPO_type=="vertical") return verticalG1;
        break;}

    case( 2 ):
    { //EmL2 gateway
        if(LPO_type=="halo") return haloG2;
        if(LPO_type=="horizontal") return horizontalG2;
        if(LPO_type=="vertical") return verticalG2;
        break;}

    case( 3 ):
    case( 4 ):
    case( 5 ):
    case( 6 ):
    { //other LP gateways
        std::cout << "\nLP orbits for this gateway has not yet been implemented (error in getOrbitIDbounds)";
        break;}

    case( 7 ):
    case( 8 ):
    case( 9 ):
    { //centralOrbit gateways

        int upper = Nsteps * Nsteps * Nsteps;
        return std::make_pair(lower, upper);
        break;}

    default: {

        std::cout << "\nError: Gateway not found (getOrbitIDbounds)";
        break;
    }
    }

    return std::make_pair(lower, upper);


}


std::vector<std::string> getTrajectoryOptions (std::string Gateway, std::string Segment) {

    std::vector<std::string> TrajectoryOptions{};

    Gateway.erase(std::remove(Gateway.begin(), Gateway.end(), 'G'), Gateway.end());
    int GW_int = std::stoi( Gateway );

    switch(GW_int)
    {
    case( 1 ): {
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG1", "EdmdG1", "EmG1", "EG1c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G1dmdM", "G1mEM", "G1dmdEdM", "G1mM",  "G1Mc"};
        }
        break;
    }
    case( 2 ): {
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG2", "EdmdG2", "EmG2", "EG2c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G2dmdEdMM", "G2dM", "G2mEM", "G2Mc"};
        }
        break;
    }
    case( 5 ): {
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG5", "EdmdG5", "EmG5", "EmmG5", "EG5c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G5dM", "G5dEdM", "G5mEM", "G5dmdEdM", "G5dEdmdEdM", "G5dmdEdEM", "G5Mc"};
        }
        break;
    }
    case( 6 ): {
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG6", "EdMdG6", "EMG6", "EdmdMdG6", "EmMG6", "EG6c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G6dM"};
        }
        break;
    }
    case( 7 ): { //centralOrbit Earth
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG7", "EG7",  "EG7c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G7dM", "G7dEdM", "G7dmdM", "G7dEdmdM", "G7Mc"};
        }
        break;
    }
    case( 8 ): { //centralOrbit Moon
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG8", "EdEdG8", "EdmdG8"  "EG8c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G8mM", "G8dmdM", "G8mEM", "G8dmdEdM", "G8Mc"};
        }
        break;
    }
    case( 9 ): { //centralOrbit Mars
        if(Segment=="EG"){
        TrajectoryOptions = {"EdG9", "EdMdG9", "EMG9", "EdmdMdG9", "EmMG9", "EG9c" };
        }
        if(Segment=="GM"){
        TrajectoryOptions = {"G9dM", "G9M", "G9Mc"};
        }
        break;
    }
    default: {
        std::cout << "\nError: Gateway-Segment combination is not available for Analysis B (getTrajectoryOptions)";
    }
    }

    return TrajectoryOptions;


}











//std::make_shared<std::unordered_map<int, std::map<double, Eigen::VectorXd>>>




