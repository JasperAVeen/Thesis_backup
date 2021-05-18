#include "tt_functions.h"
#include <sstream>
#include<iostream>
#include<iomanip>
#include<vector>
#include<cctype>
#include<cstring>
#include<string>
#include <fstream>
#include <algorithm>
#include "lagrangepoints.h"

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
using namespace tudat;



std::vector< int > createFlyBySeq (std::string flyby_seq){
    std::vector< int > flybySequence;
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'd'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'G'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'c'), flyby_seq.end());

    size_t numLegs = flyby_seq.length();

   // std::cout << "Length of flybyseq is: " << numLegs << std::endl;
    std::string planetLet  {"xVEmMJSUNP123456789"};
    std::vector<int> planetNum  {1, 2, 3, 0, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 3, 4, 3, 0, 4 };
    for (size_t i {0}; i < numLegs; i++){
        char ser = flyby_seq.at(i);
        size_t position{};
        position = planetLet.find(ser);
        int value {planetNum.at(position)};
        flybySequence.push_back( value );

    }
    return flybySequence;
}

Eigen::VectorXd createMinPerRad (std::vector< int > flybySequence){
    size_t numberOfLegs_{flybySequence.size()};
    Eigen::VectorXd minimumPericenterRadii_;
    minimumPericenterRadii_.resize( numberOfLegs_ );
    //double factor{1.04 * 1e3};

for(int i = 0; i < numberOfLegs_; i++)
{
    switch(flybySequence[ i ])
    {
    case( 1 ):{ //1.082
        minimumPericenterRadii_[ i ] = 1.082 * 1e3 *(spice_interface::getBodyProperties("MERCURY", "RADII",3).at(0)+spice_interface::getBodyProperties("MERCURY", "RADII",3).at(1)+spice_interface::getBodyProperties("MERCURY", "RADII",3).at(2))/3;//2639.7E3;
        break;}
    case( 2 ):{//1.047
        minimumPericenterRadii_[ i ] = 1.047 * 1e3 *(spice_interface::getBodyProperties("VENUS", "RADII",3).at(0)+spice_interface::getBodyProperties("VENUS", "RADII",3).at(1)+spice_interface::getBodyProperties("VENUS", "RADII",3).at(2))/3;//6251.8E3;
        break;}
    case( 3 ): { //1.048
        minimumPericenterRadii_[ i ] = 1.048 * 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;//6578.1E3;
        break;}
    case( 4 ):{
        minimumPericenterRadii_[ i ] = 1.076*1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;//1395.0E3;
        break;}
    case( 5 ):{
        minimumPericenterRadii_[ i ] = 1.6*1e3*(spice_interface::getBodyProperties("JUPITER", "RADII",3).at(0)+spice_interface::getBodyProperties("JUPITER", "RADII",3).at(1)+spice_interface::getBodyProperties("JUPITER", "RADII",3).at(2))/3;//3596.2E3;
        break;}
    case( 6 ):{
        minimumPericenterRadii_[ i ] = 1.342*1e3*(spice_interface::getBodyProperties("SATURN", "RADII",3).at(0)+spice_interface::getBodyProperties("SATURN", "RADII",3).at(1)+spice_interface::getBodyProperties("SATURN", "RADII",3).at(2))/3;//72000E3;
        break;}
    case( 7 ):{
        minimumPericenterRadii_[ i ] = 4.190*1e3*(spice_interface::getBodyProperties("URANUS", "RADII",3).at(0)+spice_interface::getBodyProperties("URANUS", "RADII",3).at(1)+spice_interface::getBodyProperties("URANUS", "RADII",3).at(2))/3;//61000.0E3;
        break;}
    case( 8 ):{
        minimumPericenterRadii_[ i ] = 1.181*1e3*(spice_interface::getBodyProperties("NEPTUNE", "RADII",3).at(0)+spice_interface::getBodyProperties("NEPTUNE", "RADII",3).at(1)+spice_interface::getBodyProperties("NEPTUNE", "RADII",3).at(2))/3;//26000.0E3;
        break;}
    case( 9 ):{
        minimumPericenterRadii_[ i ] = 9.34*1e3*(spice_interface::getBodyProperties("PLUTO", "RADII",3).at(0)+spice_interface::getBodyProperties("PLUTO", "RADII",3).at(1)+spice_interface::getBodyProperties("PLUTO", "RADII",3).at(2))/3;//25000.0E3;
        break;}
    case( 0 ):
    {
        minimumPericenterRadii_[ i ] = 1.01*1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;//2000.0E3;
        break;
    }
    default:
        std::cerr<<"Planet in flyby sequence is not defined.";
    }

}
return minimumPericenterRadii_;
}

Eigen::VectorXd createMaxPerRad (std::vector< int > flybySequence){
    size_t numberOfLegs_{flybySequence.size()};
    Eigen::VectorXd maximumPericenterRadii_;
    maximumPericenterRadii_.resize( numberOfLegs_ );
    //double factor{1.04 * 1e3};

for(int i = 0; i < numberOfLegs_; i++)
{
    switch(flybySequence[ i ])
    {
    case( 1 ):{ //1.082
        maximumPericenterRadii_[ i ] = 100 * 1e3 *(spice_interface::getBodyProperties("MERCURY", "RADII",3).at(0)+spice_interface::getBodyProperties("MERCURY", "RADII",3).at(1)+spice_interface::getBodyProperties("MERCURY", "RADII",3).at(2))/3;//2639.7E3;
        break;}
    case( 2 ):{//1.047
        maximumPericenterRadii_[ i ] = 100 * 1e3 *(spice_interface::getBodyProperties("VENUS", "RADII",3).at(0)+spice_interface::getBodyProperties("VENUS", "RADII",3).at(1)+spice_interface::getBodyProperties("VENUS", "RADII",3).at(2))/3;//6251.8E3;
        break;}
    case( 3 ): { //1.048
        maximumPericenterRadii_[ i ] = 100 * 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;//6578.1E3;
        break;}
    case( 4 ):{
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;//1395.0E3;
        break;}
    case( 5 ):{
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("JUPITER", "RADII",3).at(0)+spice_interface::getBodyProperties("JUPITER", "RADII",3).at(1)+spice_interface::getBodyProperties("JUPITER", "RADII",3).at(2))/3;//3596.2E3;
        break;}
    case( 6 ):{
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("SATURN", "RADII",3).at(0)+spice_interface::getBodyProperties("SATURN", "RADII",3).at(1)+spice_interface::getBodyProperties("SATURN", "RADII",3).at(2))/3;//72000E3;
        break;}
    case( 7 ):{
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("URANUS", "RADII",3).at(0)+spice_interface::getBodyProperties("URANUS", "RADII",3).at(1)+spice_interface::getBodyProperties("URANUS", "RADII",3).at(2))/3;//61000.0E3;
        break;}
    case( 8 ):{
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("NEPTUNE", "RADII",3).at(0)+spice_interface::getBodyProperties("NEPTUNE", "RADII",3).at(1)+spice_interface::getBodyProperties("NEPTUNE", "RADII",3).at(2))/3;//26000.0E3;
        break;}
    case( 9 ):{
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("PLUTO", "RADII",3).at(0)+spice_interface::getBodyProperties("PLUTO", "RADII",3).at(1)+spice_interface::getBodyProperties("PLUTO", "RADII",3).at(2))/3;//25000.0E3;
        break;}
    case( 0 ):
    {
        maximumPericenterRadii_[ i ] = 100*1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;//2000.0E3;
        break;
    }
    default:
        std::cerr<<"Planet in flyby sequence is not defined.";
    }

}
return maximumPericenterRadii_;
}

std::vector< std::vector< double > > createBounds (std::vector< int > flybySequence, double start_date, double launch_window, int num_DSMs, std::vector< LegTypeTransfer > LegTypeVec){
    size_t numLegs = flybySequence.size()-1;
    size_t numberOfParameters = flybySequence.size()+(num_DSMs*4);
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );

    bounds[ 0 ][ 0 ] = start_date; //MJD2000
    bounds[ 1 ][ 0 ] = start_date + launch_window; //MJD2000
   // std::cout << "Bounds: " << 0 << std::endl << bounds[0][0] <<"  " << bounds[1][0] << std::endl;

    // create traditional ToF Bounds
    enum LegOptions{ EE, Em, EV, EM, mm, mE, mV, mM, VV, VE, Vm, VM, MM, EJ, JS, SN, Vx, ME, Mm};
    int casepointer {};

    std::stringstream flybySequenceSS;
    std::copy(flybySequence.begin(), flybySequence.end(), std::ostream_iterator<int>(flybySequenceSS));
    std::string flybySequenceString = flybySequenceSS.str();

    for (int i {0}; i < numLegs; i++){
        int countLeg{i+1};

       // std::cout << "CountLeg = " << countLeg << std::endl;

        std::string legString {flybySequenceString.substr(i,2)};
        if( legString == "33" ) casepointer = EE;
        if( legString == "30" ) casepointer = Em;
        if( legString == "34" ) casepointer = EM;

        if( legString == "03" ) casepointer = mE;
        if( legString == "00" ) casepointer = mm;
        if( legString == "04" ) casepointer = mM;

        if( legString == "43" ) casepointer = ME;
        if( legString == "40" ) casepointer = Mm;
        if( legString == "44" ) casepointer = MM;





        if( legString == "32" ) casepointer = EV;
        if( legString == "02" ) casepointer = mV;
        if( legString == "22" ) casepointer = VV;
        if( legString == "23" ) casepointer = VE;
        if( legString == "20" ) casepointer = Vm;
        if( legString == "24" ) casepointer = VM;
        if( legString == "35" ) casepointer = EJ;
        if( legString == "56" ) casepointer = JS;
        if( legString == "68" ) casepointer = SN;
        if( legString == "21" ) casepointer = Vx;



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


        case (mV):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 500;
            break;

        case (VV):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 400;
            break;
        case (VE):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 250;
            break;
        case (Vm):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 250;
            break;
        case (VM):
            bounds[ 0 ][ countLeg] = 100;
            bounds[ 1 ][ countLeg ] = 1000;
            break;
        case (EJ):
            bounds[ 0 ][ countLeg ] = 400;
            bounds[ 1 ][ countLeg ] = 2000;
            break;
        case (JS):
            bounds[ 0 ][ countLeg ] = 1000;
            bounds[ 1 ][ countLeg ] = 6000;
            break;
        case (SN):
            bounds[ 0 ][ countLeg ] = 1000;
            bounds[ 1 ][ countLeg ] = 5000;
            break;
        case (Vx):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 400;
            break;
        case (EV):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 300;
            break;
        default: std::cout << "Leg input could not be recognized\n";
        }

        //std::cout << "Bounds: " << countLeg << std::endl << bounds[0][countLeg] <<std::endl << bounds[1][i+1] << std::endl;
            }
    // create additional DSM Bounds
    size_t counter{numLegs+1};
    std::vector <int> DSMpos;
     for(int i{0}; i < LegTypeVec.size(); i++){
         if(LegTypeVec.at(i)==6) {DSMpos.push_back(i);}
          }

    if(LegTypeVec[0]==5){
        bounds[ 0 ][ counter ] = 0.01;                                                                       //1st additional var - fraction of DSM
        bounds[ 1 ][ counter ] = 0.99;
        bounds[ 0 ][ counter+1 ] = 0;                                                                     //2nd additional var - V_inf at departure
        bounds[ 1 ][ counter+1 ] = 5000;
        bounds[ 0 ][ counter+2 ] = 0 * 2 * mathematical_constants::getPi< double >( );                       //3rd additional var - in-plane angle for the hyperbolic excess velocity
        bounds[ 1 ][ counter+2 ] = 1 * 2 * mathematical_constants::getPi< double >( );
        bounds[ 0 ][ counter+3 ] = std::acos(  2 * 1 - 1. ) - mathematical_constants::getPi< double >( )/2; //4th additional var - out-of-plane angle for the hyperbolic excess velocity
        bounds[ 1 ][ counter+3 ] = std::acos(  2 * 0 - 1. ) - mathematical_constants::getPi< double >( )/2;
        counter+=4;

        for (int i{1}; i<num_DSMs;i++){
            int pos = DSMpos.at(i-1);
            std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(flybySequence.at(pos));
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
    if(LegTypeVec[0]!=5){

    for (int i{0}; i<num_DSMs;i++){
        int pos = DSMpos.at(i-1);
        std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(flybySequence.at(pos));
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

    return bounds;

}
std::vector< std::vector< double > > createCassiniBounds (std::vector< int > flybySequence, double start_date, double launch_window, int num_DSMs, std::vector< std::vector< double > > bounds, std::vector< LegTypeTransfer > LegTypeVec){
    size_t numLegs = flybySequence.size()-1;
    bounds[ 0 ][ 0 ] = start_date; //MJD2000
    bounds[ 1 ][ 0 ] = start_date + launch_window; //MJD2000
   // std::cout << "Bounds: " << 0 << std::endl << bounds[0][0] <<"  " << bounds[1][0] << std::endl;

    // create traditional ToF Bounds
    enum LegOptions{ EE, Em, EV, EM, mm, mE, mV, mM, VV, VE, Vm, VM, MM, EJ, JS, SN, Vx};
    int casepointer {};

    std::stringstream flybySequenceSS;
    std::copy(flybySequence.begin(), flybySequence.end(), std::ostream_iterator<int>(flybySequenceSS));
    std::string flybySequenceString = flybySequenceSS.str();

    for (int i {0}; i < numLegs; i++){
        int countLeg{i+1};

       // std::cout << "CountLeg = " << countLeg << std::endl;

        std::string legString {flybySequenceString.substr(i,2)};
        if( legString == "33" ) casepointer = EE;
        if( legString == "30" ) casepointer = Em;
        if( legString == "32" ) casepointer = EV;
        if( legString == "34" ) casepointer = EM;
        if( legString == "00" ) casepointer = mm;
        if( legString == "03" ) casepointer = mE;
        if( legString == "02" ) casepointer = mV;
        if( legString == "04" ) casepointer = mM;
        if( legString == "22" ) casepointer = VV;
        if( legString == "23" ) casepointer = VE;
        if( legString == "20" ) casepointer = Vm;
        if( legString == "24" ) casepointer = VM;
        if( legString == "44" ) casepointer = MM;
        if( legString == "35" ) casepointer = EJ;
        if( legString == "56" ) casepointer = JS;
        if( legString == "68" ) casepointer = SN;
        if( legString == "21" ) casepointer = Vx;



        switch(casepointer){
        case (EE):
            bounds[ 0 ][ countLeg ] = 1;
            bounds[ 1 ][ countLeg ] = 400;
            break;
        case (Em):
            bounds[ 0 ][ countLeg ] = 1;
            bounds[ 1 ][ countLeg ] = 250;
            break;
        case (EV):
            bounds[ 0 ][ countLeg ] = 100;
            bounds[ 1 ][ countLeg ] = 400;
            break;
        case (EM):
            bounds[ 0 ][ countLeg ] = 50;
            bounds[ 1 ][ countLeg ] = 360;
            break;
        case (mm):
            bounds[ 0 ][ countLeg ] = 1;
            bounds[ 1 ][ countLeg ] = 100;
            break;
        case (mE):
            bounds[ 0 ][ countLeg ] = 1;
            bounds[ 1 ][ countLeg ] = 250;
            break;
        case (mV):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 500;
            break;
        case (mM):
            bounds[ 0 ][ countLeg ] = 50;
            bounds[ 1 ][ countLeg ] = 600;
            break;
        case (VV):
            bounds[ 0 ][ countLeg ] = 100;
            bounds[ 1 ][ countLeg ] = 500;
            break;
        case (VE):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 300;
            break;
        case (Vm):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 500;
            break;
        case (VM):
            bounds[ 0 ][ countLeg] = 100;
            bounds[ 1 ][ countLeg ] = 1000;
            break;
        case (MM):
            bounds[ 0 ][ countLeg] = 1;
            bounds[ 1 ][ countLeg ] = 500;
            break;
        case (EJ):
            bounds[ 0 ][ countLeg ] = 400;
            bounds[ 1 ][ countLeg ] = 1600;
            break;
        case (JS):
            bounds[ 0 ][ countLeg ] = 800;
            bounds[ 1 ][ countLeg ] = 2200;
            break;
        case (SN):
            bounds[ 0 ][ countLeg ] = 1000;
            bounds[ 1 ][ countLeg ] = 5000;
            break;
        case (Vx):
            bounds[ 0 ][ countLeg ] = 30;
            bounds[ 1 ][ countLeg ] = 400;
            break;
        default: std::cout << "Leg input could not be recognized\n";
        }

        //std::cout << "Bounds: " << countLeg << std::endl << bounds[0][countLeg] <<std::endl << bounds[1][i+1] << std::endl;
            }
    // create additional DSM Bounds
    size_t counter{numLegs+1};
    std::vector <int> DSMpos;
     for(int i{0}; i < LegTypeVec.size(); i++){
         if(LegTypeVec.at(i)==6) {DSMpos.push_back(i);}
          }

    if(LegTypeVec[0]==5){
        bounds[ 0 ][ counter ] = 0.01;                                                                       //1st additional var - fraction of DSM
        bounds[ 1 ][ counter ] = 0.90;
        bounds[ 0 ][ counter+1 ] = 3000;                                                                     //2nd additional var - V_inf at departure
        bounds[ 1 ][ counter+1 ] = 5000;
        bounds[ 0 ][ counter+2 ] = 0 * 2 * mathematical_constants::getPi< double >( );                       //3rd additional var - in-plane angle for the hyperbolic excess velocity
        bounds[ 1 ][ counter+2 ] = 1 * 2 * mathematical_constants::getPi< double >( );
        bounds[ 0 ][ counter+3 ] = std::acos(  2 * 1 - 1. ) - mathematical_constants::getPi< double >( )/2; //4th additional var - out-of-plane angle for the hyperbolic excess velocity
        bounds[ 1 ][ counter+3 ] = std::acos(  2 * 0 - 1. ) - mathematical_constants::getPi< double >( )/2;
        counter+=4;

        for (int i{1}; i<num_DSMs;i++){
            int pos = DSMpos.at(i-1);
            std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(flybySequence.at(pos));
            Eigen::VectorXd minPerRad = createMinPerRad(SwingbyPLanet);
            Eigen::VectorXd maxPerRad = createMaxPerRad(SwingbyPLanet);

            bounds[ 0 ][ counter ] = 0.01;                                                 //1st additional var - fraction of DSM
            bounds[ 1 ][ counter ] = 0.90;
            bounds[ 0 ][ counter+1 ] = -1 * mathematical_constants::getPi< double >( ); //2nd additional var - rotation angle of the GA
            bounds[ 1 ][ counter+1 ] = 1 * mathematical_constants::getPi< double >();
            bounds[ 0 ][ counter+2 ] = minPerRad[0];                                    //3rd additional var - pericenter radius of the GA
            bounds[ 1 ][ counter+2 ] = maxPerRad[0];
            bounds[ 0 ][ counter+3 ] = 0; //4th additional var - V added for the powered GA
            bounds[ 1 ][ counter+3 ] = 2000;
            counter+=4;
        }


    }
    if(LegTypeVec[0]!=5){

    for (int i{0}; i<num_DSMs;i++){
        int pos = DSMpos.at(i-1);
        std::vector<int> SwingbyPLanet; SwingbyPLanet.push_back(flybySequence.at(pos));
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

    return bounds;

}




std::vector< LegTypeTransfer > createLegType (std::string flyby_seq){
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'G'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'c'), flyby_seq.end());

    size_t numLegs = flyby_seq.length( )-1;
    std::vector< LegTypeTransfer > legTypeVector;
    size_t count{0};
    while (count <= numLegs ){
        char next{};
        if((count)<numLegs){next=flyby_seq.at(count+1);}
        if(count==0 && next=='d') {legTypeVector.push_back(MGA1DsmVelocity_Departure); count=count+2; continue;}
        if(count==0 && next!='d'){legTypeVector.push_back(MGA_Departure);count++;continue;}
        if(count!=0 && next=='d'){legTypeVector.push_back(MGA1DsmVelocity_Swingby); count=count+2; continue;}
        if(count!=0 && count!=numLegs && next!='d'){legTypeVector.push_back(MGA_Swingby);count++;continue;}
        if(count==numLegs){legTypeVector.push_back(Capture);count++;}
    }
return legTypeVector;
}


std::vector<std::shared_ptr< ephemerides::Ephemeris >> createEphemerisVec (std::vector< int > flybySequence, std::string base  ,std::string frame, bool useBary){
    size_t numberOfLegs_{flybySequence.size()};
    std::vector<std::shared_ptr< ephemerides::Ephemeris >> ephemerisVector_;
    ephemerisVector_.resize( numberOfLegs_ );
    Eigen::VectorXd minimumPericenterRadii_;
    minimumPericenterRadii_.resize( numberOfLegs_ );
    Eigen::VectorXd gravitationalParameterVector_;
    gravitationalParameterVector_.resize( numberOfLegs_ );

for(int i = 0; i < numberOfLegs_; i++)
{
    switch(flybySequence[ i ])
    {
    case( 1 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "MERCURY_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "MERCURY" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}


    case( 2 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "VENUS_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "VENUS" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 3 ): {
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "EARTH_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "EARTH" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}


    case( 4 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "MARS_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "MARS" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 5 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "JUPITER_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "JUPITER" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 6 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "SATURN_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "SATURN" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 7 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "URANUS_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "URANUS" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 8 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "NEPTUNE_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "NEPTUNE" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 9 ):{
        if(useBary){
        std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
        std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "PLUTO_BARYCENTER" );
        ephemerisVector_[ i ] = spiceEphemeris;
        }
        else {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "PLUTO" );
            ephemerisVector_[ i ] = spiceEphemeris;
        }
        break;}

    case( 0 ):
    {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "MOON" );
            ephemerisVector_[ i ] = spiceEphemeris;

        break;}

        //Special Case
    case( 10 ):
    {
            std::shared_ptr< simulation_setup::EphemerisSettings > spiceEphemerisSettings = std::make_shared< simulation_setup::DirectSpiceEphemerisSettings >( base, frame );
            std::shared_ptr< ephemerides::Ephemeris > spiceEphemeris = createBodyEphemeris( spiceEphemerisSettings, "SUN" );
            ephemerisVector_[ i ] = spiceEphemeris;

        break;}



    default:
        std::cerr<<"Planet in flyby sequence is not defined.";
    }

}
return ephemerisVector_;
}


Eigen::VectorXd createGravParamVec (std::vector< int > flybySequence){
    size_t numberOfLegs_{flybySequence.size()};
    Eigen::VectorXd gravitationalParameterVector_;
    gravitationalParameterVector_.resize( numberOfLegs_ );

for(int i = 0; i < numberOfLegs_; i++)
{
    switch(flybySequence[ i ])
    {
    case( 1 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("MERCURY_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 2 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("VENUS_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 3 ): {
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 4 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("MARS_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 5 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("JUPITER_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 6 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("SATURN_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 7 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("URANUS_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 8 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("NEPTUNE_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 9 ):{
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("PLUTO_BARYCENTER", "GM").at(0)*1e9;
        break;}
    case( 0 ):
    {
        gravitationalParameterVector_[ i ] = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        break;
    }
    default:
        std::cerr<<"Planet in flyby sequence is not defined.";
    }

}
return gravitationalParameterVector_;
}


std::vector< std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > > getShapingBasisFunctions(
        const double timeOfFlight, const int numberOfRevolutions )
{
    Eigen::VectorXd dummyVector;

    // Get recommended base functions for the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    getRecommendedRadialVelocityBaseFunctions(
                radialVelocityFunctionComponents, dummyVector, timeOfFlight );

    // Get recommended base functions for the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    getRecommendedNormalAxialBaseFunctions(
                normalVelocityFunctionComponents, dummyVector, timeOfFlight );

    // Get recommended base functions for the axial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    getRecommendedAxialVelocityBaseFunctions(
                axialVelocityFunctionComponents, dummyVector, timeOfFlight, numberOfRevolutions );

    {
        double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
        double scaleFactor = 1.0 / timeOfFlight;

        std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
                std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );
        std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
                std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                    1.0, 0.5 * frequency, scaleFactor );

        // Add two additional base functions
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
        radialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

    }

    {
        double scaleFactor = 1.0 / timeOfFlight;

        // Create base function settings for the components of the axial velocity composite function.
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 3.0, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 4.0, scaleFactor );
        std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
                std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 5.0, scaleFactor );


        // Set components for the axial velocity function.
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, firstAxialVelocityBaseFunctionSettings ) );
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondAxialVelocityBaseFunctionSettings ) );
        axialVelocityFunctionComponents.push_back(
                    createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdAxialVelocityBaseFunctionSettings ) );

    }

    return { radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents };
}

Eigen::Vector6d getDepState( const std::function< Eigen::Vector6d( const double )> StateFunction,
                                   double Time, std::string trajectory, bool useDepOrbit, double sma ){
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char dep = trajectory.at(0);
    int switcharg{};
    if(dep=='E')switcharg=13;   if(dep=='m')switcharg=10;   if(dep=='M')switcharg=14;  //Earth & Moon & Mars
    if(dep=='1')switcharg=1;    if(dep=='2')switcharg=2;    if(dep=='3')switcharg=3;   //Lagrange points
    if(dep=='4')switcharg=4;    if(dep=='5')switcharg=5;    if(dep=='6')switcharg=6;   //Lagrange points
    if(dep!='E'&& dep!='m' && dep!='M' && dep!='1'&&dep!='2'&&dep!='3'&&dep!='4'&&dep!='5'&&dep!='6')std::cout << "\nNo switcharg was generated\n";

    Eigen::Vector6d State = StateFunction( Time );

    double R_earth = 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;
    double R_moon = 1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;
    double R_mars = 1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;

    double h_earth = 250e3;
    double h_moon = 100e3;
    double h_mars = 250e3;
    double V_earth = sqrt(celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/(R_earth+h_earth));
    double V_moon = sqrt(celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER/(R_moon+h_moon));
    double V_mars = sqrt(celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/(R_mars+h_mars));
    double P_earth = mathematical_constants::PI*2*(R_earth+h_earth)/V_earth;
    double P_moon = mathematical_constants::PI*2*(R_moon+h_moon)/V_moon;
    double P_mars = mathematical_constants::PI*2*(R_mars+h_mars)/V_mars;


    double Rsoi_earth = 1.495978e11 * pow(celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER, 0.4);
    double Rsoi_mars = 2.27940e11 * pow(celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER, 0.4);


    switch (switcharg) {
        case(13):{ //Earth - parking orbit at 250 km
        if(useDepOrbit){
            double Vinf = sqrt(-1* celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/sma);
            Eigen::Vector6d sphericalState_E;
            sphericalState_E << Rsoi_earth, mathematical_constants::PI*2/P_earth*Time  , mathematical_constants::PI/2, Vinf, 0, 0;
            Eigen::Vector6d cartesianState_E = convertSphericalToCartesianState(sphericalState_E);
            State = State + cartesianState_E;}
            break;}

        case(10):{ //Moon - parking orbit at 100 km
        if(useDepOrbit){
        Eigen::Vector6d sphericalState_m;
        sphericalState_m << (h_moon+R_moon), mathematical_constants::PI*2/P_moon*Time  , mathematical_constants::PI/2, 0, mathematical_constants::PI*2/P_moon, 0;
        Eigen::Vector6d cartesianState_m = convertSphericalToCartesianState(sphericalState_m);
        State = State + cartesianState_m;}
            break;}

        case(14):{ //Mars - parking orbit at 250 km
        if(useDepOrbit){
        double Vinf = sqrt(-celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/sma);
        Eigen::Vector6d sphericalState_M;
        sphericalState_M << Rsoi_mars, mathematical_constants::PI*2/P_mars*Time  , mathematical_constants::PI/2, Vinf, 0, 0;
        Eigen::Vector6d cartesianState_M = convertSphericalToCartesianState(sphericalState_M);
        State = State + cartesianState_M;}
            break;}

        case(1):{ //EmL1 - stationary
            State = getLibrationState ("Em", 1, Time, true );
            break;}

        case(2):{ //EmL2 - stationary
            State = getLibrationState ("Em", 2, Time, true );
            break;}

        case(3):{ //EmL4 - stationary
            State = getLibrationState ("Em", 4, Time, true );
            break;}

        case(4):{ //EmL5 - stationary
            State = getLibrationState ("Em", 5, Time, true );
            break;}

        case(5):{ //SEL2 - stationary
            State = getLibrationState ("SE", 2, Time );
            break;}

        case(6):{ //SML1 - stationary
            State = getLibrationState ("SM", 1, Time );
            break;}

    default: std::cout << "\nSwitcharg not recognized\n";


    }
    return State;
}

Eigen::Vector6d getDepartureState (const std::function< Eigen::Vector6d( const double )> StateFunction, double Time, std::string trajectory, Eigen::Vector3d Vinf_vec){

    Eigen::Vector6d boundaryState = StateFunction(Time);

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char bound = trajectory.at(0);
    int switcharg{};
    if(bound=='E')switcharg=13;   if(bound=='m')switcharg=10;   if(bound=='M')switcharg=14;  //Earth & Moon & Mars
    if(bound=='1')switcharg=1;    if(bound=='2')switcharg=2;    if(bound=='3')switcharg=3;   //Lagrange points
    if(bound=='4')switcharg=4;    if(bound=='5')switcharg=5;    if(bound=='6')switcharg=6;   //Lagrange points
    if(bound!='E'&& bound!='m' && bound!='M' && bound!='1'&&bound!='2'&&bound!='3'&&bound!='4'&&bound!='5'&&bound!='6')std::cout << "\nNo switcharg was generated\n";

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
            boundaryState = getLibrationState ("Em", 1, Time, true );
            break;}

        case(2):{ //EmL2 - stationary
            boundaryState = getLibrationState ("Em", 2, Time, true );
            break;}

        case(3):{ //EmL4 - stationary
            boundaryState = getLibrationState ("Em", 4, Time, true );
            break;}

        case(4):{ //EmL5 - stationary
            boundaryState = getLibrationState ("Em", 5, Time, true );
            break;}

        case(5):{ //SEL2 - stationary
            boundaryState = getLibrationState ("SE", 2, Time );
            break;}

        case(6):{ //SML1 - stationary
            boundaryState = getLibrationState ("SM", 1, Time );
            break;}

    default: std::cout << "\nSwitcharg not recognized\n";


    }
    return boundaryState;
}

Eigen::Vector6d getArrivalState (const std::function< Eigen::Vector6d( const double )> StateFunction, double Time, std::string trajectory, Eigen::Vector3d Vinf_vec){

    Eigen::Vector6d boundaryState = StateFunction(Time);

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char bound = trajectory.back();
    int switcharg{};
    if(bound=='E')switcharg=13;   if(bound=='m')switcharg=10;   if(bound=='M')switcharg=14;  //Earth & Moon & Mars
    if(bound=='1')switcharg=1;    if(bound=='2')switcharg=2;    if(bound=='3')switcharg=3;   //Lagrange points
    if(bound=='4')switcharg=4;    if(bound=='5')switcharg=5;    if(bound=='6')switcharg=6;   //Lagrange points
    if(bound!='E'&& bound!='m' && bound!='M' && bound!='1'&&bound!='2'&&bound!='3'&&bound!='4'&&bound!='5'&&bound!='6')std::cout << "\nNo switcharg was generated\n";

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
            boundaryState = getLibrationState ("Em", 1, Time, true );
            break;}

        case(2):{ //EmL2 - stationary
            boundaryState = getLibrationState ("Em", 2, Time, true );
            break;}

        case(3):{ //EmL4 - stationary
            boundaryState = getLibrationState ("Em", 4, Time, true );
            break;}

        case(4):{ //EmL5 - stationary
            boundaryState = getLibrationState ("Em", 5, Time, true );
            break;}

        case(5):{ //SEL2 - stationary
            boundaryState = getLibrationState ("SE", 2, Time );
            break;}

        case(6):{ //SML1 - stationary
            boundaryState = getLibrationState ("SM", 1, Time );
            break;}

    default: std::cout << "\nSwitcharg not recognized\n";


    }
    return boundaryState;
}






Eigen::Vector6d getArrState( const std::function< Eigen::Vector6d( const double )> StateFunction,
                                   double Time, std::string trajectory, bool useArrOrbit, double sma ){

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char arr = trajectory.back();
    int switcharg{};
    if(arr=='E')switcharg=13;   if(arr=='m')switcharg=10;   if(arr=='M')switcharg=14;  //Earth & Moon & Mars
    if(arr=='1')switcharg=1;    if(arr=='2')switcharg=2;    if(arr=='3')switcharg=3;   //Lagrange points
    if(arr=='4')switcharg=4;    if(arr=='5')switcharg=5;    if(arr=='6')switcharg=6;   //Lagrange points

    if(arr!='E'&&arr!='m'&&arr!='M'&&arr!='1'&&arr!='2'&&arr!='3'&&arr!='4'&&arr!='5'&&arr!='6')std::cout << "\nNo switcharg was generated\n";

    Eigen::Vector6d State = StateFunction( Time );

    double R_earth = 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;
    double R_moon = 1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;
    double R_mars = 1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;

    double h_earth = 250e3;
    double h_moon = 100e3;
    double h_mars = 250e3;
    double V_earth = sqrt(celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/(R_earth+h_earth));
    double V_moon = sqrt(celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER/(R_moon+h_moon));
    double V_mars = sqrt(celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/(R_mars+h_mars));
    double P_earth = mathematical_constants::PI*2*(R_earth+h_earth)/V_earth;
    double P_moon = mathematical_constants::PI*2*(R_moon+h_moon)/V_moon;
    double P_mars = mathematical_constants::PI*2*(R_mars+h_mars)/V_mars;

    double Rsoi_earth = 1.495978e11 * pow(celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER, 0.4);
    double Rsoi_mars = 2.27940e11 * pow(celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER, 0.4);


    switch (switcharg) {
        case(13):{ //Earth - parking orbit at 250 km
        if(useArrOrbit){
            double Vinf = sqrt(-celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/sma);
            Eigen::Vector6d sphericalState_E;
            sphericalState_E << Rsoi_earth, mathematical_constants::PI*2/P_earth*Time  , mathematical_constants::PI/2, -Vinf, 0, 0;
            Eigen::Vector6d cartesianState_E = convertSphericalToCartesianState(sphericalState_E);
            State = State + cartesianState_E;}

            break;}


        case(10):{ //Moon - parking orbit at 100 km
        if(useArrOrbit){
        Eigen::Vector6d sphericalState_m;
        sphericalState_m << (h_moon+R_moon), mathematical_constants::PI*2/P_moon*Time  , mathematical_constants::PI/2, 0, mathematical_constants::PI*2/P_moon, 0;
        Eigen::Vector6d cartesianState_m = convertSphericalToCartesianState(sphericalState_m);
        State = State + cartesianState_m;}
            break;}

        case(14):{ //Mars - parking orbit at 250 km
        if(useArrOrbit){
            double Vinf = sqrt(-celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/sma);
            Eigen::Vector6d sphericalState_M;
            sphericalState_M << Rsoi_mars, mathematical_constants::PI*2/P_mars*Time  , mathematical_constants::PI/2, -Vinf, 0, 0;
            Eigen::Vector6d cartesianState_M = convertSphericalToCartesianState(sphericalState_M);
            State = State + cartesianState_M;}
            break;}

        case(1):{ //EmL1 - stationary
            State = getLibrationState ("Em", 1, Time, true );
            break;}

        case(2):{ //EmL2 - stationary
            State = getLibrationState ("Em", 2, Time, true );
            break;}

        case(3):{ //EmL4 - stationary
            State = getLibrationState ("Em", 4, Time, true );
            break;}

        case(4):{ //EmL5 - stationary
            State = getLibrationState ("Em", 5, Time, true );
            break;}

        case(5):{ //SEL2 - stationary
            State = getLibrationState ("SE", 2, Time );
            break;}

        case(6):{ //SML1 - stationary
            State = getLibrationState ("SM", 1, Time );
            break;}

    default: std::cout << "\nSwitcharg not recognized\n";


    }
    return State;
}



Eigen::Vector6d getLibrationState (std::string system, int point, double time, bool SEm){
    Eigen::Vector6d LibrationState;

    // parameter definitions
        //Gravitational parameters
         const double sunGravitationalParameter = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
         const double earthGravitationalParameter = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
         const double EBGravitationalParameter = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
         const double marsBGravitationalParameter = spice_interface::getBodyProperties("MARS_BARYCENTER", "GM").at(0)*1e9;
         const double moonGravitationalParameter = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;

        //LibrsationPoints
        circular_restricted_three_body_problem::LagrangePoint SE ( sunGravitationalParameter, EBGravitationalParameter,
                                                                      std::make_shared< root_finders::NewtonRaphson >( 1.0e-14, 1000 ) );
        circular_restricted_three_body_problem::LagrangePoint SM ( sunGravitationalParameter, marsBGravitationalParameter,
                                                                      std::make_shared< root_finders::NewtonRaphson >( 1.0e-14, 1000 ) );
        circular_restricted_three_body_problem::LagrangePoint Em ( earthGravitationalParameter, moonGravitationalParameter,
                                                                      std::make_shared< root_finders::NewtonRaphson >( 1.0e-14, 1000 ) );

        //Ephimerides
        std::vector<int> SunEB = {3}; Eigen::Vector6d distanceSunEB = createEphemerisVec(SunEB, "SUN", "ECLIPJ2000", true)[0]->getCartesianState( time );
        std::vector<int> SunMB = {4}; Eigen::Vector6d distanceSunMB = createEphemerisVec(SunMB, "SUN", "ECLIPJ2000", true)[0]->getCartesianState( time );
        std::vector<int> SunE = {3}; Eigen::Vector6d distanceSunEarth = createEphemerisVec(SunE, "SUN", "ECLIPJ2000", false)[0]->getCartesianState( time );
        std::vector<int> vecEm = {0}; Eigen::Vector6d distanceEarthMoon = createEphemerisVec(vecEm, "EARTH", "ECLIPJ2000", false)[0]->getCartesianState( time );


        //Returning the desired state vector
        if(system == "Em"||system == "EM"){
           LibrationState = Em.getLagrangeStateVec(Em, point, distanceEarthMoon);
        }
        if(system == "SE"||system == "se"){
           LibrationState = SE.getLagrangeStateVec(SE, point, distanceSunEB);
        }
        if(system == "SM"||system == "Sm"){
           LibrationState = SM.getLagrangeStateVec(SM, point, distanceSunMB);
        }

        if(SEm){
            LibrationState+=distanceSunEarth;
        }


        return LibrationState;

}


Eigen::Vector6d convertSphericalToCartesianState( const Eigen::Vector6d& sphericalState ){

    // Create Cartesian state vector, initialized with zero entries.
    Eigen::Vector6d convertedCartesianState = Eigen::Vector6d::Zero( );

    // Create local variables.
    const double radius = sphericalState( 0 );
    const double azimuthAngle = sphericalState( 1 );
    const double elevationAngle = sphericalState( 2 );
    const double radius_dot = sphericalState( 3 );
    const double azimuthAngle_dot = sphericalState( 4 );
    const double elevationAngle_dot = sphericalState( 5 );

    // Precompute sine/cosine of angles, which has multiple usages, to save computation time.
    const double cosineOfElevationAngle = std::cos( elevationAngle );
    const double sineOfElevationAngle = std::sin( elevationAngle );
    const double cosineOfAzimuthAngle = std::cos( azimuthAngle );
    const double sineOfAzimuthAngle = std::sin( azimuthAngle );

    // Perform transformation of position coordinates.
    convertedCartesianState( 0 ) = radius * cosineOfAzimuthAngle * sineOfElevationAngle;
    convertedCartesianState( 1 ) = radius * sineOfAzimuthAngle * sineOfElevationAngle;
    convertedCartesianState( 2 ) = radius * cosineOfElevationAngle;

    convertedCartesianState( 3 ) = radius_dot * cosineOfAzimuthAngle*sineOfElevationAngle - radius * sineOfAzimuthAngle * sineOfElevationAngle*azimuthAngle_dot + radius * cosineOfAzimuthAngle * cosineOfElevationAngle * elevationAngle_dot;
    convertedCartesianState( 4 ) = radius_dot * sineOfAzimuthAngle*sineOfElevationAngle + radius * cosineOfAzimuthAngle * sineOfElevationAngle * azimuthAngle_dot + radius * sineOfAzimuthAngle * cosineOfElevationAngle * elevationAngle_dot;
    convertedCartesianState( 5 ) = radius_dot * cosineOfElevationAngle - radius * sineOfElevationAngle * elevationAngle_dot;

    return convertedCartesianState;

}

std::vector< std::vector< double > > createContBounds ( std::pair< double, double > departureTimeBounds,  std::pair< double, double > timeOfFlightBounds, bool useDepOrbit, bool useArrOrbit){

    double maxVinf = 5e3;

// No orbits
 if(!useDepOrbit && !useArrOrbit){

std::vector< std::vector< double > > bounds( 2, std::vector< double >( 7, 0.0 ) );
bounds[ 0 ][ 0 ] =  departureTimeBounds.first; //1a
bounds[ 1 ][ 0 ] = departureTimeBounds.second; //1b
bounds[ 0 ][ 1 ] = timeOfFlightBounds.first;   //2a
bounds[ 1 ][ 1 ] = timeOfFlightBounds.second;  //2b
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
return bounds;

 }

if(useDepOrbit && !useArrOrbit){
// Only Dep orbit
std::vector< std::vector< double > > bounds( 2, std::vector< double >( 10, 0.0 ) );
bounds[ 0 ][ 0 ] =  departureTimeBounds.first; //1a
bounds[ 1 ][ 0 ] = departureTimeBounds.second; //1b
bounds[ 0 ][ 1 ] = timeOfFlightBounds.first;   //2a
bounds[ 1 ][ 1 ] = timeOfFlightBounds.second;  //2b
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
bounds[ 0 ][ 7 ] = -1 * maxVinf;               //8a Dep Vinf_x
bounds[ 1 ][ 7 ] = 1 * maxVinf;                //8b
bounds[ 0 ][ 8 ] = -1 * maxVinf;               //9a Dep Vinf_y
bounds[ 1 ][ 8 ] = 1 * maxVinf;                //9b
bounds[ 0 ][ 9 ] = -1 * maxVinf;               //10a Dep Vinf_z
bounds[ 1 ][ 9 ] = 1 * maxVinf;                //10b
return bounds;

}

if(!useDepOrbit && useArrOrbit){
// Only Arr orbit
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 10, 0.0 ) );
    bounds[ 0 ][ 0 ] =  departureTimeBounds.first; //1a
    bounds[ 1 ][ 0 ] = departureTimeBounds.second; //1b
    bounds[ 0 ][ 1 ] = timeOfFlightBounds.first;   //2a
    bounds[ 1 ][ 1 ] = timeOfFlightBounds.second;  //2b
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
    bounds[ 0 ][ 7 ] = -1 * maxVinf;               //8a Arr Vinf_x
    bounds[ 1 ][ 7 ] = 1 * maxVinf;                //8b
    bounds[ 0 ][ 8 ] = -1 * maxVinf;               //9a Arr Vinf_y
    bounds[ 1 ][ 8 ] = 1 * maxVinf;                //9b
    bounds[ 0 ][ 9 ] = -1 * maxVinf;               //10a Arr Vinf_z
    bounds[ 1 ][ 9 ] = 1 * maxVinf;                //10b
return bounds;

}

if(useDepOrbit && useArrOrbit){
// Only Dep orbit
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( 13, 0.0 ) );
    bounds[ 0 ][ 0 ] =  departureTimeBounds.first; //1a
    bounds[ 1 ][ 0 ] = departureTimeBounds.second; //1b
    bounds[ 0 ][ 1 ] = timeOfFlightBounds.first;   //2a
    bounds[ 1 ][ 1 ] = timeOfFlightBounds.second;  //2b
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
    bounds[ 0 ][ 7 ] = -1 * maxVinf;               //8a Dep Vinf_x
    bounds[ 1 ][ 7 ] = 1 * maxVinf;                //8b
    bounds[ 0 ][ 8 ] = -1 * maxVinf;               //9a Dep Vinf_y
    bounds[ 1 ][ 8 ] = 1 * maxVinf;                //9b
    bounds[ 0 ][ 9 ] = -1 * maxVinf;               //10a Dep Vinf_z
    bounds[ 1 ][ 9 ] = 1 * maxVinf;                //10b
    bounds[ 0 ][ 10 ] = -1 * maxVinf;               //11a Arr Vinf_x
    bounds[ 1 ][ 10 ] = 1 * maxVinf;                //11b
    bounds[ 0 ][ 11 ] = -1 * maxVinf;               //12a Arr Vinf_y
    bounds[ 1 ][ 11 ] = 1 * maxVinf;                //12b
    bounds[ 0 ][ 12 ] = -1 * maxVinf;               //13a Arr Vinf_z
    bounds[ 1 ][ 12 ] = 1 * maxVinf;                //13b
return bounds;

}
}

double getDepDeltaV( Eigen::Vector3d Vinf_vec_dep, std::string trajectory ){
    double departureDeltaV;

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char dep = trajectory.at(0);

    double R_earth = 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;
    double R_moon = 1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;
    double R_mars = 1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;

    double h_earth = 250e3;
    double h_moon = 100e3;
    double h_mars = 250e3;
    double V_earth = sqrt(celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/(R_earth+h_earth));
    double V_moon = sqrt(celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER/(R_moon+h_moon));
    double V_mars = sqrt(celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/(R_mars+h_mars));

    double Vcnul{};
    double Vesc{};
    double Vinf{};

    if(dep=='E'){Vcnul=V_earth;Vesc=sqrt(2*celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/(R_earth+h_earth));Vinf= Vinf_vec_dep.norm();}
    if(dep=='M'){Vcnul=V_mars;Vesc=sqrt(2*celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/(R_mars+h_mars));Vinf= Vinf_vec_dep.norm();}
    if(dep=='m'){Vcnul=V_moon;Vesc=sqrt(2*celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER/(R_moon+h_moon));Vinf= Vinf_vec_dep.norm();}



    departureDeltaV = sqrt(pow(Vesc,2)+pow(Vinf,2)) - Vcnul;
    //std::cout << "\n--------------------------\n Vesc = " << Vesc << "\n Vinf = " << Vinf << "\n Vcnul = " << Vcnul << "\n----\n depDeltaV = " << departureDeltaV << "\n--------------------------\n";
    return departureDeltaV;
}

double getArrDeltaV( Eigen::Vector3d Vinf_vec_arr, std::string trajectory ){
    double arrivalDeltaV;

    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'd'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'G'), trajectory.end());
    trajectory.erase(std::remove(trajectory.begin(), trajectory.end(), 'c'), trajectory.end());

    char arr = trajectory.back();

    double R_earth = 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;
    double R_moon = 1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;
    double R_mars = 1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;
    double h_earth = 250e3;
    double h_moon = 100e3;
    double h_mars = 250e3;
    double V_earth = sqrt(celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/(R_earth+h_earth));
    double V_moon = sqrt(celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER/(R_moon+h_moon));
    double V_mars = sqrt(celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/(R_mars+h_mars));

    double Vcnul{};
    double Vesc{};
    double Vinf{};

    if(arr=='E'){Vcnul=V_earth;Vesc=sqrt(2*celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER/(R_earth+h_earth)); Vinf= Vinf_vec_arr.norm();}
    if(arr=='M'){Vcnul=V_mars;Vesc=sqrt(2*celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER/(R_mars+h_mars)); Vinf= Vinf_vec_arr.norm();}
    if(arr=='m'){Vcnul=V_moon;Vesc=sqrt(2*celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER/(R_moon+h_moon)); Vinf= Vinf_vec_arr.norm();}

    arrivalDeltaV = sqrt(pow(Vesc,2)+pow(Vinf,2)) - Vcnul;

   // std::cout << "\n---------------------------------\n Vesc : " << Vesc << "             -            Vinf : " << Vinf << "             -            Vcnul : " << Vcnul<< "\n-------------------------\n";
   // std::cout << "----- Arrival DeltaV -----" << arrivalDeltaV << std::endl;

    return arrivalDeltaV;
}

Eigen::Vector2d getEccentricities ( std::string flyby_seq ){
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'd'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'G'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'c'), flyby_seq.end());

    Eigen::Vector2d Eccentricities(2);
    double EccDep{}; double EccCap{};

    //All planetary parking orbits are circular orbits, so Ecc is always 0
    // All planetary gateway orbits are assumed to be elliptical with e = 0.5
    //We now assume that gateways are stationary in Lagrange point, this can be modelled by a circular orbit with infinite semi-major axis
    char dep{};
    for(int i{0}; i<2; i++){
        if(i==0)  dep = flyby_seq.at(0);
        if(i==1)  dep = flyby_seq.back();

        int switcharg;
        if(dep=='E')switcharg=13;   if(dep=='m')switcharg=10;   if(dep=='M')switcharg=14;  //Earth & Moon & Mars
        if(dep=='1')switcharg=1;    if(dep=='2')switcharg=2;    if(dep=='3')switcharg=3;   //Lagrange points - GW
        if(dep=='4')switcharg=4;    if(dep=='5')switcharg=5;    if(dep=='6')switcharg=6;   //Lagrange points - GW
        if(dep=='7')switcharg=7;    if(dep=='8')switcharg=8;    if(dep=='9')switcharg=9;   //Celestial body orbit - GW

        if(dep!='E'&&dep!='m'&&dep!='M'&&dep!='1'&&dep!='2'&&dep!='3'&&dep!='4'&&dep!='5'&&dep!='6'&&dep!='7'&&dep!='8'&&dep!='9')std::cout << "\nNo switcharg was generated\n";

        switch (switcharg) {
            case 13:
            case 10:
            case 14:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:{
                if(i==0)  EccDep = 0;
                if(i==1)  EccCap = 0;
                break;
                }
          }
    }

    Eccentricities << EccDep, EccCap;
    return Eccentricities;
}

Eigen::Vector2d getSMA ( std::string flyby_seq ){
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'd'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'G'), flyby_seq.end());
    flyby_seq.erase(std::remove(flyby_seq.begin(), flyby_seq.end(), 'c'), flyby_seq.end());

    Eigen::Vector2d SMAs(2);
    double SMADep{}; double SMACap{};
    double R_earth = 1e3*(spice_interface::getBodyProperties("EARTH", "RADII",3).at(0)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(1)+spice_interface::getBodyProperties("EARTH", "RADII",3).at(2))/3;
    double R_mars = 1e3*(spice_interface::getBodyProperties("MARS", "RADII",3).at(0)+spice_interface::getBodyProperties("MARS", "RADII",3).at(1)+spice_interface::getBodyProperties("MARS", "RADII",3).at(2))/3;
    double R_moon = 1e3*(spice_interface::getBodyProperties("MOON", "RADII",3).at(0)+spice_interface::getBodyProperties("MOON", "RADII",3).at(1)+spice_interface::getBodyProperties("MOON", "RADII",3).at(2))/3;

    double R_p_earth = 250e3 + R_earth;
    double R_p_mars = 250e3 + R_mars;
    double R_p_moon = 100e3 + R_moon;

    double R_a_earth = ((1+0.5)/(1-0.5))*R_p_earth;
    double R_a_mars = ((1+0.5)/(1-0.5))*R_p_mars;
    double R_a_moon = ((1+0.5)/(1-0.5))*R_p_moon;

    //All planetary parking orbits are circular orbits, so SMA is radius of orbit
    // All planetary gateway orbits are assumed to be elliptical with e = 0.5 and Rp = parking orbit radius
    //We now assume that gateways are stationary in Lagrange point, this can be modelled by a circular orbit with infinite semi-major axis
    char dep{};
    for(int i{0}; i<2; i++){
        if(i==0)  dep = flyby_seq.at(0);
        if(i==1)  dep = flyby_seq.back();

    int switcharg;
    if(dep=='E')switcharg=13;   if(dep=='m')switcharg=10;   if(dep=='M')switcharg=14;  //Earth & Moon & Mars
    if(dep=='1')switcharg=1;    if(dep=='2')switcharg=2;    if(dep=='3')switcharg=3;   //Lagrange points - GW
    if(dep=='4')switcharg=4;    if(dep=='5')switcharg=5;    if(dep=='6')switcharg=6;   //Lagrange points - GW
    if(dep=='7')switcharg=7;    if(dep=='8')switcharg=8;    if(dep=='9')switcharg=9;   //Celestial body orbit - GW

    if(dep!='E'&&dep!='m'&&dep!='M'&&dep!='1'&&dep!='2'&&dep!='3'&&dep!='4'&&dep!='5'&&dep!='6'&&dep!='7'&&dep!='8'&&dep!='9')std::cout << "\nNo switcharg was generated\n";

    switch (switcharg) {
        case 13:{
            if(i==0)  SMADep = R_p_earth;
            if(i==1)  SMACap = R_p_earth;
            break;
        }
        case 10:{
            if(i==0)  SMADep = R_p_moon;
            if(i==1)  SMACap = R_p_moon;
            break;
        }
        case 14: {
            if(i==0)  SMADep = R_p_mars;
            if(i==1)  SMACap = R_p_mars;
            break;
        }
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:{
            if(i==0)  SMADep = std::numeric_limits< double >::infinity( );
            if(i==1)  SMACap = std::numeric_limits< double >::infinity( );
            break;
            }
         }
    }

    SMAs << SMADep, SMACap;
    return SMAs;
}

Eigen::Matrix3d makeSkewed (Eigen::Vector3d axis) {
    Eigen::Matrix3d skewedMatrix;
    skewedMatrix.setZero();
    skewedMatrix << 0, -axis(2), axis(1),
        axis(2), 0, -axis(0),
        -axis(1), axis(0), 0;
    return skewedMatrix;

}


/*
Eigen::VectorXd getVariableVector (std::vector< LegTypeTransfer > LegTypeVec, std::vector< double > xv){
    std::cout << "So far so good\n";


    int numberOfLegs = LegTypeVec.size();
    int num_DSMs1 = std::count(LegTypeVec.begin(), LegTypeVec.end(), 3);
    int num_DSMs2 = std::count(LegTypeVec.begin(), LegTypeVec.end(), 4);
    int num_DSMs = num_DSMs1 + num_DSMs2;
    std::cout << "So far so good\n";

    Eigen::VectorXd variableVector(numberOfLegs+1+(4*num_DSMs));

    for(int i = 0; i < numberOfLegs ; i++){
        variableVector[ i ] = xv[ i ];

    }
    std::cout << "So far so good\n";

    variableVector[ numberOfLegs ] = 1;//dummy
    variableVector *= physical_constants::JULIAN_DAY;
    std::cout << "So far so good\n";

    //Write the variable vector for DSM variables below

    std::cout << "So far so good\n";

    if(num_DSMs>0){
        int count {numberOfLegs+1};
        for(size_t i{0}; i<num_DSMs; i++){
            variableVector[ count ] = xv[ count-1 ];
            variableVector[ count + 1 ] = xv[ count ];
            variableVector[ count + 2 ] = xv[ count + 1 ];
            variableVector[ count + 3 ] = xv[ count + 2 ];
            count+=4;
        }
    }

return variableVector;

}

*/
