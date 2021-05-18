#include <iostream>
#include <fstream>
#include <functional>

#include <boost/filesystem.hpp>
#include "G:/Thesis/applicationOutput.h"

#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"
#include "Classes/tt_functions.h"
#include "Classes/hodographicshapingproblem.h"
#include "Classes/trajectoryoptions.h"



#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/sphericalShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLegSettings.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/getRecommendedBaseFunctionsHodographicShaping.h"
#include "Tudat/SimulationSetup/optimisationSettings.h"


#include "Classes/tt_functions.h"
#include "Classes/lagrangepoints.h"
#include "Classes/trajectoryoptions.h"
#include "Classes/analysisb.h"
#include "G:/Thesis/applicationOutput.h"
#include "Classes/analysisb.h"
#include "Classes/analysisb_functions.h"

#include <fstream>

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat;
using namespace tudat::shape_based_methods;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;


int main() {

    spice_interface::loadStandardSpiceKernels();
    std::cout.precision(14);


    double Sun_GM = spice_interface::getBodyProperties("SUN", "GM").at(0)*1e9;
    double EB_GM = spice_interface::getBodyProperties("EARTH_BARYCENTER", "GM").at(0)*1e9;
    double Earth_GM = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
    double Moon_GM = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
    double Mars_GM = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;


    std::cout << "\n- SUN: " << Sun_GM;
    std::cout << "\n- EB: " << EB_GM;
    std::cout << "\n- Earth: " << Earth_GM;
    std::cout << "\n- Moon: " << Moon_GM;
    std::cout << "\n- Mars: " << Mars_GM;




    /*std::string Gateway = "G8";
    std::string planet = "Moon";
    std::string LPO_type = "halo";

    double arrivalTime = 12000 * physical_constants::JULIAN_DAY;
    double GM_stayTime = 1;

    //create Bounds Map
    std::map<std::string, std::pair<double, double>> boundsMap {};
    boundsMap["SMA"] = std::make_pair(300e3, 1000e3);
    boundsMap["Ecc"] = std::make_pair(0, 0.90);
    boundsMap["AOP"] = std::make_pair(0, 360);


    std::string filePath_prefix = tudat_applications::getOutputPath( ) + "Tester/CentralOrbits/";

    //Import orbits of interest
    int Nsteps = 3;
    std::unordered_map<int, std::map<double, Eigen::VectorXd>> orbitLibrary = createOrbitLibrary(Gateway, boundsMap, LPO_type, Nsteps);

    std::map<double,Eigen::Vector6d> GW_Mars_1 =  getCentralGateway_orbit (planet,  orbitLibrary[1],  0.0, arrivalTime, GM_stayTime, true );
    std::map<double,Eigen::Vector6d> GW_Mars_5 =  getCentralGateway_orbit (planet,  orbitLibrary[5],  0.0, arrivalTime, GM_stayTime, true );
    std::map<double,Eigen::Vector6d> GW_Mars_11 =  getCentralGateway_orbit (planet,  orbitLibrary[11],  0.0, arrivalTime, GM_stayTime, true );
    std::map<double,Eigen::Vector6d> GW_Mars_17 =  getCentralGateway_orbit (planet,  orbitLibrary[17],  0.0, arrivalTime, GM_stayTime, true );
    std::map<double,Eigen::Vector6d> GW_Mars_25 =  getCentralGateway_orbit (planet,  orbitLibrary[25],  0.0, arrivalTime, GM_stayTime, true );

    input_output::writeDataMapToTextFile(GW_Mars_1, "GW_Mars_1.dat", filePath_prefix );
    input_output::writeDataMapToTextFile(GW_Mars_5, "GW_Mars_5.dat", filePath_prefix );
    input_output::writeDataMapToTextFile(GW_Mars_11, "GW_Mars_11.dat", filePath_prefix );
    input_output::writeDataMapToTextFile(GW_Mars_17, "GW_Mars_17.dat", filePath_prefix );
    input_output::writeDataMapToTextFile(GW_Mars_25, "GW_Mars_25.dat", filePath_prefix );




    std::string planet = "Earth";
    std::pair<double, double> SMA_bounds = std::make_pair(100e3, 100000e3 );
    std::pair<double, double> eccentricity_bounds = std::make_pair(0, 0.90 );
    std::pair<double, double> AOP_bounds = std::make_pair(0, 360 );


    std::unordered_map<int, std::map<double, Eigen::VectorXd>> allCentralOrbits = generateCentralOrbits (planet,  SMA_bounds,  eccentricity_bounds, AOP_bounds );


    int orbitID = 60;
    double arrivalFraction = 0.3;
    double arrivalTime = 12000 * physical_constants::JULIAN_DAY;
    double stayTime = 9;
    double departureTime = arrivalTime + (stayTime * physical_constants::JULIAN_DAY);
    std::map<double, Eigen::VectorXd> specificCentralOrbit = allCentralOrbits[orbitID];


    //Planet frame
    Eigen::VectorXd arrivalState_planet = getCentralGateway_ArrState (planet, false,  specificCentralOrbit, arrivalFraction, arrivalTime);
    Eigen::VectorXd departureState_planet = getCentralGateway_DepState (planet, false,  specificCentralOrbit, arrivalFraction, departureTime, stayTime);
    std::map<double, Eigen::VectorXd> centralOrbit_planet = getCentralGateway_orbit(planet, specificCentralOrbit, arrivalFraction, arrivalTime, stayTime, false);

    std::map<int, Eigen::VectorXd> boundaryConditions_planet {};
    boundaryConditions_planet[1] = arrivalState_planet;
    boundaryConditions_planet[2] = departureState_planet;

    std::string filePath_prefix = tudat_applications::getOutputPath( ) + "Tester/CentralOrbits/";
    input_output::writeDataMapToTextFile(boundaryConditions_planet, "boundaryConditions_planet.dat", filePath_prefix );
    input_output::writeDataMapToTextFile(centralOrbit_planet, "centralOrbit_planet.dat", filePath_prefix );

    //Sun-centered frame
    Eigen::VectorXd arrivalState_Sun = getCentralGateway_ArrState (planet, true,  specificCentralOrbit, arrivalFraction, arrivalTime);
    Eigen::VectorXd departureState_Sun = getCentralGateway_DepState (planet, true,  specificCentralOrbit, arrivalFraction, departureTime, stayTime);
    std::map<double, Eigen::VectorXd> centralOrbit_Sun = getCentralGateway_orbit(planet, specificCentralOrbit, arrivalFraction, arrivalTime, stayTime, true);

    std::map<int, Eigen::VectorXd> boundaryConditions_Sun {};
    boundaryConditions_Sun[1] = arrivalState_Sun;
    boundaryConditions_Sun[2] = departureState_Sun;

    input_output::writeDataMapToTextFile(boundaryConditions_Sun, "boundaryConditions_Sun.dat", filePath_prefix );
    input_output::writeDataMapToTextFile(centralOrbit_Sun, "centralOrbit_Sun.dat", filePath_prefix );




    std::string system = "Em";
    std::string LPO_type = "halo";
    std::string LP = "L1";
    std::pair<int, int> orbitID_bounds = std::make_pair(2000, 2200);

    std::unordered_map<int, std::map<double, Eigen::VectorXd>> allLPOs = importLPO_cr3b (system, LPO_type, LP, orbitID_bounds);

    int orbitID = 2172;
    double arrivalFraction = 0.5046572823580753;
    double arrivalTime = 12857.99859774754 * physical_constants::JULIAN_DAY;
    double stayTime =  44.31622471933142; //days
    double departureTime = arrivalTime + (stayTime * physical_constants::JULIAN_DAY);

    std::map<double, Eigen::VectorXd> specificLPO = allLPOs[orbitID];
    std::string filePath_prefix = tudat_applications::getOutputPath( ) + "LPO-Library/Tester2/";
    input_output::writeDataMapToTextFile(specificLPO, "specificLPO.dat", filePath_prefix );

    Eigen::Vector6d Arr_state = getLPGateway_ArrState(system, specificLPO, arrivalFraction, arrivalTime );
    Eigen::Vector6d Dep_state = getLPGateway_DepState (system, specificLPO,  arrivalFraction,  departureTime,  stayTime ) ;


    std::map<int, Eigen::Vector6d> ArrDepStates;
    ArrDepStates[1] << Arr_state;
    ArrDepStates[2] << Dep_state;

    input_output::writeDataMapToTextFile(ArrDepStates, "ArrDepStates.dat", filePath_prefix );

    std::map<double, Eigen::Vector6d> Gateway_orbit = getLPGateway_orbit(system, specificLPO, arrivalFraction, arrivalTime, stayTime);
    input_output::writeDataMapToTextFile(Gateway_orbit, "Gateway_orbit.dat", filePath_prefix );




    std::string filePath_prefix = tudat_applications::getOutputPath( ) + "/Orbits/";


    std::map<double, Eigen::Vector6d> EarthOrbit = getBody_orbit("Earth", 12000 * physical_constants::JULIAN_DAY, "Sun", 1 );
    input_output::writeDataMapToTextFile(EarthOrbit, "EarthOrbit.dat", filePath_prefix );

    std::map<double, Eigen::Vector6d> MarsOrbit = getBody_orbit("Mars", 12000 * physical_constants::JULIAN_DAY, "Sun", 1 );
    input_output::writeDataMapToTextFile(MarsOrbit, "MarsOrbit.dat", filePath_prefix );

    std::map<double, Eigen::Vector6d> MoonOrbit = getBody_orbit("Moon", 12813.48708796495 * physical_constants::JULIAN_DAY, "Earth", 1 );
    input_output::writeDataMapToTextFile(MoonOrbit, "MoonOrbit.dat", filePath_prefix );



    //getDepstate
    Eigen::Vector6d  Gateway_DepState{};

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
         initialDistance = P1P2ephemeris[0]->getCartesianState( departureTime ).norm();

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

        initialDistance = P1P2ephemeris[0]->getCartesianState( departureTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Earth", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Moon", "GM").at(0)*1e9;
        }

    else if (system == "SM") {
        //distance
        std::vector<int> SunMars = {4};
        P1P2ephemeris = createEphemerisVec(SunMars, "Sun", "ECLIPJ2000", false);
        initialDistance = P1P2ephemeris[0]->getCartesianState( departureTime ).norm();

        //Gravitational parameters
        primGrav = spice_interface::getBodyProperties("Sun", "GM").at(0)*1e9;
        secGrav = spice_interface::getBodyProperties("Mars", "GM").at(0)*1e9;
        }
    else {
        std::cout << "\nError - system not recognized in FUNCTION: getLPO_cart. \n";
    }


    //Get Gateway period in days
    std::map<double, Eigen::Vector6d>::reverse_iterator  it_last = specificLPO.rbegin();
    double LP_period = tudat::circular_restricted_three_body_problem::convertDimensionlessTimeToDimensionalTime(it_last->first, primGrav, secGrav, initialDistance ) / physical_constants::JULIAN_DAY;

    double DepFraction = (orbitFraction + (stayTime/LP_period));
    std::cout << "\nDepFraction:  " << DepFraction;
    while(DepFraction>1.0) { DepFraction = DepFraction-1; }
    std::cout << "\nDepFraction:  " << DepFraction;

    int timeStepofInterest_Dep = round(DepFraction * specificLPO.size());

    std::cout << "\ntimeStepofInterest_Dep:  " << timeStepofInterest_Dep;


    //int advancement_Dep = (timeStepofInterest_Dep!=0) ? (timeStepofInterest_Dep-1) : timeStepofInterest_Dep;
    std::map<double, Eigen::Vector6d>::iterator it_dep = specificLPO.begin();
    std::advance(it_dep, timeStepofInterest_Dep);
    Eigen::Vector6d stateDep_cr3b = it_dep->second;

    if (system == "Em") {
        Eigen::Vector6d state_cart_E = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, stateDep_cr3b, P1P2ephemeris[0]->getCartesianState(departureTime));
        Gateway_DepState = state_cart_E + SunEarthephemeris[0]->getCartesianState(departureTime);
    }
    else {
        Gateway_DepState = tudat::circular_restricted_three_body_problem::convertCorotatingNormalizedToCartesianCoordinates_thesis(primGrav, secGrav, stateDep_cr3b, P1P2ephemeris[0]->getCartesianState(departureTime));

    }

*/

    std::cout << "\n\n------------------------\nProgram ended succesfully!\n--------------------------\n\n";
    return EXIT_SUCCESS;
}



