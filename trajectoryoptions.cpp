#include "trajectoryoptions.h"

TrajectoryOptions::TrajectoryOptions(std::string Gateway, std::string Segment, bool useDCT, bool useDSM,  bool useMGA, bool useMGADSM):
    Gateway_{Gateway}, Segment_{Segment} ,useDCT_{useDCT}, useDSM_{useDSM}, useMGA_{useMGA}, useMGADSM_{useMGADSM}
{};

std::vector<std::string> TrajectoryOptions::makeStringVector(){

    if(Gateway_=="All"){
        std::vector<std::string> AllOptions{};
        std::vector<std::string> G1options = TrajectoryOptions::makeG1options( );
        std::vector<std::string> G2options = TrajectoryOptions::makeG2options( );
        std::vector<std::string> G3options = TrajectoryOptions::makeG3options( );
        std::vector<std::string> G4options = TrajectoryOptions::makeG4options( );

        /*
        std::vector<std::string> G5options = TrajectoryOptions::makeG5options( );
        std::vector<std::string> G6options = TrajectoryOptions::makeG6options( );
        std::vector<std::string> G7options = TrajectoryOptions::makeG7options( );
        std::vector<std::string> G8options = TrajectoryOptions::makeG8options( );
        std::vector<std::string> G9options = TrajectoryOptions::makeG9options( );
}
*/
        AllOptions.insert( AllOptions.end(), G1options.begin(), G1options.end() );
        AllOptions.insert( AllOptions.end(), G2options.begin(), G2options.end() );
        AllOptions.insert( AllOptions.end(), G3options.begin(), G3options.end() );
        AllOptions.insert( AllOptions.end(), G4options.begin(), G4options.end() );
        if(Segment_=="GG"){std::vector<std::string> GGoptions = TrajectoryOptions::makeGGoptions( );
                      AllOptions.insert( AllOptions.end(), GGoptions.begin(), GGoptions.end() );}



        /*
        AllOptions.insert( AllOptions.end(), G5options.begin(), G5options.end() );
        AllOptions.insert( AllOptions.end(), G6options.begin(), G6options.end() );
        AllOptions.insert( AllOptions.end(), G7options.begin(), G7options.end() );
        AllOptions.insert( AllOptions.end(), G8options.begin(), G8options.end() );
        AllOptions.insert( AllOptions.end(), G9options.begin(), G9options.end() );
*/
        return AllOptions;
    }

    if(Segment_=="GG"){
        std::vector<std::string> GGoptions = TrajectoryOptions::makeGGoptions( );
        return GGoptions;
    }

    if(Gateway_=="G1"){
        std::vector<std::string> G1options = TrajectoryOptions::makeG1options( );
        return G1options;
    }

    if(Gateway_=="G2"){
        std::vector<std::string> G2options = TrajectoryOptions::makeG2options( );
        return G2options;
    }
    if(Gateway_=="G3"){
        std::vector<std::string> G3options = TrajectoryOptions::makeG3options( );
        return G3options;
    }
    if(Gateway_=="G4"){
        std::vector<std::string> G4options = TrajectoryOptions::makeG4options( );
        return G4options;
    }
    if(Gateway_=="G5"){
        std::vector<std::string> G5options = TrajectoryOptions::makeG5options( );
        return G5options;
    }
    if(Gateway_=="G6"){
        std::vector<std::string> G6options = TrajectoryOptions::makeG6options( );
        return G6options;
    }

    /*
    if(Gateway_=="G7"){
        std::vector<std::string> G7options = TrajectoryOptions::makeG7options( );
        return G7options;
    }
    if(Gateway_=="G8"){
        std::vector<std::string> G8options = TrajectoryOptions::makeG8options();
        return G8options;
    }
    if(Gateway_=="G9"){
        std::vector<std::string> G9options = TrajectoryOptions::makeG9options( );
        return G9options;
    }
    */
    else {std::cout<< "Error: Gateway input not recognized. ";}


}


std::vector<std::string> TrajectoryOptions::makeG1options( )
{

        std::vector<std::string> G1Options{};

        if(Segment_=="All"){
            std::vector<std::string> EG_G1Options{};


        if(useDCT_){ std::vector<std::string> EG_G1_DCT = {"EG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_DCT.begin(), EG_G1_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G1_DSM = {"EdG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_DSM.begin(), EG_G1_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G1_MGA = {"EmG1", "EEG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_MGA.begin(), EG_G1_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G1_MGADSM = {"EdmdG1", "EdEdG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_MGADSM.begin(), EG_G1_MGADSM.end() );
                     }

        G1Options.insert( G1Options.end(), EG_G1Options.begin(), EG_G1Options.end() );
        std::vector<std::string> GM_G1Options{};


    if(useDCT_){ std::vector<std::string> GM_G1_DCT = {"G1M"};
                 GM_G1Options.insert( GM_G1Options.end(), GM_G1_DCT.begin(), GM_G1_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> GM_G1_DSM = {"G1dM"};
                 GM_G1Options.insert( GM_G1Options.end(), GM_G1_DSM.begin(), GM_G1_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> GM_G1_MGA = {"G1mM", "G1EM", "G1mEM", "G1EmM"};
                 GM_G1Options.insert( GM_G1Options.end(), GM_G1_MGA.begin(), GM_G1_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> GM_G1_MGADSM = {"G1dmdM", "G1dEdM", "G1dmdEdM", "G1dEdmdM"};
                 GM_G1Options.insert( GM_G1Options.end(), GM_G1_MGADSM.begin(), GM_G1_MGADSM.end() );
                 }

    G1Options.insert( G1Options.end(), GM_G1Options.begin(), GM_G1Options.end() );
    std::vector<std::string> GE_G1Options{};


if(useDCT_){ std::vector<std::string> GE_G1_DCT = {"G1E"};
             GE_G1Options.insert( GE_G1Options.end(), GE_G1_DCT.begin(), GE_G1_DCT.end() );
             }
if(useDSM_){ std::vector<std::string> GE_G1_DSM = {"G1dE"};
             GE_G1Options.insert( GE_G1Options.end(), GE_G1_DSM.begin(), GE_G1_DSM.end() );
             }
if(useMGA_){ std::vector<std::string> GE_G1_MGA = {"G1mE", "G1EE"};
             GE_G1Options.insert( GE_G1Options.end(), GE_G1_MGA.begin(), GE_G1_MGA.end() );
             }

if(useMGADSM_){ std::vector<std::string> GE_G1_MGADSM = {"G1dmdE", "G1dEdE"};
             GE_G1Options.insert( GE_G1Options.end(), GE_G1_MGADSM.begin(), GE_G1_MGADSM.end() );
             }

G1Options.insert( G1Options.end(), GE_G1Options.begin(), GE_G1Options.end() );
std::vector<std::string> MG_G1Options{};


if(useDCT_){ std::vector<std::string> MG_G1_DCT = {"MG1"};
         MG_G1Options.insert( MG_G1Options.end(), MG_G1_DCT.begin(), MG_G1_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> MG_G1_DSM = {"MdG1"};
         MG_G1Options.insert( MG_G1Options.end(), MG_G1_DSM.begin(), MG_G1_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> MG_G1_MGA = {"MmG1", "MEG1", "MmEG1", "MEmG1"};
         MG_G1Options.insert( MG_G1Options.end(), MG_G1_MGA.begin(), MG_G1_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> MG_G1_MGADSM = {"MdmdG1", "MdEdG1", "MdmdEdG1", "MdEdmdG1"};
         MG_G1Options.insert( MG_G1Options.end(), MG_G1_MGADSM.begin(), MG_G1_MGADSM.end() );
         }

G1Options.insert( G1Options.end(), MG_G1Options.begin(), MG_G1Options.end() );

            return G1Options;
        }

        if(Segment_=="EG"){
            std::vector<std::string> EG_G1Options{};


        if(useDCT_){ std::vector<std::string> EG_G1_DCT = {"EG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_DCT.begin(), EG_G1_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G1_DSM = {"EdG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_DSM.begin(), EG_G1_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G1_MGA = {"EmG1", "EEG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_MGA.begin(), EG_G1_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G1_MGADSM = {"EdmdG1", "EdEdG1"};
                     EG_G1Options.insert( EG_G1Options.end(), EG_G1_MGADSM.begin(), EG_G1_MGADSM.end() );
                     }

        G1Options.insert( G1Options.end(), EG_G1Options.begin(), EG_G1Options.end() );

        }

        if(Segment_=="GM"){
            std::vector<std::string> GM_G1Options{};


        if(useDCT_){ std::vector<std::string> GM_G1_DCT = {"G1M"};
                     GM_G1Options.insert( GM_G1Options.end(), GM_G1_DCT.begin(), GM_G1_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GM_G1_DSM = {"G1dM"};
                     GM_G1Options.insert( GM_G1Options.end(), GM_G1_DSM.begin(), GM_G1_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GM_G1_MGA = {"G1mM", "G1EM", "G1mEM", "G1EmM"};
                     GM_G1Options.insert( GM_G1Options.end(), GM_G1_MGA.begin(), GM_G1_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GM_G1_MGADSM = {"G1dmdM", "G1dEdM", "G1dmdEdM", "G1dEdmdM"};
                     GM_G1Options.insert( GM_G1Options.end(), GM_G1_MGADSM.begin(), GM_G1_MGADSM.end() );
                     }

        G1Options.insert( G1Options.end(), GM_G1Options.begin(), GM_G1Options.end() );

        }

        if(Segment_=="GE"){
            std::vector<std::string> GE_G1Options{};


        if(useDCT_){ std::vector<std::string> GE_G1_DCT = {"G1E"};
                     GE_G1Options.insert( GE_G1Options.end(), GE_G1_DCT.begin(), GE_G1_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GE_G1_DSM = {"G1dE"};
                     GE_G1Options.insert( GE_G1Options.end(), GE_G1_DSM.begin(), GE_G1_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GE_G1_MGA = {"G1mE", "G1EE"};
                     GE_G1Options.insert( GE_G1Options.end(), GE_G1_MGA.begin(), GE_G1_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GE_G1_MGADSM = {"G1dmdE", "G1dEdE"};
                     GE_G1Options.insert( GE_G1Options.end(), GE_G1_MGADSM.begin(), GE_G1_MGADSM.end() );
                     }

        G1Options.insert( G1Options.end(), GE_G1Options.begin(), GE_G1Options.end() );

        }


        if(Segment_=="MG"){
            std::vector<std::string> MG_G1Options{};


        if(useDCT_){ std::vector<std::string> MG_G1_DCT = {"MG1"};
                     MG_G1Options.insert( MG_G1Options.end(), MG_G1_DCT.begin(), MG_G1_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> MG_G1_DSM = {"MdG1"};
                     MG_G1Options.insert( MG_G1Options.end(), MG_G1_DSM.begin(), MG_G1_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> MG_G1_MGA = {"MmG1", "MEG1", "MmEG1", "MEmG1"};
                     MG_G1Options.insert( MG_G1Options.end(), MG_G1_MGA.begin(), MG_G1_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> MG_G1_MGADSM = {"MdmdG1", "MdEdG1", "MdmdEdG1", "MdEdmdG1"};
                     MG_G1Options.insert( MG_G1Options.end(), MG_G1_MGADSM.begin(), MG_G1_MGADSM.end() );
                     }

        G1Options.insert( G1Options.end(), MG_G1Options.begin(), MG_G1Options.end() );

        }


        return G1Options;


}
std::vector<std::string> TrajectoryOptions::makeG2options( )
{
        std::vector<std::string> G2Options{};

        if(Segment_=="EG"){
            std::vector<std::string> EG_G2Options{};


        if(useDCT_){ std::vector<std::string> EG_G2_DCT = {"EG2"};
                     EG_G2Options.insert( EG_G2Options.end(), EG_G2_DCT.begin(), EG_G2_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G2_DSM = {"EdG2"};
                     EG_G2Options.insert( EG_G2Options.end(), EG_G2_DSM.begin(), EG_G2_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G2_MGA = {"EmG2", "EEG2"};
                     EG_G2Options.insert( EG_G2Options.end(), EG_G2_MGA.begin(), EG_G2_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G2_MGADSM = {"EdmdG2", "EdEdG2"};
                     EG_G2Options.insert( EG_G2Options.end(), EG_G2_MGADSM.begin(), EG_G2_MGADSM.end() );
                     }

        G2Options.insert( G2Options.end(), EG_G2Options.begin(), EG_G2Options.end() );

        }

        if(Segment_=="GM"){
            std::vector<std::string> GM_G2Options{};


        if(useDCT_){ std::vector<std::string> GM_G2_DCT = {"G2M"};
                     GM_G2Options.insert( GM_G2Options.end(), GM_G2_DCT.begin(), GM_G2_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GM_G2_DSM = {"G2dM"};
                     GM_G2Options.insert( GM_G2Options.end(), GM_G2_DSM.begin(), GM_G2_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GM_G2_MGA = {"G2mM", "G2EM", "G2mEM", "G2EmM"};
                     GM_G2Options.insert( GM_G2Options.end(), GM_G2_MGA.begin(), GM_G2_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GM_G2_MGADSM = {"G2dmdM", "G2dEdM", "G2dmdEdM", "G2dEdmdM"};
                     GM_G2Options.insert( GM_G2Options.end(), GM_G2_MGADSM.begin(), GM_G2_MGADSM.end() );
                     }

        G2Options.insert( G2Options.end(), GM_G2Options.begin(), GM_G2Options.end() );

        }

        if(Segment_=="GE"){
            std::vector<std::string> GE_G2Options{};


        if(useDCT_){ std::vector<std::string> GE_G2_DCT = {"G2E"};
                     GE_G2Options.insert( GE_G2Options.end(), GE_G2_DCT.begin(), GE_G2_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GE_G2_DSM = {"G2dE"};
                     GE_G2Options.insert( GE_G2Options.end(), GE_G2_DSM.begin(), GE_G2_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GE_G2_MGA = {"G2mE", "G2EE"};
                     GE_G2Options.insert( GE_G2Options.end(), GE_G2_MGA.begin(), GE_G2_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GE_G2_MGADSM = {"G2dmdE", "G2dEdE"};
                     GE_G2Options.insert( GE_G2Options.end(), GE_G2_MGADSM.begin(), GE_G2_MGADSM.end() );
                     }

        G2Options.insert( G2Options.end(), GE_G2Options.begin(), GE_G2Options.end() );

        }


        if(Segment_=="MG"){
            std::vector<std::string> MG_G2Options{};


        if(useDCT_){ std::vector<std::string> MG_G2_DCT = {"MG2"};
                     MG_G2Options.insert( MG_G2Options.end(), MG_G2_DCT.begin(), MG_G2_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> MG_G2_DSM = {"MdG2"};
                     MG_G2Options.insert( MG_G2Options.end(), MG_G2_DSM.begin(), MG_G2_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> MG_G2_MGA = {"MmG2", "MEG2", "MmEG2", "MEmG2"};
                     MG_G2Options.insert( MG_G2Options.end(), MG_G2_MGA.begin(), MG_G2_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> MG_G2_MGADSM = {"MdmdG2", "MdEdG2", "MdmdEdG2", "MdEdmdG2"};
                     MG_G2Options.insert( MG_G2Options.end(), MG_G2_MGADSM.begin(), MG_G2_MGADSM.end() );
                     }

        G2Options.insert( G2Options.end(), MG_G2Options.begin(), MG_G2Options.end() );

        }


        return G2Options;

}
std::vector<std::string> TrajectoryOptions::makeG3options( )
{
        std::vector<std::string> G3Options{};

        if(Segment_=="EG"){
            std::vector<std::string> EG_G3Options{};


        if(useDCT_){ std::vector<std::string> EG_G3_DCT = {"EG3"};
                     EG_G3Options.insert( EG_G3Options.end(), EG_G3_DCT.begin(), EG_G3_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G3_DSM = {"EdG3"};
                     EG_G3Options.insert( EG_G3Options.end(), EG_G3_DSM.begin(), EG_G3_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G3_MGA = {"EmG3", "EEG3"};
                     EG_G3Options.insert( EG_G3Options.end(), EG_G3_MGA.begin(), EG_G3_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G3_MGADSM = {"EdmdG3", "EdEdG3"};
                     EG_G3Options.insert( EG_G3Options.end(), EG_G3_MGADSM.begin(), EG_G3_MGADSM.end() );
                     }

        G3Options.insert( G3Options.end(), EG_G3Options.begin(), EG_G3Options.end() );

        }

        if(Segment_=="GM"){
            std::vector<std::string> GM_G3Options{};


        if(useDCT_){ std::vector<std::string> GM_G3_DCT = {"G3M"};
                     GM_G3Options.insert( GM_G3Options.end(), GM_G3_DCT.begin(), GM_G3_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GM_G3_DSM = {"G3dM"};
                     GM_G3Options.insert( GM_G3Options.end(), GM_G3_DSM.begin(), GM_G3_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GM_G3_MGA = {"G3mM", "G3EM", "G3mEM", "G3EmM"};
                     GM_G3Options.insert( GM_G3Options.end(), GM_G3_MGA.begin(), GM_G3_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GM_G3_MGADSM = {"G3dmdM", "G3dEdM", "G3dmdEdM", "G3dEdmdM"};
                     GM_G3Options.insert( GM_G3Options.end(), GM_G3_MGADSM.begin(), GM_G3_MGADSM.end() );
                     }

        G3Options.insert( G3Options.end(), GM_G3Options.begin(), GM_G3Options.end() );

        }

        if(Segment_=="GE"){
            std::vector<std::string> GE_G3Options{};


        if(useDCT_){ std::vector<std::string> GE_G3_DCT = {"G3E"};
                     GE_G3Options.insert( GE_G3Options.end(), GE_G3_DCT.begin(), GE_G3_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GE_G3_DSM = {"G3dE"};
                     GE_G3Options.insert( GE_G3Options.end(), GE_G3_DSM.begin(), GE_G3_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GE_G3_MGA = {"G3mE", "G3EE"};
                     GE_G3Options.insert( GE_G3Options.end(), GE_G3_MGA.begin(), GE_G3_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GE_G3_MGADSM = {"G3dmdE", "G3dEdE"};
                     GE_G3Options.insert( GE_G3Options.end(), GE_G3_MGADSM.begin(), GE_G3_MGADSM.end() );
                     }

        G3Options.insert( G3Options.end(), GE_G3Options.begin(), GE_G3Options.end() );

        }


        if(Segment_=="MG"){
            std::vector<std::string> MG_G3Options{};


        if(useDCT_){ std::vector<std::string> MG_G3_DCT = {"MG3"};
                     MG_G3Options.insert( MG_G3Options.end(), MG_G3_DCT.begin(), MG_G3_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> MG_G3_DSM = {"MdG3"};
                     MG_G3Options.insert( MG_G3Options.end(), MG_G3_DSM.begin(), MG_G3_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> MG_G3_MGA = {"MmG3", "MEG3", "MmEG3", "MEmG3"};
                     MG_G3Options.insert( MG_G3Options.end(), MG_G3_MGA.begin(), MG_G3_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> MG_G3_MGADSM = {"MdmdG3", "MdEdG3", "MdmdEdG3", "MdEdmdG3"};
                     MG_G3Options.insert( MG_G3Options.end(), MG_G3_MGADSM.begin(), MG_G3_MGADSM.end() );
                     }

        G3Options.insert( G3Options.end(), MG_G3Options.begin(), MG_G3Options.end() );

        }


        return G3Options;

}
std::vector<std::string> TrajectoryOptions::makeG4options()
{
        std::vector<std::string> G4Options{};

        if(Segment_=="EG"){
            std::vector<std::string> EG_G4Options{};


        if(useDCT_){ std::vector<std::string> EG_G4_DCT = {"EG4"};
                     EG_G4Options.insert( EG_G4Options.end(), EG_G4_DCT.begin(), EG_G4_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G4_DSM = {"EdG4"};
                     EG_G4Options.insert( EG_G4Options.end(), EG_G4_DSM.begin(), EG_G4_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G4_MGA = {"EmG4", "EEG4"};
                     EG_G4Options.insert( EG_G4Options.end(), EG_G4_MGA.begin(), EG_G4_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G4_MGADSM = {"EdmdG4", "EdEdG4"};
                     EG_G4Options.insert( EG_G4Options.end(), EG_G4_MGADSM.begin(), EG_G4_MGADSM.end() );
                     }

        G4Options.insert( G4Options.end(), EG_G4Options.begin(), EG_G4Options.end() );

        }

        if(Segment_=="GM"){
            std::vector<std::string> GM_G4Options{};


        if(useDCT_){ std::vector<std::string> GM_G4_DCT = {"G4M"};
                     GM_G4Options.insert( GM_G4Options.end(), GM_G4_DCT.begin(), GM_G4_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GM_G4_DSM = {"G4dM"};
                     GM_G4Options.insert( GM_G4Options.end(), GM_G4_DSM.begin(), GM_G4_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GM_G4_MGA = {"G4mM", "G4EM", "G4mEM", "G4EmM"};
                     GM_G4Options.insert( GM_G4Options.end(), GM_G4_MGA.begin(), GM_G4_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GM_G4_MGADSM = {"G4dmdM", "G4dEdM", "G4dmdEdM", "G4dEdmdM"};
                     GM_G4Options.insert( GM_G4Options.end(), GM_G4_MGADSM.begin(), GM_G4_MGADSM.end() );
                     }

        G4Options.insert( G4Options.end(), GM_G4Options.begin(), GM_G4Options.end() );

        }

        if(Segment_=="GE"){
            std::vector<std::string> GE_G4Options{};


        if(useDCT_){ std::vector<std::string> GE_G4_DCT = {"G4E"};
                     GE_G4Options.insert( GE_G4Options.end(), GE_G4_DCT.begin(), GE_G4_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GE_G4_DSM = {"G4dE"};
                     GE_G4Options.insert( GE_G4Options.end(), GE_G4_DSM.begin(), GE_G4_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GE_G4_MGA = {"G4mE", "G4EE"};
                     GE_G4Options.insert( GE_G4Options.end(), GE_G4_MGA.begin(), GE_G4_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GE_G4_MGADSM = {"G4dmdE", "G4dEdE"};
                     GE_G4Options.insert( GE_G4Options.end(), GE_G4_MGADSM.begin(), GE_G4_MGADSM.end() );
                     }

        G4Options.insert( G4Options.end(), GE_G4Options.begin(), GE_G4Options.end() );

        }


        if(Segment_=="MG"){
            std::vector<std::string> MG_G4Options{};


        if(useDCT_){ std::vector<std::string> MG_G4_DCT = {"MG4"};
                     MG_G4Options.insert( MG_G4Options.end(), MG_G4_DCT.begin(), MG_G4_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> MG_G4_DSM = {"MdG4"};
                     MG_G4Options.insert( MG_G4Options.end(), MG_G4_DSM.begin(), MG_G4_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> MG_G4_MGA = {"MmG4", "MEG4", "MmEG4", "MEmG4"};
                     MG_G4Options.insert( MG_G4Options.end(), MG_G4_MGA.begin(), MG_G4_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> MG_G4_MGADSM = {"MdmdG4", "MdEdG4", "MdmdEdG4", "MdEdmdG4"};
                     MG_G4Options.insert( MG_G4Options.end(), MG_G4_MGADSM.begin(), MG_G4_MGADSM.end() );
                     }

        G4Options.insert( G4Options.end(), MG_G4Options.begin(), MG_G4Options.end() );

        }


        return G4Options;

}

std::vector<std::string> TrajectoryOptions::makeG5options()
{
        std::vector<std::string> G5Options{};

        if(Segment_=="EG"){
            std::vector<std::string> EG_G5Options{};


        if(useDCT_){ std::vector<std::string> EG_G5_DCT = {"EG5"};
                     EG_G5Options.insert( EG_G5Options.end(), EG_G5_DCT.begin(), EG_G5_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G5_DSM = {"EdG5"};
                     EG_G5Options.insert( EG_G5Options.end(), EG_G5_DSM.begin(), EG_G5_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G5_MGA = {"EmG5", "EEG5", "EmEG5", "EmmG5", "EEmG5"};
                     EG_G5Options.insert( EG_G5Options.end(), EG_G5_MGA.begin(), EG_G5_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G5_MGADSM = {"EdmdG5", "EdEdG5", "EdmdEdG5", "EdmdmdG5", "EdEdmdG5"};
                     EG_G5Options.insert( EG_G5Options.end(), EG_G5_MGADSM.begin(), EG_G5_MGADSM.end() );
                     }

        G5Options.insert( G5Options.end(), EG_G5Options.begin(), EG_G5Options.end() );

        }

        if(Segment_=="GM"){
            std::vector<std::string> GM_G5Options{};


        if(useDCT_){ std::vector<std::string> GM_G5_DCT = {"G5M"};
                     GM_G5Options.insert( GM_G5Options.end(), GM_G5_DCT.begin(), GM_G5_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GM_G5_DSM = {"G5dM"};
                     GM_G5Options.insert( GM_G5Options.end(), GM_G5_DSM.begin(), GM_G5_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GM_G5_MGA = {"G5mM", "G5EM", "G5mEM", "G5EmM", "G5mEmM", "G5EmEM", "G5EmmM", "G5EEmM", "G5mEEM"};
                     GM_G5Options.insert( GM_G5Options.end(), GM_G5_MGA.begin(), GM_G5_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GM_G5_MGADSM = {"G5dmdM", "G5dEdM", "G5dmdEdM", "G5dEdmdM", "G5dmdEdmdM", "G5dEdmdEdM", "G5dEdmdmdM", "G5dEdEdmdM", "G5dmdEdEM"};
                     GM_G5Options.insert( GM_G5Options.end(), GM_G5_MGADSM.begin(), GM_G5_MGADSM.end() );
                     }

        G5Options.insert( G5Options.end(), GM_G5Options.begin(), GM_G5Options.end() );

        }

        if(Segment_=="GE"){
            std::vector<std::string> GE_G5Options{};


        if(useDCT_){ std::vector<std::string> GE_G5_DCT = {"G5E"};
                     GE_G5Options.insert( GE_G5Options.end(), GE_G5_DCT.begin(), GE_G5_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GE_G5_DSM = {"G5dE"};
                     GE_G5Options.insert( GE_G5Options.end(), GE_G5_DSM.begin(), GE_G5_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GE_G5_MGA = {"G5mE", "G5EE", "G5mEE", "G5mmE", "G5EmE"};
                     GE_G5Options.insert( GE_G5Options.end(), GE_G5_MGA.begin(), GE_G5_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GE_G5_MGADSM = {"G5dmdE", "G5dEdE", "G5dmdEdE", "G5dmdmdE", "G5dEdmdE"};
                     GE_G5Options.insert( GE_G5Options.end(), GE_G5_MGADSM.begin(), GE_G5_MGADSM.end() );
                     }

        G5Options.insert( G5Options.end(), GE_G5Options.begin(), GE_G5Options.end() );

        }


        if(Segment_=="MG"){
            std::vector<std::string> MG_G5Options{};


        if(useDCT_){ std::vector<std::string> MG_G5_DCT = {"MG5"};
                     MG_G5Options.insert( MG_G5Options.end(), MG_G5_DCT.begin(), MG_G5_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> MG_G5_DSM = {"MdG5"};
                     MG_G5Options.insert( MG_G5Options.end(), MG_G5_DSM.begin(), MG_G5_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> MG_G5_MGA = {"MmG5", "MEG5", "MmEG5", "MEmG5"};
                     MG_G5Options.insert( MG_G5Options.end(), MG_G5_MGA.begin(), MG_G5_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> MG_G5_MGADSM = {"MdmdG5", "MdEdG5", "MdmdEdG5", "MdEdmdG5"};
                     MG_G5Options.insert( MG_G5Options.end(), MG_G5_MGADSM.begin(), MG_G5_MGADSM.end() );
                     }

        G5Options.insert( G5Options.end(), MG_G5Options.begin(), MG_G5Options.end() );

        }


        return G5Options;

}

std::vector<std::string> TrajectoryOptions::makeG6options()
{
        std::vector<std::string> G6Options{};

        if(Segment_=="EG"){
            std::vector<std::string> EG_G6Options{};


        if(useDCT_){ std::vector<std::string> EG_G6_DCT = {"EG6"};
                     EG_G6Options.insert( EG_G6Options.end(), EG_G6_DCT.begin(), EG_G6_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> EG_G6_DSM = {"EdG6"};
                     EG_G6Options.insert( EG_G6Options.end(), EG_G6_DSM.begin(), EG_G6_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> EG_G6_MGA = {"EmG6", "EEG6", "EmEG6", "EEmG6", "EmEmG6", "EEmEG6", "EEmmG6", "EEEmG6", "EmEEG6", "EMG6", "EmMG6", "EEMG6"};
                     EG_G6Options.insert( EG_G6Options.end(), EG_G6_MGA.begin(), EG_G6_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> EG_G6_MGADSM =  {"EdmdG6", "EdEdG6", "EdmdEdG6", "EdEdmdG6", "EdmdEdmdG6", "EdEdmdEdG6", "EdEdmdmdG6", "EdEdEdmdG6", "EdmdEdEdG6", "EdMdG6", "EdmdMdG6", "EdEdMdG6"}; //"EdmdG6", "EdEdG6", "EdmdEdG6", "EdEdmdG6", "EdmdEdmdG6", NOG DOEN: "EdEdmdEdG6", "EdEdmdmdG6", "EdEdEdmdG6", "EdmdEdEdG6"
                     EG_G6Options.insert( EG_G6Options.end(), EG_G6_MGADSM.begin(), EG_G6_MGADSM.end() );
                     }

        G6Options.insert( G6Options.end(), EG_G6Options.begin(), EG_G6Options.end() );

        }

        if(Segment_=="GM"){
            std::vector<std::string> GM_G6Options{};


        if(useDCT_){ std::vector<std::string> GM_G6_DCT = {"G6M"};
                     GM_G6Options.insert( GM_G6Options.end(), GM_G6_DCT.begin(), GM_G6_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GM_G6_DSM = {"G6dM"};
                     GM_G6Options.insert( GM_G6Options.end(), GM_G6_DSM.begin(), GM_G6_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GM_G6_MGA = {"G6MM"};
                     GM_G6Options.insert( GM_G6Options.end(), GM_G6_MGA.begin(), GM_G6_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GM_G6_MGADSM = {"G6dMdM"};
                     GM_G6Options.insert( GM_G6Options.end(), GM_G6_MGADSM.begin(), GM_G6_MGADSM.end() );
                     }

        G6Options.insert( G6Options.end(), GM_G6Options.begin(), GM_G6Options.end() );

        }

        if(Segment_=="GE"){
            std::vector<std::string> GE_G6Options{};


        if(useDCT_){ std::vector<std::string> GE_G6_DCT = {"G6E"};
                     GE_G6Options.insert( GE_G6Options.end(), GE_G6_DCT.begin(), GE_G6_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> GE_G6_DSM = {"G6dE"};
                     GE_G6Options.insert( GE_G6Options.end(), GE_G6_DSM.begin(), GE_G6_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> GE_G6_MGA = {"G6mE", "G6EE", "G6mEE", "G6EmE", "G6mEmE", "G6EmEE", "G6EmmE", "G6EEmE", "G6ME", "G6MmE", "G6MEE", "G6MmEE"};
                     GE_G6Options.insert( GE_G6Options.end(), GE_G6_MGA.begin(), GE_G6_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> GE_G6_MGADSM = {"G6dmdE", "G6dEdE", "G6dmdEdE", "G6dEdmdE", "G6dmdEdmE", "G6dEdmdEdE", "G6dEdmdmdE", "G6dEdEdmdE", "G6dMdE", "G6dMdmdE", "G6dMdEdE", "G6dMdmdEdE"};
                     GE_G6Options.insert( GE_G6Options.end(), GE_G6_MGADSM.begin(), GE_G6_MGADSM.end() );
                     }

        G6Options.insert( G6Options.end(), GE_G6Options.begin(), GE_G6Options.end() );

        }


        if(Segment_=="MG"){
            std::vector<std::string> MG_G6Options{};


        if(useDCT_){ std::vector<std::string> MG_G6_DCT = {"MG6"};
                     MG_G6Options.insert( MG_G6Options.end(), MG_G6_DCT.begin(), MG_G6_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> MG_G6_DSM = {"MdG6"};
                     MG_G6Options.insert( MG_G6Options.end(), MG_G6_DSM.begin(), MG_G6_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> MG_G6_MGA = {"MMG6"};
                     MG_G6Options.insert( MG_G6Options.end(), MG_G6_MGA.begin(), MG_G6_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> MG_G6_MGADSM = {"MdMdG6"};
                     MG_G6Options.insert( MG_G6Options.end(), MG_G6_MGADSM.begin(), MG_G6_MGADSM.end() );
                     }

        G6Options.insert( G6Options.end(), MG_G6Options.begin(), MG_G6Options.end() );

        }


        return G6Options;

}





std::vector<std::string> TrajectoryOptions::makeGGoptions(){

    std::vector<std::string> GGOptions{};

    if(Gateway_=="All"){
        std::vector<std::string> G1_GGOptions{};


    if(useDCT_){ std::vector<std::string> G1_GG_DCT = {"G1G6", "G1G9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_DCT.begin(), G1_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G1_GG_DSM = {"G1dG6", "G1dG9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_DSM.begin(), G1_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G1_GG_MGA = {"G1mG6", "G1EG6", "G1mEG6", "G1EmG6", "G1mG9", "G1EG9", "G1mEG9", "G1EmG9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_MGA.begin(), G1_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G1_GG_MGADSM = {"G1dmdG6", "G1dEdG6", "G1dmdEdG6", "G1dEdmdG6", "G1dmdG9", "G1dEdG9", "G1dmdEdG9", "G1dEdmdG9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_MGADSM.begin(), G1_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G1_GGOptions.begin(), G1_GGOptions.end() );

    std::vector<std::string> G2_GGOptions{};


if(useDCT_){ std::vector<std::string> G2_GG_DCT = {"G2G6", "G2G9"};
             G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_DCT.begin(), G2_GG_DCT.end() );
             }
if(useDSM_){ std::vector<std::string> G2_GG_DSM = {"G2dG6", "G2dG9"};
             G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_DSM.begin(), G2_GG_DSM.end() );
             }
if(useMGA_){ std::vector<std::string> G2_GG_MGA = {"G2mG6", "G2EG6", "G2mEG6", "G2EmG6", "G2mG9", "G2EG9", "G2mEG9", "G2EmG9"};
             G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_MGA.begin(), G2_GG_MGA.end() );
             }

if(useMGADSM_){ std::vector<std::string> G2_GG_MGADSM = {"G2dmdG6", "G2dEdG6", "G2dmdEdG6", "G2dEdmdG6", "G2dmdG9", "G2dEdG9", "G2dmdEdG9", "G2dEdmdG9"};
             G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_MGADSM.begin(), G2_GG_MGADSM.end() );
             }

GGOptions.insert( GGOptions.end(), G2_GGOptions.begin(), G2_GGOptions.end() );

std::vector<std::string> G3_GGOptions{};


if(useDCT_){ std::vector<std::string> G3_GG_DCT = {"G3G6", "G3G9"};
         G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_DCT.begin(), G3_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G3_GG_DSM = {"G3dG6", "G3dG9"};
         G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_DSM.begin(), G3_GG_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> G3_GG_MGA = {"G3mG6", "G3EG6", "G3mEG6", "G3EmG6", "G3mG9", "G3EG9", "G3mEG9", "G3EmG9"};
         G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_MGA.begin(), G3_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G3_GG_MGADSM = {"G3dmdG6", "G3dEdG6", "G3dmdEdG6", "G3dEdmdG6", "G3dmdG9", "G3dEdG9", "G3dmdEdG9", "G3dEdmdG9"};
         G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_MGADSM.begin(), G3_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G3_GGOptions.begin(), G3_GGOptions.end() );

std::vector<std::string> G4_GGOptions{};


if(useDCT_){ std::vector<std::string> G4_GG_DCT = {"G4G6", "G4G9"};
         G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_DCT.begin(), G4_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G4_GG_DSM = {"G4dG6", "G4dG9"};
         G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_DSM.begin(), G4_GG_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> G4_GG_MGA = {"G4mG6", "G4EG6", "G4mEG6", "G4EmG6", "G4mG9", "G4EG9", "G4mEG9", "G4EmG9"};
         G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_MGA.begin(), G4_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G4_GG_MGADSM = {"G4dmdG6", "G4dEdG6", "G4dmdEdG6", "G4dEdmdG6", "G4dmdG9", "G4dEdG9", "G4dmdEdG9", "G4dEdmdG9"};
         G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_MGADSM.begin(), G4_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G4_GGOptions.begin(), G4_GGOptions.end() );

std::vector<std::string> G5_GGOptions{};


if(useDCT_){ std::vector<std::string> G5_GG_DCT = {"G5G6", "G5G9"};
         G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_DCT.begin(), G5_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G5_GG_DSM = {"G5dG6", "G5dG9"};
         G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_DSM.begin(), G5_GG_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> G5_GG_MGA = {"G5mG6", "G5EG6", "G5mEG6", "G5EmG6", "G5mG9", "G5EG9", "G5mEG9", "G5EmG9"};
         G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_MGA.begin(), G5_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G5_GG_MGADSM = {"G5dmdG6", "G5dEdG6", "G5dmdEdG6", "G5dEdmdG6", "G5dmdG9", "G5dEdG9", "G5dmdEdG9", "G5dEdmdG9"};
         G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_MGADSM.begin(), G5_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G5_GGOptions.begin(), G5_GGOptions.end() );

std::vector<std::string> G6_GGOptions{};


if(useDCT_){ std::vector<std::string> G6_GG_DCT = {"G6G1", "G6G2", "G6G3", "G6G4", "G6G5", "G6G7", "G6G8"};
         G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_DCT.begin(), G6_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G6_GG_DCT = {"G6dG1", "G6dG2", "G6dG3", "G6dG4", "G6dG5", "G6dG7", "G6dG8"};
         G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_DCT.begin(), G6_GG_DCT.end() );
         }
if(useMGA_){ std::vector<std::string> G6_GG_MGA = {"G6MG1", "G6MG2", "G6MG3", "G6MG4", "G6MG5", "G6MG7", "G6EG8", "G6EG1", "G6EG2", "G6EG3", "G6EG4", "G6EG5", "G6EG7", "G6EG8"};
         G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_MGA.begin(), G6_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G6_GG_MGADSM = {"G6dMdG1", "G6dMdG2", "G6dMdG3", "G6dMdG4", "G6dMdG5", "G6dMdG7", "G6dEdG8", "G6dEdG1", "G6dEdG2", "G6dEdG3", "G6dEdG4", "G6dEdG5", "G6dEdG7", "G6dEdG8"};
         G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_MGADSM.begin(), G6_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G6_GGOptions.begin(), G6_GGOptions.end() );

std::vector<std::string> G7_GGOptions{};


if(useDCT_){ std::vector<std::string> G7_GG_DCT = {"G7G6", "G7G9"};
         G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_DCT.begin(), G7_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G7_GG_DSM = {"G7dG6", "G7dG9"};
         G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_DSM.begin(), G7_GG_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> G7_GG_MGA = {"G7mG6", "G7EG6", "G7mEG6", "G7EmG6", "G7mG9", "G7EG9", "G7mEG9", "G7EmG9"};
         G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_MGA.begin(), G7_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G7_GG_MGADSM = {"G7dmdG6", "G7dEdG6", "G7dmdEdG6", "G7dEdmdG6", "G7dmdG9", "G7dEdG9", "G7dmdEdG9", "G7dEdmdG9"};
         G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_MGADSM.begin(), G7_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G7_GGOptions.begin(), G7_GGOptions.end() );

std::vector<std::string> G8_GGOptions{};


if(useDCT_){ std::vector<std::string> G8_GG_DCT = {"G8G6", "G8G9"};
         G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_DCT.begin(), G8_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G8_GG_DSM = {"G8dG6", "G8dG9"};
         G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_DSM.begin(), G8_GG_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> G8_GG_MGA = {"G8mG6", "G8EG6", "G8mEG6", "G8EmG6", "G8mG9", "G8EG9", "G8mEG9", "G8EmG9"};
         G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_MGA.begin(), G8_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G8_GG_MGADSM = {"G8dmdG6", "G8dEdG6", "G8dmdEdG6", "G8dEdmdG6", "G8dmdG9", "G8dEdG9", "G8dmdEdG9", "G8dEdmdG9"};
         G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_MGADSM.begin(), G8_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G8_GGOptions.begin(), G8_GGOptions.end() );

std::vector<std::string> G9_GGOptions{};


if(useDCT_){ std::vector<std::string> G9_GG_DCT = {"G9G1", "G9G2", "G9G3", "G9G4", "G9G5", "G9G7", "G9G8"};
         G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_DCT.begin(), G9_GG_DCT.end() );
         }
if(useDSM_){ std::vector<std::string> G9_GG_DSM = {"G9dG1", "G9dG2", "G9dG3", "G9dG4", "G9dG5", "G9dG7", "G9dG8"};
         G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_DSM.begin(), G9_GG_DSM.end() );
         }
if(useMGA_){ std::vector<std::string> G9_GG_MGA = {"G9MG1", "G9MG2", "G9MG3", "G9MG4", "G9MG5", "G9MG7", "G9EG8", "G9EG1", "G9EG2", "G9EG3", "G9EG4", "G9EG5", "G9EG7", "G9EG8"};
         G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_MGA.begin(), G9_GG_MGA.end() );
         }

if(useMGADSM_){ std::vector<std::string> G9_GG_MGADSM = {"G9dMdG1", "G9dMdG2", "G9dMdG3", "G9dMdG4", "G9dMdG5", "G9dMdG7", "G9dEdG8", "G9dEdG1", "G9dEdG2", "G9dEdG3", "G9dEdG4", "G9dEdG5", "G9dEdG7", "G9dEdG8"};
         G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_MGADSM.begin(), G9_GG_MGADSM.end() );
         }

GGOptions.insert( GGOptions.end(), G9_GGOptions.begin(), G9_GGOptions.end() );

return GGOptions;

    }



    if(Gateway_=="G1"){
        std::vector<std::string> G1_GGOptions{};


    if(useDCT_){ std::vector<std::string> G1_GG_DCT = {"G1G6", "G1G9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_DCT.begin(), G1_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G1_GG_DSM = {"G1dG6", "G1dG9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_DSM.begin(), G1_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G1_GG_MGA = {"G1mG6", "G1EG6", "G1mEG6", "G1EmG6", "G1mG9", "G1EG9", "G1mEG9", "G1EmG9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_MGA.begin(), G1_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G1_GG_MGADSM = {"G1dmdG6", "G1dEdG6", "G1dmdEdG6", "G1dEdmdG6", "G1dmdG9", "G1dEdG9", "G1dmdEdG9", "G1dEdmdG9"};
                 G1_GGOptions.insert( G1_GGOptions.end(), G1_GG_MGADSM.begin(), G1_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G1_GGOptions.begin(), G1_GGOptions.end() );

    }


    if(Gateway_=="G2"){
        std::vector<std::string> G2_GGOptions{};


    if(useDCT_){ std::vector<std::string> G2_GG_DCT = {"G2G6", "G2G9"};
                 G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_DCT.begin(), G2_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G2_GG_DSM = {"G2dG6", "G2dG9"};
                 G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_DSM.begin(), G2_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G2_GG_MGA = {"G2mG6", "G2EG6", "G2mEG6", "G2EmG6", "G2mG9", "G2EG9", "G2mEG9", "G2EmG9"};
                 G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_MGA.begin(), G2_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G2_GG_MGADSM = {"G2dmdG6", "G2dEdG6", "G2dmdEdG6", "G2dEdmdG6", "G2dmdG9", "G2dEdG9", "G2dmdEdG9", "G2dEdmdG9"};
                 G2_GGOptions.insert( G2_GGOptions.end(), G2_GG_MGADSM.begin(), G2_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G2_GGOptions.begin(), G2_GGOptions.end() );

    }

    if(Gateway_=="G3"){
        std::vector<std::string> G3_GGOptions{};


    if(useDCT_){ std::vector<std::string> G3_GG_DCT = {"G3G6", "G3G9"};
                 G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_DCT.begin(), G3_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G3_GG_DSM = {"G3dG6", "G3dG9"};
                 G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_DSM.begin(), G3_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G3_GG_MGA = {"G3mG6", "G3EG6", "G3mEG6", "G3EmG6", "G3mG9", "G3EG9", "G3mEG9", "G3EmG9"};
                 G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_MGA.begin(), G3_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G3_GG_MGADSM = {"G3dmdG6", "G3dEdG6", "G3dmdEdG6", "G3dEdmdG6", "G3dmdG9", "G3dEdG9", "G3dmdEdG9", "G3dEdmdG9"};
                 G3_GGOptions.insert( G3_GGOptions.end(), G3_GG_MGADSM.begin(), G3_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G3_GGOptions.begin(), G3_GGOptions.end() );

    }


    if(Gateway_=="G4"){
        std::vector<std::string> G4_GGOptions{};


    if(useDCT_){ std::vector<std::string> G4_GG_DCT = {"G4G6", "G4G9"};
                 G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_DCT.begin(), G4_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G4_GG_DSM = {"G4dG6", "G4dG9"};
                 G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_DSM.begin(), G4_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G4_GG_MGA = {"G4mG6", "G4EG6", "G4mEG6", "G4EmG6", "G4mG9", "G4EG9", "G4mEG9", "G4EmG9"};
                 G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_MGA.begin(), G4_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G4_GG_MGADSM = {"G4dmdG6", "G4dEdG6", "G4dmdEdG6", "G4dEdmdG6", "G4dmdG9", "G4dEdG9", "G4dmdEdG9", "G4dEdmdG9"};
                 G4_GGOptions.insert( G4_GGOptions.end(), G4_GG_MGADSM.begin(), G4_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G4_GGOptions.begin(), G4_GGOptions.end() );

    }


    if(Gateway_=="G5"){
        std::vector<std::string> G5_GGOptions{};


    if(useDCT_){ std::vector<std::string> G5_GG_DCT = {"G5G6", "G5G9"};
                 G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_DCT.begin(), G5_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G5_GG_DSM = {"G5dG6", "G5dG9"};
                 G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_DSM.begin(), G5_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G5_GG_MGA = {"G5mG6", "G5EG6", "G5mEG6", "G5EmG6", "G5mG9", "G5EG9", "G5mEG9", "G5EmG9"};
                 G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_MGA.begin(), G5_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G5_GG_MGADSM = {"G5dmdG6", "G5dEdG6", "G5dmdEdG6", "G5dEdmdG6", "G5dmdG9", "G5dEdG9", "G5dmdEdG9", "G5dEdmdG9"};
                 G5_GGOptions.insert( G5_GGOptions.end(), G5_GG_MGADSM.begin(), G5_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G5_GGOptions.begin(), G5_GGOptions.end() );

    }


    if(Gateway_=="G7"){
        std::vector<std::string> G7_GGOptions{};


    if(useDCT_){ std::vector<std::string> G7_GG_DCT = {"G7G6", "G7G9"};
                 G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_DCT.begin(), G7_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G7_GG_DSM = {"G7dG6", "G7dG9"};
                 G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_DSM.begin(), G7_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G7_GG_MGA = {"G7mG6", "G7EG6", "G7mEG6", "G7EmG6", "G7mG9", "G7EG9", "G7mEG9", "G7EmG9"};
                 G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_MGA.begin(), G7_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G7_GG_MGADSM = {"G7dmdG6", "G7dEdG6", "G7dmdEdG6", "G7dEdmdG6", "G7dmdG9", "G7dEdG9", "G7dmdEdG9", "G7dEdmdG9"};
                 G7_GGOptions.insert( G7_GGOptions.end(), G7_GG_MGADSM.begin(), G7_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G7_GGOptions.begin(), G7_GGOptions.end() );

    }

    if(Gateway_=="G8"){
        std::vector<std::string> G8_GGOptions{};


    if(useDCT_){ std::vector<std::string> G8_GG_DCT = {"G8G6", "G8G9"};
                 G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_DCT.begin(), G8_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G8_GG_DSM = {"G8dG6", "G8dG9"};
                 G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_DSM.begin(), G8_GG_DSM.end() );
                 }
    if(useMGA_){ std::vector<std::string> G8_GG_MGA = {"G8mG6", "G8EG6", "G8mEG6", "G8EmG6", "G8mG9", "G8EG9", "G8mEG9", "G8EmG9"};
                 G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_MGA.begin(), G8_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G8_GG_MGADSM = {"G8dmdG6", "G8dEdG6", "G8dmdEdG6", "G8dEdmdG6", "G8dmdG9", "G8dEdG9", "G8dmdEdG9", "G8dEdmdG9"};
                 G8_GGOptions.insert( G8_GGOptions.end(), G8_GG_MGADSM.begin(), G8_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G8_GGOptions.begin(), G8_GGOptions.end() );

    }



    if(Gateway_=="G6"){
        std::vector<std::string> G6_GGOptions{};


    if(useDCT_){ std::vector<std::string> G6_GG_DCT = {"G6G1", "G6G2", "G6G3", "G6G4", "G6G5", "G6G7", "G6G8"};
                 G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_DCT.begin(), G6_GG_DCT.end() );
                 }
    if(useDSM_){ std::vector<std::string> G6_GG_DCT = {"G6dG1", "G6dG2", "G6dG3", "G6dG4", "G6dG5", "G6dG7", "G6dG8"};
                 G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_DCT.begin(), G6_GG_DCT.end() );
                 }
    if(useMGA_){ std::vector<std::string> G6_GG_MGA = {"G6MG1", "G6MG2", "G6MG3", "G6MG4", "G6MG5", "G6MG7", "G6EG8", "G6EG1", "G6EG2", "G6EG3", "G6EG4", "G6EG5", "G6EG7", "G6EG8"};
                 G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_MGA.begin(), G6_GG_MGA.end() );
                 }

    if(useMGADSM_){ std::vector<std::string> G6_GG_MGADSM = {"G6dMdG1", "G6dMdG2", "G6dMdG3", "G6dMdG4", "G6dMdG5", "G6dMdG7", "G6dEdG8", "G6dEdG1", "G6dEdG2", "G6dEdG3", "G6dEdG4", "G6dEdG5", "G6dEdG7", "G6dEdG8"};
                 G6_GGOptions.insert( G6_GGOptions.end(), G6_GG_MGADSM.begin(), G6_GG_MGADSM.end() );
                 }

    GGOptions.insert( GGOptions.end(), G6_GGOptions.begin(), G6_GGOptions.end() );

    }

    if(Gateway_=="G9"){
            std::vector<std::string> G9_GGOptions{};


        if(useDCT_){ std::vector<std::string> G9_GG_DCT = {"G9G1", "G9G2", "G9G3", "G9G4", "G9G5", "G9G7", "G9G8"};
                     G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_DCT.begin(), G9_GG_DCT.end() );
                     }
        if(useDSM_){ std::vector<std::string> G9_GG_DSM = {"G9dG1", "G9dG2", "G9dG3", "G9dG4", "G9dG5", "G9dG7", "G9dG8"};
                     G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_DSM.begin(), G9_GG_DSM.end() );
                     }
        if(useMGA_){ std::vector<std::string> G9_GG_MGA = {"G9MG1", "G9MG2", "G9MG3", "G9MG4", "G9MG5", "G9MG7", "G9EG8", "G9EG1", "G9EG2", "G9EG3", "G9EG4", "G9EG5", "G9EG7", "G9EG8"};
                     G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_MGA.begin(), G9_GG_MGA.end() );
                     }

        if(useMGADSM_){ std::vector<std::string> G9_GG_MGADSM = {"G9dMdG1", "G9dMdG2", "G9dMdG3", "G9dMdG4", "G9dMdG5", "G9dMdG7", "G9dEdG8", "G9dEdG1", "G9dEdG2", "G9dEdG3", "G9dEdG4", "G9dEdG5", "G9dEdG7", "G9dEdG8"};
                     G9_GGOptions.insert( G9_GGOptions.end(), G9_GG_MGADSM.begin(), G9_GG_MGADSM.end() );
                     }

        GGOptions.insert( GGOptions.end(), G9_GGOptions.begin(), G9_GGOptions.end() );

        }


    return GGOptions;



}









