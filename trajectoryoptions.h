#ifndef TRAJECTORYOPTIONS_H
#define TRAJECTORYOPTIONS_H

#include<vector>
#include<string>
#include<iostream>

class TrajectoryOptions
{
public:
    TrajectoryOptions(std::string Gateway = "All", std::string Segment = "All", bool useDCT=true, bool useDSM=true, bool useMGA=true, bool useMGADSM=true);

    std::vector<std::string> makeStringVector();

    std::vector<std::string> makeG1options();
    std::vector<std::string> makeG2options();
    std::vector<std::string> makeG3options();
    std::vector<std::string> makeG4options();
    std::vector<std::string> makeG5options();
    std::vector<std::string> makeG6options();
    std::vector<std::string> makeG7options();
    std::vector<std::string> makeG8options();
    std::vector<std::string> makeG9options();
    std::vector<std::string> makeGGoptions();





private:
    std::string Gateway_;
    std::string Segment_;
    bool useDCT_;
    bool useDSM_;
    bool useMGA_;
    bool useMGADSM_;
    bool returnFlight_;


};
#endif // TRAJECTORYOPTIONS_H
