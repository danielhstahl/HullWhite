#ifndef __REALWORLDMEASURE_H_INCLUDED__
#define __REALWORLDMEASURE_H_INCLUDED__
#include <cmath>
#include "CurveFeatures.h"

double findHistoricalMean(std::vector<SpotValue>&);
auto generateVasicek(
    const auto&, //current value
    const auto&, //time to next simulation
    const auto&, //a
    const auto&, //b
    const auto&,//sigma
    const auto& //normal simulation
); 
#include "RealWorldMeasure.hpp"

#endif