#ifndef __YIELDIO_H_INCLUDED__
#define __YIELDIO_H_INCLUDED__
#include <string>
#include <iostream>
#include "RealWorldMeasure.h"
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson


std::vector<SpotValue> populateYieldFromExternalSource(Date&, auto&, auto&);
#include "YieldIO.hpp"

#endif