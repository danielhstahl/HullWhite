std::vector<SpotValue> populateYieldFromExternalSource(Date& dt, auto& yldClass, auto& daysPlus){
    std::string yieldData;
    std::string historicalData;

    for (yieldData; std::getline(std::cin, yieldData);) {
        break;
    }
    for (historicalData; std::getline(std::cin, historicalData);) {
        break;
    }

    rapidjson::Document dYield;
    rapidjson::Document dHistorical;
 
    dYield.Parse(yieldData.c_str());//yield data
    dHistorical.Parse(historicalData.c_str()); //historical data

    yieldData.clear();
    historicalData.clear();

    int n=dYield.Size();
    int m=dHistorical["data"]["observations"].Size();
  //Deal with yield
    std::vector<SpotValue> SwapRates;
    std::vector<SpotValue> LiborRates;
    dt.setScale("day");
    for(int i=0; i<n; ++i){
        if(strcmp(dYield[i]["type"].GetString(), "Swap")==0){
            SwapRates.push_back(SpotValue(dt+dYield[i]["daysPlus"].GetInt(), std::stod(dYield[i]["value"].GetString())*.01));
        }
        else{
            LiborRates.push_back(SpotValue(dt+dYield[i]["daysPlus"].GetInt(), std::stod(dYield[i]["value"].GetString())*.01));
        }
    }

    //end deal with yield
    //Deal with historical
    std::vector<SpotValue> historical;
    std::string val;
    for(int i=0; i<m; ++i){
        val=dHistorical["data"]["observations"][i]["value"].GetString();
        
        if(val!="."&&val!=""){//nulls show up as "."
            historical.push_back(SpotValue(dHistorical["data"]["observations"][i]["date"].GetString(), std::stod(val)*.01));
        }
    }
    yldClass.computeSimpleSwapSpline(LiborRates, SwapRates, dt);
    daysPlus=dHistorical["daysPlus"].GetInt()/360.0;
    return historical;
    
    //b=findHistoricalMean(historical, dHistorical["daysPlus"].GetInt()/360.0, a);
}
    