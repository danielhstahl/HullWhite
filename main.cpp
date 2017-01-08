#define _USE_MATH_DEFINES
#include <iostream>
#include "BlackScholes.h"
#include <cfloat>
#include "HullWhiteEngine.h"
#include <string>
#include "YieldIO.h"
#include "RealWorldMeasure.h"
#include "CurveFeatures.h"
#include "YieldSpline.h"
#include <unordered_map>
#include "HandlePath.h" //for creating sample paths
#include "MC.h" //monte carlo
#include "Histogram.h" //bins data
#include "SimulateNorm.h"
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson

 int main(){

  /* int n=29;//number of test yields
    std::vector<SpotValue> testYield;
    double currRate=.02;
    double sig=.02;
    double a=.3;
    double b=.04;

    auto bndV=[&](double r, double a, double b, double sigma, double t){
        double at=(1-exp(-a*t))/a;
        double ct=sigma*sigma;
        ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
        return exp(-at*r+ct);
    };


    Date currDate;
    for(int i=0; i<n; ++i){
        testYield.push_back(SpotValue(currDate+(i+1), (1/bndV(currRate, a, b, sig, (currDate+(i+1))-currDate)-1)/((currDate+(i+1))-currDate)));
      //std::cout<<"price: "<<bndV(currRate, a, b, sig, (currDate+(i+1))-currDate)<<std::endl;
  }
  YieldSpline yld(testYield, currDate, currRate);
    //std::cout<<exp(-yld.Yield(2))<<std::endl;
  std::vector<double> couponTimes(5);
  couponTimes[0]=.5;
  couponTimes[1]=1;
  couponTimes[2]=1.5;
  couponTimes[3]=2;
  couponTimes[4]=2.5;
    //std::cout<<bndV(currRate, a, b, sig, 1)<<std::endl;
    //std::cout<<exp(-yld.Yield(1))<<std::endl;
    //std::cout<<yld.Forward(0)<<std::endl;
    //std::cout<<Bond_Price(currRate, a, sig, 0, 1, yld)<<std::endl;
    //double bndt=Bond_Price(1, yld);
    //std::cout<<"Price: "<<Bond_Price(1.0, yld)<<std::endl;

    //double sig=.02;
    double strike=.04;
    double futureTime=.5;
    double swpMaturity=5.5;
    double optMaturity=1.5;
    double delta=.25;
     AutoDiff curr(currRate, 1.0);
    //std::cout<<Swaption(currRate, a, sig, strike, futureTime, swpMaturity, optMaturity, delta, yld)<<std::endl;
    std::cout<<Swaption(currRate, a, sig, strike, 0.0, 4.5, .5, delta, yld)<<std::endl;
    //std::cout<<AmericanSwaption(currRate, a, sig, strike, futureTime, swpMaturity, optMaturity, delta, yld)<<std::endl;

     std::cout<<AmericanSwaption(currRate, a, sig, strike, 0.0, 4.5, .5, delta, yld)<<std::endl;

     AutoDiff ev=Swaption(curr, a, sig, strike, 0.0, 4.5, .5, delta, yld);

    std::cout<<ev.getStandard()<<std::endl;
    std::cout<<ev.getDual()<<std::endl;


    AutoDiff v=AmericanSwaption(curr, a, sig, strike, 0.0, 4.5, .5, delta, yld);
    std::cout<<v.getStandard()<<std::endl;
    std::cout<<v.getDual()<<std::endl;
    */


    Date currDate;   
    YieldSpline yld;
    double b;//long run average
    double daysDiff;//years from now that libor rate goes till (typically 7 days divided by 360)
    std::vector<SpotValue> historical=populateYieldFromExternalSource(currDate, yld, daysDiff);//this will wait from input from external source
    yld.getSpotCurve();//send data to node
    yld.getForwardCurve(); //send data ot node
    HullWhiteEngine<double> HW;
    double r0=yld.getShortRate(); //note we can change this here to an AutoDiff if we want sensitivities
    SimulateNorm rNorm;
    MC<double> monteC;
    auto runParameters=[&](std::string& parameters){
        rapidjson::Document parms;
        parms.Parse(parameters.c_str());//yield data
        parameters.clear();
        std::vector<AssetFeatures> portfolio;
        AssetFeatures asset;
        //parse the string that came in
        currDate.setScale("year");
        if(parms.FindMember("T")!=parms.MemberEnd()){
            asset.Maturity=currDate+parms["T"].GetDouble();
        }
        if(parms.FindMember("k")!=parms.MemberEnd()){
            asset.Strike=parms["k"].GetDouble();
        }
        if(parms.FindMember("delta")!=parms.MemberEnd()){
            asset.Tenor=parms["delta"].GetDouble();
        }
        if(parms.FindMember("Tm")!=parms.MemberEnd()){
            asset.UnderlyingMaturity=currDate+parms["Tm"].GetDouble();
        }
        asset.type=parms["asset"].GetInt();
        double a=parms["a"].GetDouble(); //can be made autodiff too
        double sigma=parms["sigma"].GetDouble(); //can be made autodiff too
        HW.setSigma(sigma);
        HW.setReversion(a);
        b=findHistoricalMean(historical, daysDiff, a);

        currDate.setScale("day");
        Date PortfolioMaturity;
        if(parms.FindMember("t")!=parms.MemberEnd()){
            PortfolioMaturity=currDate+parms["t"].GetInt();
        }
        int m=0;
        if(parms.FindMember("n")!=parms.MemberEnd()){
            m=parms["n"].GetInt();
        }
        monteC.setM(m);
        portfolio.push_back(asset);
        std::vector<Date> path=getUniquePath(portfolio, PortfolioMaturity);

        monteC.simulateDistribution([&](){
            return executePortfolio(portfolio, currDate,
                [&](const auto& currVal, const auto& time){
                    double vl=rNorm.getNorm();
                    return generateVasicek(currVal, time, a, b, sigma, vl);
                },
                r0,
                path,
                [&](AssetFeatures& asset, auto& rate, Date& maturity,   Date& asOfDate){
                    return HW.HullWhitePrice(asset, rate, maturity, asOfDate, yld);
                }
            );
        });

        std::vector<double> dist=monteC.getDistribution();
        double min=DBL_MAX; //purposely out of order because actual min and max are found within the function
        double max=DBL_MIN;
        binAndSend(min, max, dist); //send histogram to node
    };
    while(true){
        std::string parameters;
        for (parameters; std::getline(std::cin, parameters);) {
            break;
        }
        runParameters(parameters);
    }



    /*
  std::unordered_map<std::string, AutoDiff> parameters;
  parameters.insert({"Underlying", AutoDiff(50, 0)});
  parameters.insert({"Strike", AutoDiff(50, 0)});
  parameters.insert({"Maturity", AutoDiff(1.0, 0)});
  parameters.insert({"Sigma", AutoDiff(0.3, 0)});
  parameters.insert({"R", AutoDiff(.03, 0)});

  bool is_first_iteration = true;
  std::cout<<"The following are the parameters used in this demonstration:"<<std::endl;
  for(const auto & parameterPair : parameters) {
     std::cout << parameterPair.first << ": ";
     std::cout<<parameterPair.second.getStandard();
     std::cout << std::endl;
  }
  std::cout<<"Choose one of the parameters to find the price and derivative with respect to that parameter: ";
  std::string response;
  std::cin>>response;
 // std::cout<<response<<std::endl;
  while(parameters.find(response)==parameters.end()){
    std::cout<<"Parameter doesn't exist!  Choose a valid parameter:"<<std::endl;
    std::cin>>response;
  }
  parameters.find(response)->second.setDual(1);
  AutoDiff discount=exp(-parameters.find("R")->second*parameters.find("Maturity")->second);
  AutoDiff Call=BSCall(parameters.find("Underlying")->second, discount, parameters.find("Strike")->second, parameters.find("Sigma")->second*sqrt(parameters.find("Maturity")->second));
  AutoDiff Put=BSPut(parameters.find("Underlying")->second, discount, parameters.find("Strike")->second, parameters.find("Sigma")->second*sqrt(parameters.find("Maturity")->second));
  std::cout<<"Call Price: "<<Call.getStandard()<<" Partial Derivative with respect to "<<response<<": "<<Call.getDual()<<std::endl;
  std::cout<<"Put Price: "<<Put.getStandard()<<" Partial Derivative with respect to "<<response<<": "<<Put.getDual()<<std::endl;*/
}
