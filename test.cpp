#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include "HullWhite.h"
#include "AutoDiff.h"

const double currRate=.02;
const double sig=.02;
const double a=.3;
const double b=.04;
const double delta=.25;
const double strike=.04;
const double futureTime=.5;
const double swpMaturity=5.5;
const double optMaturity=1.5;
/*auto bndV=[&](double r, double a, double b, double sigma, double t){
    double at=(1-exp(-a*t))/a;
    double ct=sigma*sigma;
    ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
    return exp(-at*r+ct);
};*/
auto square=[](const auto& val){
    return val*val;
};
auto yield=[&](const auto& t){
    auto at=(1-exp(-a*t))/a;
    auto ct=(b-sig*sig/(2*a*a))*(at-t)-sig*sig*at*at/(4*a);
    return (at*currRate-ct);
};
auto forward=[&](const auto& t){
    return b+exp(-a*t)*(currRate-b)-(sig*sig/(2*a*a))*square(1-exp(-a*t));
};

std::vector<double> couponTimes(5);

TEST_CASE("Test stuff", "[HullWhite]"){
    couponTimes[0]=.5;
    couponTimes[1]=1;
    couponTimes[2]=1.5;
    couponTimes[3]=2;
    couponTimes[4]=2.5;
    //std::cout<<yield(AutoDiff<double>(5.0, 1.0)).getDual()<<std::endl;
    /*const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& strike,
    const FirstFutureTime& t, 
    const SecondFutureTime& T, 
    const ThirdFutureTime& TM, 
    const Delta& delta, 
    const GetYield& yield,
    const GetInstantaneousForward& forward*/
    REQUIRE(hullwhite::Bond_Price(currRate, a, sig, futureTime, optMaturity, yield, forward)==hullwhite::Bond_Price(optMaturity-futureTime, yield));
    //std::cout<<hullwhite::Bond_Price(currRate, a, sig, futureTime, optMaturity, yield, forward)<<std::endl;
    //std::cout<<hullwhite::Bond_Price(optMaturity-futureTime, yield)<<std::endl;

    std::cout<<hullwhite::Swaption(currRate, a, sig, strike, 0.0, 4.5, .5, delta, yield, forward)<<std::endl;
    //std::cout<<AmericanSwaption(currRate, a, sig, strike, futureTime, swpMaturity, optMaturity, delta, yld)<<std::endl;
    /*AutoDiff<double> curr(currRate, 1.0);
    std::cout<<hullwhite::AmericanSwaption(currRate, a, sig, strike, 0.0, 4.5, .5, delta, yield, forward)<<std::endl;

    auto ev=hullwhite::Swaption(curr, a, sig, strike, 0.0, 4.5, .5, delta, yield, forward);

    std::cout<<ev.getStandard()<<std::endl;
    std::cout<<ev.getDual()<<std::endl;


    auto v=hullwhite::AmericanSwaption(curr, a, sig, strike, 0.0, 4.5, .5, delta, yield, forward);
    std::cout<<v.getStandard()<<std::endl;
    std::cout<<v.getDual()<<std::endl;*/
    //REQUIRE(futilities::for_each_parallel(std::move(testV), squareTestV)==std::vector<int>({25, 36, 49}));
}
