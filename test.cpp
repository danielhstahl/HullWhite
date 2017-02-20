#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include "HullWhite.h"
#include "AutoDiff.h"


TEST_CASE("Test Bond Price", "[HullWhite]"){
    const double currRate=.02;
    const double sig=.02;
    const double a=.3;
    const double b=.04;
    const double delta=.25;
    const double strike=.04;
    const double futureTime=.5;
    const double swpMaturity=5.5;
    const double optMaturity=1.5;
    auto square=[](auto&& val){
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
    REQUIRE(hullwhite::Bond_Price(currRate, a, sig, futureTime, optMaturity, yield, forward)==hullwhite::Bond_Price(optMaturity-futureTime, yield));
}
TEST_CASE("Test swap", "[HullWhite]"){
    const double currRate=.02;
    const double sig=.02;
    const double a=.3;
    const double b=.04;
    const double delta=.25;
    const double strike=.04;
    const double futureTime=.5;
    const double swpMaturity=5.5;
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
    REQUIRE(
        hullwhite::Swap_Price(
            currRate, a, sig, 
            futureTime, swpMaturity, delta,
            hullwhite::Swap_Rate(
               currRate,
               a, 
               sig,
               futureTime,
               swpMaturity,
               delta,
               yield,
               forward 
            ),
            yield, forward
        )==Approx(0.0)
    );
}
TEST_CASE("Swaption", "[HullWhite]"){
    const double currRate=.02;
    const double sig=.02;
    const double a=.3;
    const double b=.04;
    const double delta=.25;
    const double swaprate=.04;
    const double futureTime=.5;
    const double swpTenor=5;
    const double optMaturity=1.5;
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
    //std::cout<<"American "<<hullwhite::AmericanSwaption(currRate, a, sig, swaprate, 0.0, swpTenor, optMaturity, delta, yield, forward)<<std::endl;
    //auto analytic=hullwhite::Payer_Swaption(currRate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
    for(int i=0; i<5; ++i){
        auto rate=.000001+i*.1;
        auto analytic=hullwhite::Payer_Swaption(rate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
        auto tree=hullwhite::SwaptionWithTree(rate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
        std::cout<<"rate "<<rate<<", value "<<analytic<<", tree "<<tree<<", ratio "<<tree/analytic<<std::endl;
    }
    /*REQUIRE(
        hullwhite::PayerSwaption(currRate, a, sig, swaprate, 0.0, swpTenor, optMaturity, delta, yield, forward)==
        Approx(
            hullwhite::SwaptionWithTree(currRate, a, sig, swaprate, 0.0, swpTenor, optMaturity, delta, yield, forward)
        )
    );*/
}
TEST_CASE("ZCB option Reference", "[HullWhite]"){
    //http://www.quantcalc.net/BondOption_Vasicek.html
    const double currRate=.01;
    const double sig=.03;
    const double a=.05;
    const double b=.04;
    const double delta=.25;
    const double strike=.96;
    const double futureTime=0;
    const double bondMaturity=3;
    const double optMaturity=2;
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
    REQUIRE(
        hullwhite::Bond_Call(currRate, a, sig, futureTime, optMaturity, bondMaturity,strike, yield, forward)==
        Approx(
            .033283
        )
    );
}
