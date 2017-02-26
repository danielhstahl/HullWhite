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
TEST_CASE("Test forward swap", "[HullWhite]"){
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
        hullwhite::Swap_Rate(
            currRate,
            a, 
            sig,
            futureTime,
            swpMaturity,
            delta,
            yield,
            forward 
        )==Approx(
        hullwhite::Forward_Swap_Rate(
            currRate, a, sig, 
            futureTime, futureTime, swpMaturity, delta,
            yield, forward
        ))
    );
}
TEST_CASE("Payer Swaption", "[HullWhite]"){
    const double currRate=.05;
    const double sig=.01;
    const double a=.05;
    const double b=.05;
    const double delta=.25;
    const double futureTime=0;
    const double swpTenor=5;
    const double optMaturity=1;
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
    const double swaprate=hullwhite::Forward_Swap_Rate(
        currRate,
        a, 
        sig,
        futureTime,
        optMaturity,
        swpTenor+optMaturity,
        delta,
        yield,
        forward 
    );    
    auto analytic=hullwhite::Payer_Swaption(currRate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
    auto tree=hullwhite::PayerSwaptionWithTree(currRate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
    REQUIRE(analytic==Approx(tree).epsilon(.0001));
}
TEST_CASE("Receiver Swaption", "[HullWhite]"){
    const double currRate=.05;
    const double sig=.01;
    const double a=.05;
    const double b=.05;
    const double delta=.25;
    const double futureTime=0;
    const double swpTenor=5;
    const double optMaturity=1;
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
    const double swaprate=hullwhite::Forward_Swap_Rate(
        currRate,
        a, 
        sig,
        futureTime,
        optMaturity,
        swpTenor+optMaturity,
        delta,
        yield,
        forward 
    );    
    auto analytic=hullwhite::Receiver_Swaption(currRate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
    auto tree=hullwhite::ReceiverSwaptionWithTree(currRate, a, sig, swaprate, futureTime, swpTenor, optMaturity, delta, yield, forward);
    REQUIRE(analytic==Approx(tree).epsilon(.0001));
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
TEST_CASE("ZCB to CB option", "[HullWhite]"){
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
    std::vector<double> fakeCoupon={bondMaturity};
    REQUIRE(
        hullwhite::Bond_Call(currRate, a, sig, futureTime, optMaturity, bondMaturity,strike, yield, forward)==
        hullwhite::Coupon_Bond_Call(currRate, a, sig, strike, futureTime, optMaturity, fakeCoupon,0.0, yield, forward)
    );
}
