#ifndef __HULLWHITEENGINE_H_INCLUDED__
#define __HULLWHITEENGINE_H_INCLUDED__
#include <cmath>
#include "HullWhite.h"
#include "Date.h"
#include "CurveFeatures.h" 
#include "YieldSpline.h"

#define BOND 0
#define EURODOLLARFUTURE 1
#define COUPONBOND 2
#define BONDCALL 3
#define BONDPUT 4
#define COUPONBONDCALL 5
#define COUPONBONDPUT 6
#define CAPLET 7
#define SWAP 8 
#define SWAPTION 9
#define AMERICANSWAPTION 10

 /*struct AssetFeatures{
    Date Maturity;
    Date UnderlyingMaturity;
    double Strike;
    double Tenor;
    const int type; //types are defined where?
     
 };*/


template<typename Number> //autodiff or double
class HullWhiteEngine{
    private:
        //YieldSpline *yld;//pointer to yield class
        Number a;
        Number sigma;
        Number r;
    public:
        HullWhiteEngine();    
        //HullWhiteEngine(Number&, Number&);  
        void setSigma(Number&);
        void setReversion(Number&);
        //void setYield(auto*);
        //void setShortRate(auto&);
        auto HullWhitePrice(
            AssetFeatures&, 
            auto&, //rate
            Date&, //future time,
            Date&, //asOfDate
            YieldSpline& //yield class
        );
    
};
#include "HullWhiteEngine.hpp"





#endif