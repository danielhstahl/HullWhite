#ifndef __HULLWHITE_H_INCLUDED__
#define __HULLWHITE_H_INCLUDED__
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <vector>
#include <unordered_map>
#include "BTree.h"
#include "BlackScholes.h" //for option pricing
#include "Newton.h"
#include "FunctionalUtilities.h"

constexpr double prec1=.0000001;
constexpr double prec2=.0000001;
constexpr double r_init=.03;
constexpr int max_iter=50;
/*Note: the fundamental times here are (0, t, T, TM).  0 is current time (and is reflective of the current yield curve), t is some future time that we may want to price options at given the underlying at that time, T is an "initial" maturity and TM a "Final" maturity.  While it is natural to think of (0<t<T<TM), I only require 0<t and 0<T<TM. Note that ALL TIMES ARE WITH RESPECT TO 0!  */
namespace hullwhite{
  template<typename MeanRevertSpeed, typename DiffFutureTimes>
  auto A( /*A(T-t) from the Hull White PDE*/
    const MeanRevertSpeed& a,/*a*/
    const DiffFutureTimes& tDiff/*T-t*/
  ){
    return (1.0-exp(-a*tDiff))/a;
  }
  template<typename MeanRevertSpeed, typename FirstFutureTime, typename SecondFutureTime>
  auto A( /*A(t, T) from the Hull White PDE*/
    const MeanRevertSpeed& a,/*a*/
    const FirstFutureTime& t,/*t*/
    const SecondFutureTime& T /*T */
  ){
    return A(a, T-t);
  }
  
  template<typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename GetYield, typename GetInstantaneousForward>
  auto C( /*C(t, T) from the Hull White PDE...this may require taking a reference to a yield class!*/
    const MeanRevertSpeed& a, /*a */
    const Volatility& sigma, /*a */
    const FirstFutureTime& t,/*t*/
    const SecondFutureTime& T, /*T */
    const GetYield& yield, /*does yield*/
    const GetInstantaneousForward& forward /* does instantaenous froward*/
  ){
    auto sqr=exp(-a*T)-exp(-a*t);
    return yield(t)-yield(T)+forward(t)*A(a, t, T)-sigma*sigma*sqr*sqr*(exp(2.0*a*t)-1.0)/(4.0*a*a*a);
  }
  template<typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime>
  auto T_Forward_Bond_Volatility( /*This is the volatility of a forward bond B(t, TM)/(B(t, T)), TM>T and can be used for zero coupon bond pricing using BS formula.  TM may equal T+delta, but must be manually added before calling this function*/
    const MeanRevertSpeed& a,/*Speed of mean reversion ("a")*/
    const Volatility& sigma, /*interest rate volatility ("sigma") */
    const FirstFutureTime& t, /*future time t */
    const SecondFutureTime& T, /*Initial Bond maturity T*/
    const ThirdFutureTime& TM/*Final Bond maturity TM*/
  ){
      auto expD=1-exp(-a*(TM-T));
      auto expT=1-exp(-2*a*(T-t));
      return sigma*sqrt(expT/(2.0*a*a*a))*expD;
  }

  template<typename MeanRevertSpeed, typename Volatility,typename SecondFutureTime, typename GetInstantaneousForward>
  auto phiT(
      const MeanRevertSpeed& a, /*speed of mean reversion*/
      const Volatility& sigma, /*interst rate volatility*/
      const SecondFutureTime& T, /*expectation horizon*/
      const GetInstantaneousForward& forward/*functions "Forward" should be the instantanoues forward rate */
  ){
      auto expT=1-exp(-a*T);
      return forward(T)+sigma*sigma*expT*expT/(2*a*a);
      
  }

  template<typename R, typename MeanRevertSpeed, typename Volatility,typename FirstFutureTime, typename SecondFutureTime, typename GetInstantaneousForward>
  auto muR(/*the expected value of r_t under risk neutral measure: E[r_T|t]*/
    const R& r_t,
    const MeanRevertSpeed& a, /*speed of mean reversion*/
    const Volatility& sigma, /*interst rate volatility*/
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*expectation horizon*/
    const GetInstantaneousForward& forward/*functions "Forward" should be the instantanoues forward rate */
  ){
      return phiT(a, sigma, T, forward)+(r_t-phiT(a, sigma, t, forward))*exp(-a*(T-t));
  }
  template<typename MeanRevertSpeed, typename Volatility,typename FirstFutureTime, typename SecondFutureTime>
  auto varianceR(/*the variance of r_t under risk neutral measure */
    const MeanRevertSpeed& a, /*speed of mean reversion*/
    const Volatility& sigma, /*interst rate volatility*/
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T /*expectation horizon*/
  ){
    return sigma*sigma*(1.0-exp(-2.0*a*(T-t)))/(2.0*a);
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename GetYield, typename GetInstantaneousForward>
  auto Bond_Price(/*The zero coupon bond price under Hull White*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*Bond expiration*/
    const GetYield& yield,/*continuously compounded zero coupon yield...multiplied by T*/
    const GetInstantaneousForward& forward //instantanoues forward rate 
  ){
    return exp(-r_t*A(a, t, T)+C(a, sigma, t, T, yield, forward));
  }
  template<typename SecondFutureTime,typename GetYield>  
  auto Bond_Price(/*The zero coupon bond price under Hull White for t=0*/
    const SecondFutureTime& T, /*Bond expiration*/
    const GetYield& yield/*continuously compounded zero coupon (yield*t) at time 0*/
  ){
    return exp(-yield(T));
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename Coupon, typename GetYield, typename GetInstantaneousForward>
  auto Coupon_Bond_Price(/*The coupon bond price under Hull White*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    std::vector<SecondFutureTime>& couponTimes, /*these are coupon times FROM 0!  these should be in order but dont have to be*/
    const Coupon& couponRate,
    const GetYield& yield,/*the continuously compounded zero coupon yield*/
    const GetInstantaneousForward& forward //instantanoues forward rate
  ){
    std::sort(couponTimes.begin(), couponTimes.end());
    return futilities::sum(couponTimes, [&](const auto& val, const auto& index){
      return val>t?Bond_Price(r_t, a, sigma, t, val, yield, forward)*couponRate:0.0;
    })+Bond_Price(r_t, a, sigma, t, couponTimes.back(), yield, forward);
  }
  template<typename SecondFutureTime, typename Coupon, typename GetYield>
  auto Coupon_Bond_Price(/*The coupon bond price under Hull White at time 0*/
    std::vector<SecondFutureTime>& couponTimes, /*these are coupon times FROM 0!  these should be in order*/
    const Coupon& couponRate,
    const GetYield& yield/*continuously compounded zero coupon yield*/
  ){
    std::sort(couponTimes.begin(), couponTimes.end());
    return futilities::sum(couponTimes, [&](const auto& val, const auto& index){
      return Bond_Price(val, yield)*couponRate;
    })+Bond_Price(couponTimes.back(), yield);
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime,  typename Strike, typename GetYield, typename GetInstantaneousForward>
  auto Bond_Call(/*The price of a call option on zero coupon bond under Hull White*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*option maturity*/
    const ThirdFutureTime& TM, /*bond maturity*/
    const Strike& strike,
    const GetYield& yield, //continuously compounded zero coupon yield
    const GetInstantaneousForward& forward //instantanoues forward rate 
  ){
    return BSCall(
      Bond_Price(r_t, a, sigma, t, TM, yield, forward), /*underlying*/
      Bond_Price(r_t, a, sigma, t, T, yield, forward), /*discount factor*/
      strike,
      T_Forward_Bond_Volatility(a, sigma, t, T, TM)/*volatility of underlying*/
    );
  }
  template<typename MeanRevertSpeed, typename Volatility, typename SecondFutureTime, typename ThirdFutureTime, typename Strike, typename GetYield>
  auto Bond_Call(/*The price of a call option on zero coupon bond under Hull White at t=0*/
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const SecondFutureTime& T, /*option maturity*/
    const ThirdFutureTime& TM, /*bond maturity*/
    const Strike& strike,
    const GetYield& yield/*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
  ){
    return BSCall(
      Bond_Price(TM, yield), /*underlying*/
      Bond_Price(T, yield), /*discount factor*/
      strike,
      T_Forward_Bond_Volatility(a, sigma, 0, T, TM) /*volatility of underlying*/
    );
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Coupon, typename Strike, typename GetYield, typename GetInstantaneousForward>
  auto Coupon_Bond_Call(/*The price of a call option on coupon bond under Hull White...uses jamshidian's trick*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& strike,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*option maturity*/
    std::vector<ThirdFutureTime>& couponTimes, /*these are coupon times FROM 0!  these should be in order but dont have to be*/
    const Coupon& couponRate,
    const GetYield& yield, //continuously compounded zero coupon yield
    const GetInstantaneousForward& forward //instantanoues forward rate
  ){
    std::sort(couponTimes.begin(), couponTimes.end());
    auto myOptimalR=newton::zeros([&](auto &r){
      return Coupon_Bond_Price(r, a, sigma, T, couponTimes, couponRate, yield, forward)-strike; //T is "future" time since bond is priced at opion maturity 
    }, r_init, prec1, prec2, max_iter);

    return futilities::sum(couponTimes, [&](const auto& val, const auto& index){
      return val>T?couponRate*Bond_Call(r_t, a, sigma, t, T, val, Bond_Price(myOptimalR, a, sigma, T, val, yield, forward), yield, forward):0.0;
    })+(couponTimes.back()>T?Bond_Call(r_t, a, sigma, t, T, couponTimes.back(), Bond_Price(myOptimalR, a, sigma, T, couponTimes.back(), yield, forward), yield, forward):0.0);
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Strike, typename GetYield, typename GetInstantaneousForward>
  auto Bond_Put(/*The price of a Put option on zero coupon bond under Hull White*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*option maturity*/
    const ThirdFutureTime& TM, /*bond maturity*/
    const Strike& strike,
    const GetYield& yield,/*the continuously compounded zero coupon yield*/
    const GetInstantaneousForward& forward //instantanoues forward rate and 
  ){
    return BSPut(
      Bond_Price(r_t, a, sigma, t, TM, yield, forward), /*underlying*/
      Bond_Price(r_t, a, sigma, t, T, yield, forward), /*discount factor*/
      strike,
      T_Forward_Bond_Volatility(a, sigma, t, T, TM) /*volatility of underlying*/
    );
  }
  template<typename MeanRevertSpeed, typename Volatility, typename SecondFutureTime, typename ThirdFutureTime, typename Strike, typename GetYield>
  auto Bond_Put(/*The price of a Put option on zero coupon bond under Hull White at t=0*/
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const SecondFutureTime& T, /*option maturity*/
    const ThirdFutureTime& TM, /*bond maturity*/
    const Strike& strike,
    const GetYield& yield/*continuously compounded zero coupon yield*/
  ){
    return BSPut(
      Bond_Price(TM, yield), /*underlying*/
      Bond_Price(T, yield), /*discount factor*/
      strike,
      T_Forward_Bond_Volatility(a, sigma, 0, T, TM) /*volatility of underlying*/
    );
  }


  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Coupon, typename Strike, typename GetYield, typename GetInstantaneousForward>
  auto Coupon_Bond_Put(/*The price of a call option on coupon bond under Hull White...uses jamshidian's trick*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& strike,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*option maturity*/
    std::vector<ThirdFutureTime>& couponTimes, /*these are coupon times FROM 0!  these should be in order but dont have to be*/
    const Coupon& couponRate,
    const GetYield& yield, //continuously compounded zero coupon yield
    const GetInstantaneousForward& forward //instantanoues forward rate
  ){
    std::sort(couponTimes.begin(), couponTimes.end());
    auto myOptimalR=newton::zeros([&](auto &r){
      return Coupon_Bond_Price(r, a, sigma, T, couponTimes, couponRate, yield, forward)-strike; //T is "future" time since bond is priced at opion maturity
    }, r_init, prec1, prec2, max_iter);
    return (futilities::sum(couponTimes, [&](const auto& val, const auto& index){
      return val>T?Bond_Put(r_t, a, sigma, t, T, val, Bond_Price(myOptimalR, a, sigma, T, val, yield, forward), yield, forward):0.0;
    })*couponRate)+(couponTimes.back()>T?Bond_Put(r_t, a, sigma, t, T, couponTimes.back(), Bond_Price(myOptimalR, a, sigma, T, couponTimes.back(), yield, forward), yield, forward):0.0);
  }

  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename Delta,  typename Strike, typename GetYield, typename GetInstantaneousForward>
  auto Caplet(/*price of a caplet on simple bond yield*/
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*option maturity*/
    const Delta& delta, /*tenor of the simple yield*/
    const Strike& strike,
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    return (strike*delta+1.0)*Bond_Put(r_t, a, sigma, t, T, T+delta, 1.0/(delta*strike+1.0), yield, forward);
  }

  template<typename MeanRevertSpeed, typename Volatility,  typename SecondFutureTime, typename Delta,  typename Strike, typename GetYield, typename GetInstantaneousForward>
  auto Caplet(/*price of a caplet on simple bond yield at time t=0*/
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const SecondFutureTime& T, /*option maturity*/
    const Delta& delta, /*tenor of the simple yield*/
    const Strike& strike,
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    return (strike*delta+1)*Bond_Put(a, sigma, T, T+delta, 1/(delta*strike+1), yield, forward);
  }

  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto Euro_Dollar_Future(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*maturity*/
    const Delta& delta, /*tenor of the simple yield*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
      auto expT=exp(-a*(T-t));
      auto expD=exp(-a*delta);
      auto gamma=(sigma*sigma/(a*a*a))*(1.0-expD)*((1.0-expT)-expD*.5*(1.0-expT*expT));
      return ((Bond_Price(r_t, a, sigma, t, T, yield, forward)/Bond_Price(r_t, a, sigma, t, T+delta, yield, forward))*exp(gamma)-1.0)/delta;
  }

  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto Swap_Rate(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*swap maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    int numPayments=floor((T-t)/delta)+1; //this should be an integer!  remember, T-t is the total swap length
    return (1.0-Bond_Price(r_t, a, sigma, t, T+delta, yield, forward))/(futilities::sum(1, numPayments+1, [&](const auto& index){
      return Bond_Price(r_t, a, sigma, t, t+delta*index, yield, forward);
    })*delta);
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto Forward_Swap_Rate(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*swap initiation*/
    const ThirdFutureTime& Tm, /*swap maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    int numPayments=floor((Tm-T)/delta)+1; //this should be an integer!  remember, T-t is the total swap length
    return (Bond_Price(r_t, a, sigma, t, T, yield, forward)-Bond_Price(r_t, a, sigma, t, Tm+delta, yield, forward))/(futilities::sum(1, numPayments+1, [&](const auto& index){
      return Bond_Price(r_t, a, sigma, t, T+delta*index, yield, forward);
    })*delta);
  }

  template<typename R, typename MeanRevertSpeed, typename Volatility, typename FirstFutureTime, typename SecondFutureTime, typename Delta, typename SwapRate, typename GetYield, typename GetInstantaneousForward>
  auto Swap_Price( 
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& T, /*swap maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const SwapRate& rate,
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    int numPayments=floor((T-t)/delta)+1; //this should be an integer if t lands on a tenor date
    auto firstExchangeDate=T-(numPayments-1)*delta;    
    return Bond_Price(r_t, a, sigma, t, firstExchangeDate, yield, forward)-futilities::sum(1, numPayments, [&](const auto& index){
      return Bond_Price(r_t, a, sigma, t, firstExchangeDate+delta*index, yield, forward)*rate*delta;
    })-(1.0+rate*delta)*Bond_Price(r_t, a, sigma, t, firstExchangeDate+delta*(numPayments), yield, forward);

  }


  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto Payer_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor...length of time for swap at TM*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    assert(TM>t);
    int numPayments=floor(SwapTenor/delta)+1; 
    auto couponTimes=futilities::for_each_parallel(1, numPayments+1, [&](const auto& index){
      return TM+index*delta;
    });
    return Coupon_Bond_Put(r_t, a, sigma, 1.0, t, TM, couponTimes, swapRate*delta, yield, forward);//swaption is equal to put on coupon bond with coupon=swaption swapRate*delta and strike 1.
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto Receiver_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor...length of time for swap at TM*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    assert(TM>t);
    int numPayments=floor(SwapTenor/delta)+1; 
    auto couponTimes=futilities::for_each_parallel(1, numPayments+1, [&](const auto& index){
      return TM+index*delta;
    });
    return Coupon_Bond_Call(r_t, a, sigma, 1.0, t, TM, couponTimes, swapRate*delta, yield, forward);//swaption is equal to put on coupon bond with coupon=swaption swapRate*delta and strike 1.
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto American_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward,
    bool isPayer,
    int num_steps
  ){
      assert(TM>t);
      auto alphaFunction=[&](double t1, const auto& currVal, double dt, int j){
          return -(a*currVal)/sigma;
      };
      auto sigmaFunction=[&](double t1, const auto& currVal, double dt, int j){
          return 0.0;//interesting  
      };
      auto fInv=[&](double t1, const auto& y, double dt, int j){
          return sigma*y;//+phi(t+t1);
      };
      auto checkTreeStep=[&](int step, int j, double t1, double phiAtT){
        return step!=j?phiT(a, sigma, t1, forward):phiAtT;
      };
      int treeStep=-1;
      double phiAtT=0;
      auto payoff=[&](double t1, const auto& currVal, double dt, int j){
          t1=t1+t;
          phiAtT=checkTreeStep(treeStep, j, t1, phiAtT);
          //phiAtT=phiT(a, sigma, t1, forward);
          treeStep=j;
          auto swp=Swap_Price(currVal+phiAtT, a, sigma, t1, SwapTenor+t1, delta, swapRate, yield, forward);
          return isPayer?(swp>0?swp:0.0):(swp<0?0.0:swp);
      };
      auto discount=[&](double t1, const auto& currVal, double dt, int j){
          phiAtT=checkTreeStep(treeStep, j, t+t1, phiAtT);
          treeStep=j;
          return exp(-(currVal+phiAtT)*dt);  
      };
      return btree::computePrice(
        alphaFunction, 
        sigmaFunction, 
        fInv, 
        payoff, 
        discount, 
        (r_t-phiT(a, sigma, t, forward))/sigma,//initial "y"
        num_steps, 
        TM-t
      );
  }


  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto American_Payer_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward,
    int num_steps
  ){
    return American_Swaption(r_t, a, sigma, swapRate, t, SwapTenor, TM, delta, yield, forward, true, num_steps);
  }
  /**Defaults to 100 steps*/
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto American_Payer_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    return American_Payer_Swaption(r_t, a, sigma, swapRate, t, SwapTenor, TM, delta, yield, forward, 100);
  }
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto American_Receiver_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward,
    int num_steps
  ){
    return American_Swaption(r_t, a, sigma, swapRate, t, SwapTenor, TM, delta, yield, forward, false, num_steps);
  }
   /**Defaults to 100 steps*/
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto American_Receiver_Swaption(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
    return American_Receiver_Swaption(r_t, a, sigma, swapRate, t, SwapTenor, TM, delta, yield, forward, 100);
  }


  /**Intended for testing only!!*/
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto PayerSwaptionWithTree(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
      assert(TM>t);
      auto alphaFunction=[&](double t1, const auto& currVal, double dt, int j){
          return -(a*currVal)/sigma;
      };
      auto sigmaFunction=[&](double t1, const auto& currVal, double dt, int j){
          return 0.0;//interesting  
      };
      auto fInv=[&](double t1, const auto& y, double dt, int j){
          return sigma*y;//+phi(t+t1);
      };
      //to save recomputing along all nodes
      auto checkTreeStep=[&](int step, int j, double t, double phiAtT){
        return step!=j?phiT(a, sigma, t, forward):phiAtT;
      };
      int treeStep=-1;
      double phiAtT=0;
      auto payoff=[&](double t1, const auto& currVal, double dt, int width){
          t1=t1+t;
          phiAtT=checkTreeStep(treeStep, width, t1, phiAtT);
          //phiAtT=phiT(a, sigma, t1, forward);
          treeStep=width;
          auto swp=Swap_Price(currVal+phiAtT, a, sigma, t1, t1+SwapTenor, delta, swapRate, yield, forward);
          return swp>0?swp:0.0;
      };
      auto discount=[&](double t1, const auto& currVal, double dt, int width){
          phiAtT=checkTreeStep(treeStep, width, t+t1, phiAtT);
          //phiAtT=phiT(a, sigma, t1+t, forward);
          treeStep=width;
          return exp(-(currVal+phiAtT)*dt);  
      };
      return btree::computePrice(
        alphaFunction, 
        sigmaFunction, 
        fInv, 
        payoff, 
        discount, 
        (r_t-phiT(a, sigma, t, forward))/sigma,//initial "y"
        5000, 
        TM-t,
        false
      );
  }
  /**Intended for testing only!!*/
  template<typename R, typename MeanRevertSpeed, typename Volatility, typename Strike, typename FirstFutureTime, typename SecondFutureTime, typename ThirdFutureTime, typename Delta, typename GetYield, typename GetInstantaneousForward>
  auto ReceiverSwaptionWithTree(
    const R& r_t,
    const MeanRevertSpeed& a,
    const Volatility& sigma,
    const Strike& swapRate,
    const FirstFutureTime& t, /*future time*/
    const SecondFutureTime& SwapTenor, /*swap tenor*/
    const ThirdFutureTime& TM, /*option maturity*/
    const Delta& delta, /*tenor of the floating rate*/
    const GetYield& yield,
    const GetInstantaneousForward& forward
  ){
      assert(TM>t);
      auto alphaFunction=[&](double t1, const auto& currVal, double dt, int j){
          return -(a*currVal)/sigma;
      };
      auto sigmaFunction=[&](double t1, const auto& currVal, double dt, int j){
          return 0.0;//interesting  
      };
      auto fInv=[&](double t1, const auto& y, double dt, int j){
          return sigma*y;//+phi(t+t1);
      };
      //to save recomputing along all nodes
      auto checkTreeStep=[&](int step, int j, double t, double phiAtT){
        return step!=j?phiT(a, sigma, t, forward):phiAtT;
      };
      int treeStep=-1;
      double phiAtT=0;
      auto payoff=[&](double t1, const auto& currVal, double dt, int width){
          t1=t1+t;
          phiAtT=checkTreeStep(treeStep, width, t1, phiAtT);
          //phiAtT=phiT(a, sigma, t1, forward);
          treeStep=width;
          auto swp=Swap_Price(currVal+phiAtT, a, sigma, t1, t1+SwapTenor, delta, swapRate, yield, forward);
          return swp<0?0.0:swp;
      };
      auto discount=[&](double t1, const auto& currVal, double dt, int width){
          phiAtT=checkTreeStep(treeStep, width, t+t1, phiAtT);
          //phiAtT=phiT(a, sigma, t1+t, forward);
          treeStep=width;
          return exp(-(currVal+phiAtT)*dt);  
      };
      return btree::computePrice(
        alphaFunction, 
        sigmaFunction, 
        fInv, 
        payoff, 
        discount, 
        (r_t-phiT(a, sigma, t, forward))/sigma,//initial "y"
        5000, 
        TM-t,
        false
      );
  }
  /**generates a vasicek given a, b, sigma, and normal random number*/
  template<typename R, typename MeanRevertSpeed, typename LongRunMean, typename Volatility, typename FirstFutureTime, typename Simul>
  auto generateVasicek(const R& currVal,  const MeanRevertSpeed& a, const LongRunMean& b, const Volatility& sigma, const FirstFutureTime& nextTime, const Simul& simul){
      auto tmp=exp(-a*nextTime);
      return b*(1-tmp)+currVal*tmp+sigma*sqrt((1-exp(-2*a*nextTime))/(2*a))*simul;
  }
}
#endif
