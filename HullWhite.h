#ifndef __HULLWHITE_H_INCLUDED__
#define __HULLWHITE_H_INCLUDED__
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <vector>
#include <unordered_map>
//#include <functional>
#include "Tree.h"
#include "AutoDiff.h"
#include "BondUtilities.h"
#include "BlackScholes.h" //for option pricing
#include "Newton.h"

/*Note: the fundamental times here are (0, t, T, TM).  0 is current time (and is reflective of the current yield curve), t is some future time that we may want to price options at given the underlying at that time, T is an "initial" maturity and TM a "Final" maturity.  While it is natural to think of (0<t<T<TM), I only require 0<t and 0<T<TM. Note that ALL TIMES ARE WITH RESPECT TO 0!  */

auto T_Forward_Bond_Volatility( /*This is the volatility of a forward bond B(t, TM)/(B(t, T)), TM>T and can be used for zero coupon bond pricing using BS formula.  TM may equal T+delta, but must be manually added before calling this function*/
  const auto&,/*Speed of mean reversion ("a")*/
  const auto&, /*interest rate volatility ("sigma") */
  const auto&, /*future time t */
  const auto&, /*Initial Bond maturity T*/
  const auto& /*Final Bond maturity TM*/
);
auto A( /*A(t, T) from the Hull White PDE*/
  const auto&,/*a*/
  const auto&,/*t*/
  const auto& /*T */
);
auto A( /*A(T-t) from the Hull White PDE*/
  const auto&,/*a*/
  const auto& /*T-t*/
);
auto C( /*C(t, T) from the Hull White PDE...this may require taking a reference to a yield class!*/
  const auto&, /*a */
  const auto&, /*a */
  const auto&,/*t*/
  const auto&, /*T */
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto phiT(
    const auto&, /*speed of mean reversion*/
    const auto&, /*interst rate volatility*/
    const auto&, /*expectation horizon*/
    auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto muR(/*the expected value of r_t under risk neutral measure: E[r_T|t]*/
  const auto&, /*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*interst rate volatility*/
  const auto&, /*future time*/
  const auto&, /*expectation horizon*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto varianceR(/*the variance of r_t under risk neutral measure */
  const auto&, /*speed of mean reversion*/
  const auto&, /*interst rate volatility*/
  const auto&, /*future time*/
  const auto& /*expectation horizon*/
);
auto Bond_Price(/*The zero coupon bond price under Hull White*/
  const auto&, /*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*future time*/
  const auto&, /*Bond expiration*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Bond_Price(/*The zero coupon bond price under Hull White for t=0*/
  const auto&, /*Bond expiration*/
  auto& /*This is a yield class passed here...this should include member functions "Forward(t)" and "Yield(t)" and these should be the instantanoues forward rate and the continuously compounded zero coupon (yield*t) at time 0*/
);
auto Coupon_Bond_Price(/*The coupon bond price under Hull White*/
  const auto&, /*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*future time*/
  const std::vector<auto>&, /*these are coupon times FROM 0!  these should be in order but dont have to be*/
  const auto&, /*coupon rate*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Coupon_Bond_Price(/*The coupon bond price under Hull White at time 0*/
  const std::vector<auto>& , /*these are coupon times FROM 0!  these should be in order though not required*/
  const auto&, /*coupon rate*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Bond_Call(/*The price of a call option on zero coupon bond under Hull White*/
  const auto&, /*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*future time*/
  const auto&, /*option maturity*/
  const auto&, /*bond maturity*/
  const auto&, /*strike*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Bond_Call(/*The price of a call option on zero coupon bond under Hull White at t=0*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*option maturity*/
  const auto&, /*bond maturity*/
  const auto&, /*strike*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Coupon_Bond_Call(/*The price of a call option on coupon bond under Hull White...uses jamshidian's trick*/
  const auto&, /*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*strike*/
  const auto&, /*future time*/
    const auto&, //option maturity
  const std::vector<auto>&, /*these are coupon times FROM 0!  these should be in order but dont have to be*/
  const auto&, /* coupon rate */
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Bond_Put(/*The price of a Put option on zero coupon bond under Hull White*/
  const auto&, /*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*future time*/
  const auto&, /*option maturity*/
  const auto&, /*bond maturity*/
  const auto&, /*strike*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Bond_Put(/*The price of a Put option on zero coupon bond under Hull White at t=0*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*option maturity*/
  const auto&, /*bond maturity*/
  const auto&, /*strike*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Coupon_Bond_Put(/*The price of a put option on coupon bond under Hull White...uses jamshidian's trick*/
  const auto&,/*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*strike*/
  const auto&, /*future time*/
    const auto&, //option maturity
  const std::vector<auto>&, /*these are coupon times FROM 0!  these should be in order but dont have to be*/
  const auto&, /*coupon rate*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Caplet(/*price of a caplet on simple bond yield*/
  const auto&,/*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility*/
  const auto&, /*future time*/
  const auto&, /*option maturity*/
  const auto&, /*tenor of the simple yield*/
  const auto&, /*strike*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Caplet(/*price of a caplet on simple bond yield at time t=0*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility */
  const auto&, /*option maturity*/
  const auto&, /*tenor of the simple yield*/
  const auto&, /*strike*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto EuroDollarFuture(
  const auto&,/*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility */
  const auto&, /*future time*/
  const auto&, /*maturity*/
  const auto&, /*tenor of the simple yield*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Swap_Rate(
  const auto&,/*r_t*/
  const auto&, /*speed of mean reversion*/
  const auto&, /*volatility */
  const auto&, /*future time*/
  const auto&, /*swap maturity*/
  const auto&, /*tenor of the floating rate*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
template<typename optionMaturity>
auto Swaption(
  const auto&, /*r_t*/
  const auto&,/*speed of mean reversion*/
  const auto&, /*volatility */
  const auto&, /*strike*/
  const auto&, /*future time*/
  const auto&, /*swap maturity*/
  const optionMaturity&, /*option maturity*/
  const auto&, /*tenor of the floating rate*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);
auto Swap_Price( 
  const auto&, /*r_t*/
  const auto&,/*speed of mean reversion*/
  const auto&, /*volatility */
  const auto&, /*future time*/
  const auto&, /*swap maturity*/
  const auto&, /*tenor of the floating rate*/
  const auto&, /*swap rate*/
  auto& /*This is a yield class passed here...this should include member functions "Forward" and "Yield" and these should be the instantanoues forward rate and the continuously compounded zero coupon yield*/
);

template<typename optionMaturity> /* */
auto AmericanSwaption(
  const auto&,
  const auto&,
  const auto&,
  const auto&,
  const auto&, /*future time*/
  const auto&, /*swap maturity*/
  const optionMaturity&, /*option maturity*/
  const auto&, /*tenor of the floating rate*/
  auto&
);

#include "HullWhite.hpp"



#endif
