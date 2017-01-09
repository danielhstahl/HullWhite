#ifndef __BLACKSCHOLES_H_INCLUDED__
#define __BLACKSCHOLES_H_INCLUDED__
#include <cmath>

template<typename S, typename Discount, typename K, typename Sigma>
auto BSPut(const S &S0, const Discount &discount, const K &k, const Sigma &sigma){
    if(sigma>0){
        double s=sqrt(2.0);
        auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
        return S0*(.5*erf(d1/s)-.5)+k*discount*(.5-.5*(erf((d1-sigma)/s)));
    }
    else{
        return k>S0?(k-S0)*discount:0.0;
    }
}
template<typename S, typename Discount, typename K, typename Sigma>
auto BSCall(const S &S0, const Discount &discount, const K &k, const Sigma &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
    if(sigma>0){
        double s=sqrt(2.0);
        auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
        return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
    }
    else{
        return S0>k?(S0-k)*discount:0.0;
    }
}

#endif
