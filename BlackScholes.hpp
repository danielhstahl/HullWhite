auto BSPut(const auto &S0, const auto &discount, const auto &k, const auto &sigma){
    if(sigma>0){
        double s=sqrt(2.0);
        auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
        return S0*(.5*erf(d1/s)-.5)+k*discount*(.5-.5*(erf((d1-sigma)/s)));
    }
    else{
        if(k>S0){
            return (k-S0)*discount;
        }
        else{
            return 0.0;
        }
    }
}
auto BSCall(const auto &S0, const auto &discount, const auto &k, const auto &sigma){ //note that sigma includes sqrt(t) term so in vanilla BS sigma is equal to volatility*sqrt(T)
    if(sigma>0){
        double s=sqrt(2.0);
        auto d1=log(S0/(discount*k))/(sigma)+sigma*.5;
        return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma)/s)));
    }
    else{
        if(S0>k){
            return (S0-k)*discount;
        }
        else{
            return 0.0;
        }
    }
}
