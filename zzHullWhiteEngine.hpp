template<typename Number> 
HullWhiteEngine<Number>::HullWhiteEngine(){
    
}
/*HullWhiteEngine::HullWhiteEngine(Number& a_, Number& sigma_){
    a=a_;
    sigma=sigma_;
}*/
template<typename Number> 
void HullWhiteEngine<Number>::setSigma(Number& sigma_){
    sigma=sigma_;
}
template<typename Number> 
void HullWhiteEngine<Number>::setReversion(Number& a_){
    a=a_;
}
/*template<typename Number> 
void HullWhiteEngine<Number>::setYield(auto *yld_){
    yld=yld_;
    
}*/
/*template<typename Number> 
void HullWhiteEngine<Number>::setShortRate(auto& rate){
    r=rate;
}*/
template<typename Number> 
auto HullWhiteEngine<Number>::HullWhitePrice(
        AssetFeatures& asset, 
        auto& rate, 
        Date& futureTime,
        Date& asOfDate,
        YieldSpline& yld
    ){
    std::vector<double> coupons;
    //std::cout<<rate<<std::endl;
    switch(asset.type){
            
        case BOND:
            return Bond_Price(rate, a, sigma, futureTime-asOfDate, asset.Maturity-asOfDate, yld);
        case EURODOLLARFUTURE:
            return EuroDollarFuture(rate, a, sigma, futureTime-asOfDate, asset.Maturity-asOfDate, asset.Tenor, yld);
        case COUPONBOND:
           // std::vector<double> coupons;
            for(auto& it:asset.Coupons){
                coupons.push_back(it-asOfDate);
            }
            return Coupon_Bond_Price(rate, a, sigma, futureTime-asOfDate, coupons, asset.CouponRate, yld);
        case BONDCALL:
            return Bond_Call(rate, a, sigma, futureTime-asOfDate, asset.Maturity-asOfDate, asset.UnderlyingMaturity-asOfDate, asset.Strike, yld);
        case BONDPUT:
            return Bond_Put(rate, a, sigma, futureTime-asOfDate, asset.Maturity-asOfDate, asset.UnderlyingMaturity-asOfDate, asset.Strike, yld);
        case COUPONBONDCALL:
           // std::vector<double> coupons;
            for(auto& it:asset.Coupons){
                coupons.push_back(it-asOfDate);
            }
            return Coupon_Bond_Call(rate, a, sigma, asset.Strike, futureTime-asOfDate, asset.UnderlyingMaturity-asOfDate, coupons, asset.CouponRate, yld);
        case COUPONBONDPUT:
            //std::vector<double> coupons;
            for(auto& it:asset.Coupons){
                coupons.push_back(it-asOfDate);
            }
            return Coupon_Bond_Put(rate, a, sigma, asset.Strike, futureTime-asOfDate, asset.UnderlyingMaturity-asOfDate, coupons, asset.CouponRate, yld);
        case CAPLET:
            return Caplet(rate, a, sigma, futureTime-asOfDate, asset.Maturity-asOfDate, asset.Tenor, asset.Strike, yld);
        case SWAP:
            return Swap_Price(rate, a, sigma, futureTime-asOfDate, asset.Maturity-asOfDate, asset.Tenor, asset.Strike, yld);//swap rate is strike here...
        case SWAPTION:
            return Swaption(rate, a, sigma, asset.Strike, futureTime-asOfDate, asset.UnderlyingMaturity-asOfDate, asset.Maturity-asOfDate, asset.Tenor, yld);
        case AMERICANSWAPTION:
            return AmericanSwaption(rate, a, sigma, asset.Strike, futureTime-asOfDate, asset.UnderlyingMaturity-asOfDate, asset.Maturity-asOfDate, asset.Tenor, yld);
            
            
    }
}
