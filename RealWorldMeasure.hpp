double findHistoricalMean(std::vector<SpotValue>& historicalRates, double shortRateTime, double a){//assumes that real world rate process follows Vasicek
    int n=historicalRates.size();
    double dt=0;
    double b=0;
    //double shortRateTime=7.0/360.0;
    historicalRates[0].value=convertLiborToContinuous(historicalRates[0].value, shortRateTime);
    for(int i=0; i<(n-1); ++i){
        historicalRates[i+1].date.setScale("year");
        dt=historicalRates[i+1].date-historicalRates[i].date;
        historicalRates[i+1].value=convertLiborToContinuous(historicalRates[i+1].value, shortRateTime); //convert to continuous time
        b+=(historicalRates[i+1].value-exp(-    dt*a)*historicalRates[i].value)/(1.0-exp(-dt*a));
    }
    b=b/(n-1);
    return b;
  //std::cout<<"This is b: "<<b<<std::endl;
}
auto generateVasicek(const auto& currVal, const auto& nextTime, const auto& a, const auto& b, const auto& sigma, auto& simul){
    auto tmp=exp(-a*nextTime);
    return b*(1-tmp)+currVal*tmp+sigma*sqrt((1-exp(-2*a*nextTime))/(2*a))*simul;
}