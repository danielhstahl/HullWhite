LDFLAGS=-L../NewtonOptimization -lNewton -L../DateUtilities -lDate  -L../eigen   -L../AutoDiff -lAutoDiff 
INCLUDES=-I../NewtonOptimization -I../DateUtilities -I../FixedIncomeUtilities -I../eigen -I../BinomialTree -I../rapidjson -I../AutoDiff -I../MonteCarlo

OptionPricing: main.o
	g++ -std=c++14 -O3  -w -fPIC main.o  $(LDFLAGS) $(INCLUDES) -o OptionPricing -fopenmp

main.o: main.cpp BlackScholes.h BlackScholes.hpp HullWhite.h HullWhite.hpp HullWhiteEngine.h HullWhiteEngine.hpp YieldIO.h YieldIO.hpp RealWorldMeasure.h RealWorldMeasure.hpp
	g++ -std=c++14 -O3  -w -c -fPIC main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

clean:
	-rm *.o OptionPricing
