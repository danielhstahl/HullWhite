INCLUDES=-I../GaussNewton -I../FixedIncomeUtilities  -I../BinomialTree -I../rapidjson -I../AutoDiff -I../MonteCarlo -I../FunctionalUtilities

test: test.o
	g++ -std=c++14 -O3  -w -fPIC test.o  $(LDFLAGS) $(INCLUDES) -o test -fopenmp

test.o: test.cpp BlackScholes.h HullWhite.h
	g++ -std=c++14 -O3  -w -c -fPIC test.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

clean:
	-rm *.o test
