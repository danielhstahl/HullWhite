INCLUDES=-I../GaussNewton  -I../BinomialTree -I../AutoDiff -I../FunctionalUtilities

test: test.o
	g++ -std=c++14 -O3 -pthread --coverage -w -fPIC test.o  $(LDFLAGS) $(INCLUDES) -o test -fopenmp

test.o: test.cpp BlackScholes.h HullWhite.h
	g++ -std=c++14 -O3  -pthread --coverage -w -c -fPIC test.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

clean:
	-rm *.o test
