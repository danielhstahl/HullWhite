INCLUDES=-I../GaussNewton  -I../BinomialTree -I../AutoDiff -I../FunctionalUtilities -I../TupleUtilities

GCCVAL=g++

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	GCCVAL=g++-7
endif
test: test.o
	$(GCCVAL) -std=c++14 -O3 -pthread --coverage -w -fPIC test.o  $(LDFLAGS) $(INCLUDES) -o test -fopenmp

test.o: test.cpp BlackScholes.h HullWhite.h
	$(GCCVAL) -std=c++14 -O3  -pthread --coverage -w -c -fPIC test.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

clean:
	-rm *.o test *.gcno
