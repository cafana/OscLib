# Extremely crude stop-gap makefile

CFLAGS := -std=c++17 -I.. -I${ROOT_INC} -I${EIGEN_INC} -I${BOOST_INC} -I${GSL_INC} -g -O3 -Wall -Werror

LDFLAGS := -L${ROOTSYS}/lib -lCore -L${GSL_LIB} -lgsl -lgslcblas

LDFLAGS_BINS := -L${ROOTSYS}/lib -lCore -lGraf -lGpad -lHist -lMathCore -lRIO -Llib -lOscLib

# Only try to use stan if the product is set up
ifdef STAN_MATH_INC
CFLAGS += -DOSCLIB_STAN -I${STAN_MATH_INC} -Wno-unused-function
endif

SRCS := $(wildcard *.cxx)
HDRS := $(wildcard *.h)
OBJS := $(patsubst %.cxx,tmp/%.o,$(SRCS))
BIN_SRCS := $(wildcard test/test*.cc)
BINS := $(patsubst test/test%.cc,bin/test%,$(BIN_SRCS))

all: prereqs ${OBJS} ${BINS} lib/libOscLib.so

tmp/%.o: %.cxx ${HDRS}
	g++ $< ${CFLAGS} -c -fpic -o $@

bin/test%: test/test%.cc ${HDRS} lib/libOscLib.so
	g++ $< ${CFLAGS} ${LDFLAGS_BINS} -o $@

lib/libOscLib.so: ${OBJS}
	g++ -shared -o $@ tmp/*.o ${LDFLAGS}

prereqs:
	@echo Checking all necessary env vars are set
	test ${ROOT_INC} # ROOT_INC
	test ${EIGEN_INC} # EIGEN_INC
	test ${BOOST_INC} # BOOST_INC
	test ${GSL_INC} # GSL_INC
	test ${ROOTSYS} # ROOTSYS
	test ${GSL_LIB} # GSL_LIB
	mkdir -p tmp
	mkdir -p lib
	mkdir -p bin

clean:
	rm tmp/*.o lib/*.so bin/*
