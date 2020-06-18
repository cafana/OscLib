# Extremely crude stop-gap makefile

CFLAGS := -std=c++17 -I.. -I${ROOT_INC} -I${EIGEN_INC} -I${BOOST_INC} -I${GSL_INC} -g -O3

LDFLAGS := -L${ROOTSYS}/lib -lCore -L${GSL_LIB} -lgsl -lgslcblas

# Uncomment to enable (broken) stan build
#CFLAGS += -DOSCLIB_STAN -I${STAN_MATH_INC}

SRCS := $(wildcard *.cxx)
HDRS := $(wildcard *.h)
OBJS := $(patsubst %.cxx,%.o,$(SRCS))

%.o: %.cxx ${HDRS}
	g++ $< ${CFLAGS} -c -fpic

all: prereqs ${OBJS}
	g++ -shared -o libOscLib.so *.o ${LDFLAGS}
	+make -C test

prereqs:
	@echo Checking all necessary env vars are set
	test ${ROOT_INC} # ROOT_INC
	test ${EIGEN_INC} # EIGEN_INC
	test ${BOOST_INC} # BOOST_INC
	test ${GSL_INC} # GSL_INC
	test ${ROOTSYS} # ROOTSYS
	test ${GSL_LIB} # GSL_LIB

clean:
	rm *.o *.so
