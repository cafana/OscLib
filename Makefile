# Extremely crude stop-gap makefile

CFLAGS := -I.. $(shell root-config --cflags) -I${EIGEN_INC} -I${BOOST_INC} -I${GSL_INC} -g -O3

LDFLAGS := $(shell root-config --libs) -L${GSL_LIB} -lgsl -lgslcblas

# Uncomment to enable (broken) stan build
#CFLAGS += -DOSCLIB_STAN -I${STAN_MATH_INC}

SRCS := $(wildcard *.cxx)
HDRS := $(wildcard *.h)
OBJS := $(patsubst %.cxx,%.o,$(SRCS))

%.o: %.cxx ${HDRS}
	g++ $< ${CFLAGS} -c -fpic

all: ${OBJS}
	g++ -shared -o libOscLib.so *.o ${LDFLAGS}
	+make -C test

clean:
	rm *.o *.so
