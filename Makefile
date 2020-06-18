# Extremely crude stop-gap makefile

CFLAGS := -I.. `root-config --libs --cflags` -I${EIGEN_INC} -I${BOOST_INC} -I${GSL_INC} -L${GSL_LIB} -lgsl -lgslcblas -g -O3

# Uncomment to enable (broken) stan build
#CFLAGS += -DOSCLIB_STAN -I${STAN_MATH_INC}

SRCS := $(wildcard *.cxx)
HDRS := $(wildcard *.h)
OBJS := $(patsubst %.cxx,%.o,$(SRCS))

%.o: %.cxx ${HDRS}
	g++ $< ${CFLAGS} -c -fpic

all: ${OBJS}
	g++ -shared -o libOscLib.so *.o ${CFLAGS}
	+make -C test

clean:
	rm *.o *.so
