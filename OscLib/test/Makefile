# NB, the resulting executables must be run from the directory above so they
# can find libOscLib

CFLAGS := -std=c++17 -I../.. -I${ROOT_INC} -I${EIGEN_INC} -g -O3
LDFLAGS := -L${ROOTSYS}/lib -lCore -lGraf -lGpad -lHist -lMathCore -lRIO -L.. -lOscLib

# Uncomment to enable (broken) stan build
#CFLAGS += -DOSCLIB_STAN -I${STAN_MATH_INC}

all:
	g++ testOsc.cc ${CFLAGS} ${LDFLAGS} -o testOsc
	g++ testOscRndm.cc ${CFLAGS} ${LDFLAGS} -o testOscRndm
	g++ testSterile.cc ${CFLAGS} ${LDFLAGS} -o testSterile
	g++ testEarthModel.cc ${CFLAGS} ${LDFLAGS} -o testEarthModel
