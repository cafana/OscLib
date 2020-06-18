# NB, the resulting executables must be run from the directory above so they
# can find libOscLib

CFLAGS := -I../.. `root-config --cflags --libs` -I${EIGEN_INC} -L.. -lOscLib -g -O3

# Uncomment to enable (broken) stan build
#CFLAGS += -DOSCLIB_STAN -I${STAN_MATH_INC}

all:
	g++ testOsc.cc ${CFLAGS} -o testOsc
	g++ testOscRndm.cc ${CFLAGS} -o testOscRndm
	g++ testSterile.cc ${CFLAGS} -o testSterile
	g++ testEarthModel.cc ${CFLAGS} -o testEarthModel
