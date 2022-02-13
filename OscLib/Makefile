# BASH-isms are used in shell directives so set the shell to be bash
SHELL = bash

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    DSOEXT = so
endif
ifeq ($(UNAME_S),Darwin)
    DSOEXT = dylib
endif

# Extremely crude stop-gap makefile

# Record what flags were used. Can dump with
# readelf -p .GCC.command.line libOscLib.$(DSOEXT)
CFLAGS := -std=c++17 -I.. -I${ROOT_INC} -I${EIGEN_INC} -I${BOOST_INC} -I${GSL_INC} -g -Wall -Wpedantic -Wextra -Werror -grecord-gcc-switches -fvisibility=default

CFLAGS += $(shell if [[ x${QUALIFIER} == *debug* ]]; then echo ''; else echo '-O3'; fi )

CXX := $(shell if [[ x${QUALIFIER} == *c7* ]]; then echo 'clang++'; else echo 'g++'; fi )

LDFLAGS := -L${ROOT_LIB} -lCore -lHist -L${GSL_LIB} -lgsl -lgslcblas

LDFLAGS_BINS := -lOscLib -L${ROOT_LIB} -lCore -lGraf -lGpad -lHist -lMathCore -lRIO -Llib

# Only try to use stan if the product is set up
ifdef STAN_MATH_INC
CFLAGS += -DOSCLIB_STAN -D_REENTRANT -I${SUNDIALS_INC} -I${STAN_MATH_INC} -I${TBB_INC} -Wno-unused-function

# If we can't find "2019" in the version, ie it's newer than that
ifeq (,$(findstring v2019,${TBB_UPS_VERSION}))
CFLAGS += -DTBB_INTERFACE_NEW
endif

LDFLAGS += -L${TBB_LIB} -ltbb
LDFLAGS_BINS += -L${TBB_LIB} -ltbb
endif

SRCS := $(wildcard *.cxx)
HDRS := $(wildcard *.h)
OBJS := $(patsubst %.cxx,tmp/%.o,$(SRCS))
BIN_SRCS := $(wildcard test/test*.cc)
BINS := $(patsubst test/test%.cc,bin/test%,$(BIN_SRCS))

all: prereqs ${OBJS} ${BINS} lib/libOscLib.$(DSOEXT)

tmp/%.o: %.cxx ${HDRS}
	${CXX} $< ${CFLAGS} -c -fpic -o $@

bin/test%: test/test%.cc ${HDRS} lib/libOscLib.$(DSOEXT)
	${CXX} $< ${CFLAGS} ${LDFLAGS_BINS} -o $@

lib/libOscLib.$(DSOEXT): ${OBJS}
	${CXX} -shared -o $@ tmp/*.o ${LDFLAGS}

prereqs:
	@echo Checking all necessary env vars are set
	test ${ROOT_INC} # ROOT_INC
	test ${EIGEN_INC} # EIGEN_INC
	test ${BOOST_INC} # BOOST_INC
	test ${GSL_INC} # GSL_INC
	test ${ROOTSYS} # ROOTSYS
	test ${GSL_LIB} # GSL_LIB
	test ${ROOT_LIB} # ROOT_LIB
	mkdir -p tmp
	mkdir -p lib
	mkdir -p bin

clean:
	rm tmp/*.o lib/*.$(DSOEXT) bin/*
