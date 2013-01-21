#Build parameters to be changed
CXXFLAGS += -g -std=c++11 -fpermissive # Compiler flags: -g for debug, -O for optimization
LDFLAGS  += # Linker flags
LIBS      = -lm -lgsl -lgslcblas #common libraries
INCLUDE   = # Common include path