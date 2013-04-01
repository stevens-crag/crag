#Build parameters to be changed
CXXFLAGS += -O3 -std=c++11 -fpermissive -DNDEBUG -DBOOST_POOL_NO_MT# Compiler flags: -g for debug, -O for optimization
LDFLAGS  += # Linker flags
LIBS      = -lm -lgsl -lgslcblas #common libraries
INCLUDE   = # Common include path
