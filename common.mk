##         
##
##   Makefile


## (by sasha) -fpermissive was added to CFLAGS 


###########################################################################
##
## GLOBAL VARIABLES
##
###########################################################################
THISLIB = $(shell basename `pwd`)
##
##
##  DIRECTORIES
##
OBJ_DIR = src/obj
BIN_DIR = bin
MAIN_DIR = main
LIB_DIR = lib
SRC_DIR = src
TESTS_DIR = test

##
##  SET SEARCH PATHs
##

vpath %.cpp $(SRC_DIR) $(MAIN_DIR)
vpath %.o $(OBJ_DIR)
vpath %.a $(LIB_DIR)
##
##
##
##  COMPILER PARAMETERS
##
##

include ../build-config.mk

LDFLAGS += -L $(LIB_DIR)/  $(foreach dir,$(DEPEND_ON),-L ../$(dir)/$(LIB_DIR))
LIBS    += -l$(THISLIB) $(foreach name,$(DEPEND_ON),-l$(name)) $(LOCAL_LIBS)
INCLUDE += -I include/ $(foreach dir,$(DEPEND_ON),-I../$(dir)/include) $(LOCAL_INCLUDE)
DEPENDONLIBS = $(DEPEND_ON)
OBJFILES = $(foreach file,$(SRC),$(OBJ_DIR)/$(file).o)
LIBTESTS = $(wildcard $(TESTS_DIR)/*.cpp)

#############################################################################
##
## COMPILATION
##
#############################################################################

#   $@ stands for the target name (i.e., the resulting executable)
#   $? stands for the dependency list (i.e., the .o files)
#

.PHONY:  all lib  $(THISLIB) lib$(THISLIB)  clean $(MAIN) tests check


# These are all the files to be compiled.
ALL     = lib  $(MAIN)
all:    $(ALL)
lib lib$(THISLIB):	$(LIB_DIR)/lib$(THISLIB).a

## Make local library
$(LIB_DIR)/lib$(THISLIB).a : $(OBJFILES)
	if [ ! -d $(LIB_DIR) ]; then mkdir $(LIB_DIR); fi
	$(AR) cr $@ $(OBJFILES)
	ranlib $@


# Compile the executable targets
$(MAIN) : % :   $(OBJ_DIR)/%.o  $(LIB_DIR)/lib$(THISLIB).a $(DEPENDONLIBS)
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	@echo "$(DEPENDONLIBS)"
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) $< $(LDFLAGS) $(LIBS) -o $(BIN_DIR)/$@
	@echo
	@echo " ./$@ : compiled sucessfully."
	@echo

# Make libs that current library depends on
$(DEPENDONLIBS): % :
	cd ../$*; $(MAKE) -f Makefile lib

## compile object files
$(OBJ_DIR)/%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c $< -o $@

## Build dependencies
$(OBJ_DIR)/%.d: %.cpp
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	$(SHELL) -ec '$(CXX) $(CPPFLAGS) -MM $(INCLUDE) $(CXXFLAGS) $< | sed "s:$*.o:$(OBJ_DIR)/& $@:g"' > $@


# Google test framework integration
GTEST_DIR = ../gtest

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h \
                ../general/include/gmp_boost_pool_allocator.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.

$(GTEST_DIR)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR)/include -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $(GTEST_DIR)/gtest-all.o

$(GTEST_DIR)/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR)/include -I$(GTEST_DIR) -I../general/include $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc -o $(GTEST_DIR)/gtest_main.o

$(GTEST_DIR)/gtest_main.a : $(GTEST_DIR)/gtest-all.o $(GTEST_DIR)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

## compile object files for unit tests
$(OBJ_DIR)/unittest_%.o : $(TESTS_DIR)/%.cpp $(GTEST_HEADERS) 
	$(CXX) $(CPPFLAGS) $(INCLUDE) -I$(GTEST_DIR)/include $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/unittest_%.d : $(TESTS_DIR)/%.cpp $(GTEST_HEADERS)
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	$(SHELL) -ec '$(CXX) $(CPPFLAGS) -MM $(INCLUDE) -I$(GTEST_DIR)/include $(CXXFLAGS) $< | sed "s:$*.o:$(OBJ_DIR)/unittest_& $@:g"' > $@

	
tests: $(BIN_DIR)/all_unittests
$(BIN_DIR)/all_unittests: $(LIBTESTS:$(TESTS_DIR)/%.cpp=$(OBJ_DIR)/unittest_%.o) $(GTEST_DIR)/gtest_main.a \
                      $(LIB_DIR)/lib$(THISLIB).a $(DEPENDONLIBS) general
	@echo "found tests: $(LIBTESTS:$(TESTS_DIR)/%.cpp=%)"
	@if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	$(CXX) $(CPPFLAGS) $(INCLUDE)-I$(GTEST_DIR)/include $(CXXFLAGS) $(LIBTESTS:$(TESTS_DIR)/%.cpp=$(OBJ_DIR)/unittest_%.o) $(GTEST_DIR)/gtest_main.a \
	    $(LDFLAGS) -lpthread -lgeneral $(LIBS) -o $@
	@echo
	@echo " ./$@ : compiled sucessfully."
	@echo

check: tests
ifdef FILTER
	$(BIN_DIR)/all_unittests --gtest_filter="$(FILTER)"
else
	$(BIN_DIR)/all_unittests
endif
	
	

# Clean target to remove backup, object, and core files
clean:
	rm -f  $(OBJ_DIR)/*.d  $(OBJ_DIR)/*.o $(LIB_DIR)/*.a $(BIN_DIR)/*
	
# Include all the dependency files, this will  automatically remake them if they are
# out of date
include $(foreach file,$(SRC) $(MAIN) $(LIBTESTS:$(TESTS_DIR)/%.cpp=unittest_%),$(OBJ_DIR)/$(file).d)
