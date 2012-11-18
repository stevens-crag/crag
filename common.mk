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

##
##  SET SEARCH PATHs
##

GSL_PATH = /usr/pkg/lib/


vpath %.cpp $(SRC_DIR) $(MAIN_DIR)
vpath %.o $(OBJ_DIR)
vpath %.a $(LIB_DIR)
##
##
##
##  COMPILER PARAMETERS
##
##

CC      = g++   # 
CFLAGS  = -g  -fpermissive # Compiler flags: -g for debug, -O for optimization
LDFLAGS = -L$(GSL_PATH) -L $(LIB_DIR)/  $(foreach dir,$(DEPEND_ON),-L../$(dir)/$(LIB_DIR))      # Linker flags
LIBS    = -lm -lgsl -lgslcblas  -l$(THISLIB) $(foreach name,$(DEPEND_ON),-l$(name)) $(LOCAL_LIBS)
INCLUDE = -I include/ $(foreach dir,$(DEPEND_ON),-I../$(dir)/include) $(LOCAL_INCLUDE)
DEPENDONLIBS = $(DEPEND_ON)
OBJFILES = $(foreach file,$(SRC),$(OBJ_DIR)/$(file).o)


#############################################################################
##
## COMPILATION
##
#############################################################################

#   $@ stands for the target name (i.e., the resulting executable)
#   $? stands for the dependency list (i.e., the .o files)
#

.PHONY:  all lib  $(THISLIB) lib$(THISLIB)  clean $(MAIN)


# These are all the files to be compiled.
ALL     = lib  $(MAIN)
all:    $(ALL)
lib lib$(THISLIB):	$(LIB_DIR)/lib$(THISLIB).a

## Make local library
$(LIB_DIR)/lib$(THISLIB).a : $(OBJFILES)
	if [ ! -d lib ]; then mkdir lib; fi
	$(AR) cr $@ $(OBJFILES)
	ranlib $@


# Compile the executable targets
$(MAIN) : % :   $(OBJ_DIR)/%.o  $(LIB_DIR)/lib$(THISLIB).a $(DEPENDONLIBS)
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	@echo "$(DEPENDONLIBS)"
	$(CC) $(CFLAGS) $(INCLUDE) $< $(LDFLAGS) $(LIBS)    -o $(BIN_DIR)/$@
	@echo
	@echo " ./$@ : compiled sucessfully."
	@echo

# Make libs that current library depends on
$(DEPENDONLIBS): % :
	cd ../$*; $(MAKE) -f Makefile lib


## compile object files
$(OBJ_DIR)/%.o : %.cpp
	${CC} ${CFLAGS} $(INCLUDE) -c $< -o $@


## Build dependencies
$(OBJ_DIR)/%.d: %.cpp
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	$(SHELL) -ec '$(CC) -MM $(INCLUDE) $< | sed "s:$*.o:$(OBJ_DIR)/& $@:g"' > $@

# Include all the dependency files, this will  automatically remake them if they are
# out of date
include $(foreach file,$(SRC) $(MAIN),$(OBJ_DIR)/$(file).d)


# Clean target to remove backup, object, and core files
clean:
	rm -f  $(OBJ_DIR)/*.d  $(OBJ_DIR)/*.o $(LIB_DIR)/*.a

