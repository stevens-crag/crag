##         Clustering Algorithm
##
##   Makefile

FLEX    = flex++
BISON   = bison
CC      = g++   # 
CFLAGS  = -g  -static  # Compiler flags: -g for debug, -O for optimization
LDFLAGS =  -L lib/ -L ../Group     # Linker flags
LIBS    =    

OBJ_DIR = src/obj/
BIN_DIR = bin/

#   $@ stands for the target name (i.e., the resulting executable)
#   $? stands for the dependency list (i.e., the .o files)
#
COMPILE_EXECUTABLE = $(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDE) -o $(BIN_DIR)$@ $? $(LIBS)

# These are all the files to be compiled.
ALL     = plib 
all:    $(ALL)
lib:	plib
# List of executable targets
#mainParser: src/mainParser.cpp  $(OBJ_DIR) WordBisonGrammar.o  $(OBJ_DIR)WordFlex.o
#	$(COMPILE_EXECUTABLE) 

#$(OBJ_DIR)WordBisonGrammar.o:   src/HierarchicalAlgorithms.c
#	${CC} ${CFLAGS} $(INCLUDE) -c $< -o $@
#$(OBJ_DIR)WordFlex.o:   src/WordFlex.c
#	${CC} ${CFLAGS} $(INCLUDE) -c $< -o $@


## PARSER FILES
plib:	
	$(BISON) -d src/WordBisonGrammar.y -o src/WordBisonGrammar.c;
	$(FLEX)  -osrc/WordFlex.cpp src/WordFlex.l; 
	cd $(OBJ_DIR); ${CC} ${CFLAGS} -I../../include/ -I../ -c ../WordBisonGrammar.c ../WordFlex.cpp; 
	cd ../..
	ar cr lib/libParsers.a $(OBJ_DIR)WordBisonGrammar.o $(OBJ_DIR)WordFlex.o
	ranlib lib/libParsers.a

# Clean target to remove backup, object, and core files
clean:
	rm -f  $(OBJ_DIR)*.o $(BIN_DIR)*  lib/*.a
