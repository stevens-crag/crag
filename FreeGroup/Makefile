
##
## Authors:Alexei Miasnikov
##
## 
##


##
## DEPEND_ON: list all dependencies here
##

DEPEND_ON = Maps Alphabet Graph SbgpFG  Elt general  ranlib

#
# If you need any extra DEFINES, put them here

override DEFINE +=

#
# Set debug level here

# DEBUG =

#
# if you want to pass extra libraries to the linker which are not
# used by default
#

LOCAL_LIBS = -L/lib/ -Wl,-Bstatic -lgmpxx -lgmp -Wl,-Bdynamic


# 
# SRC: lists all source files
SRC = WhiteheadGraph FreeGroup StraightLineProgramWord StraightLineProgramWord_operators StraightLineProgramWord_truncate StraightLineProgramWord_split SLPSet slp


# 
# MAIN lists all source files with main function:

MAIN = test_cs


## All the compiler and linker stuff here:                           
#                                                                    
include ../common.mk                                          



