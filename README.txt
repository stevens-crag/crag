  
	CRyptography And Groups (CRAG) C++ and Python Library  
 
Copyright (C) 2004-2009 The Algebraic Cryptography Center at Stevens. 
 
  
GENERAL DESCRIPTION 
  
The Cryptography And Groups (CRAG) Library provides an environment to test   
cryptographic protocols constructed from non-commutative groups, for example the   
braid group. The Library is written in C++ and provides an interface and   
routines for computations. There are implementations of basic algebraic objects   
like words, maps and subgroups. We plan to continually expand the list of group-  
theoretic algorithms implemented in the library. In addition the Library will   
contain classes and routines implementing non-classical heuristic approaches and   
tools to perform statistical and exploratory analysis of algebraic data.   
Together with the C++ source code CRAG will include interface to Python   
scripting language.   
Third party contributions and cooperation in the areas appropriate to the goals   
of the project are most welcome.  
  
  
 
 
REQUIREMENTS  
 
You need the following utilities installed:   
 
1.	g++ version 2.3.x or later. Type g++ --version to verify this. This   
package can be downloaded from http://www.gnu.org/.  
 
2.	GNU make (any version). Type make --version to verify this. On some   
systems this command may be called gmake or gnumake. This package can be   
downloaded from http://www.gnu.org/.  
 
3.	The GNU Scientific Library (GSL) .  This package can be downloaded from   
http://www.gnu.org/software/gsl/.  

4.      GNU Multiple Precision Arithmetic Library (GMP). This library is not required
to compile basic directories but if you want complete functionality, you should install 
it either using you favourite package manager or directly from here:
http://gmplib.org/
 
The CRAG library is written to comply with ISO C++ standard so it should be   
successfully compiled using most standard C++ compilers.  
  
 
 
 
INSTALLATION 
 
1.	Download the archive and extract files using 'tar' command :  
	tar xfvz crag_x.xx.tgz  
You should have the directory  crag_x.xx containing the library   
source created.  
 
2.	Check that the GSL (GNU Scientific Library) distribution is available.  
Change GSL_PATH variable in common.mk if needed.  
 
3.	Type   'make' (or gmake ).   It should compile all the executables  
available in the current distribution.    
 
 
 
 
LIBRARY ORGANIZATION 
  
The library is partitioned into several directories grouping files according to  
the specific tasks they perform or objects they implement. Each directory has  
the following structure:  
	Makefile	-	makefile which defines the source files and 		 
				dependencies on other directories;  
	include/	-	directory with include (.h) files;  
	src/		-  	directory with source (.cpp) files;  
	main/		- 	directory with source code of main functions; 
	bin/ 		-	directory where executables are created;  
	lib/		-	directory where the library file is created. 
  
  
To compile all executables in a particular directory type 'make' inside that  
directory.  
 
To delete all object files, the library file and executables type 'make clean'.  
  
To add a new source file: copy the .cpp file into src/ directory; add its  name   
without extension to the list of src files in Makefile (variable SRC).   
For example, if you have foo1.cpp and foo2.cpp in your src/ directory, then in  
Makefile you should have the following line:  
  
SRC = foo1 foo2  
  
To add a new executable: put the file with source of the main function into  
main/ directory; add the corresponding file name without extension to the list  
of executables in Makefile (variable MAIN).  
For example, to add main1.cpp add:  
MAIN = main  
  
To compile a particular executable,  type 'make executable_name',  where   
executable_name is the name in the list of the variable MAIN.  
  
In case you want to use classes or functions from another directory, you should  
add the name of that directory to the dependence list (variable DEPEND_ON in the  
Makefile).   
    
  
 
 
LIBRARY CONTENTS:  
  
Doc		- Contains class documentation and FAQ  
Examples	- Directory contains a number of examples of using the library  
general         - A set of general purpose classes 
Elt             - Definition and implementation of basic operations of a group word 
Alphabet        - Definitions of finite and infinite alphabets 
Group           - General definition of Finitely Presented group, Advanced  
                  Denh's algorithm  
FreeGroup      	- Definition of a free group, Whitehead's graph.  
Maps            - Definition of maps, random automorphism generation.  
StringSimilarity- String similarity mesures (Hamming, Editing ... )  
BraidGroup    	- Garside normal forms, Birman-Ko-Lee normal forms, Dehornoy 
		  form, definition of Braid group.  
TheGrigorchukGroup - Definition and algorithms for the Grigorchuk group.
ThompsonGroup   - Definition and algorithms for the Thompson group.
FreeMetabelianGroup - Definition and algorithms for Free Metabelian Groups.
HigmanGroup 	- implementation of power circuits and solution to the world problem in Higman's group
Graph           - Definitions and algorithms for graphs and automata.   
Graphics      	- Algorithms for visualization of algebraic objects.  
Equation        - Definition of equations over finitely presented group 
CryptoAAG    	- Implementation of the Arithmetica key exchange protocol 
CryptoShftConj  - Implementation of the Shifted Conjugacy authentication protocol.  
CryptoKL        - Implementation of the Ko-Lee et. al.  Key exchange  protocol.  
CryprtAE        - Implementation of the TTP attack on the  Algebraic Eraser key exchange
SbgpFG         	- Algorithms for subgroups of free groups (membership test,   
		  Nielsen generators e.t.c)   
ranlib          - Wrapper for the ranlib library (Random number generators,   
		  probability distributions)   
 

