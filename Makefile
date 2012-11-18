
.PHONY: all lib clean

all: 
	cd ranlib; $(MAKE) all
	cd general; $(MAKE) all
	cd Elt; $(MAKE) all
	cd Alphabet; $(MAKE) all
	cd Graph; $(MAKE) all
	cd Group; $(MAKE) all
	cd SbgpFG; $(MAKE) all
	cd FreeGroup; $(MAKE) all
	cd HigmanGroup; $(MAKE) all
	cd BraidGroup; $(MAKE) all
	cd CryptoKL; $(MAKE) all
	cd CryptoAAG; $(MAKE) all
	cd CryptoAE; $(MAKE) all
	cd CryptoShftConj; $(MAKE) all
	cd CryptoTripleDecomposition; $(MAKE) all
	cd Equation; $(MAKE) all
	cd Maps; $(MAKE) all
	cd Examples; $(MAKE) all
	cd Graphics; $(MAKE) all 
	cd StringSimilarity; $(MAKE) all
	cd TheGrigorchukGroup; $(MAKE) all
lib:
	cd ranlib; $(MAKE) lib
	cd general; $(MAKE) lib
	cd Elt; $(MAKE) lib
	cd Alphabet; $(MAKE) lib
	cd Graph; $(MAKE) lib
	cd Group; $(MAKE) all
	cd SbgpFG; $(MAKE) lib
	cd FreeGroup; $(MAKE) lib
	cd HigmanGroup; $(MAKE) all
	cd BraidGroup; $(MAKE) lib
	cd CryptoKL; $(MAKE) lib
	cd CryptoShftConj; $(MAKE) lib
	cd CryptoAAG; $(MAKE) lib
	cd CryptoAE; $(MAKE) lib
	cd CryptoTripleDecomposition; $(MAKE) lib
	cd Equation; $(MAKE) lib
	cd Examples; $(MAKE) lib
	cd Graphics; $(MAKE) lib 
	cd StringSimilarity; $(MAKE) lib
	cd Maps; $(MAKE) lib
	cd TheGrigorchukGroup; $(MAKE) lib
clean:
	rm -f -R Release
	if [ -d ranlib ]; then cd ranlib && $(MAKE) clean; fi
	if [ -d general ]; then cd general && $(MAKE) clean; fi
	if [ -d Elt ]; then cd Elt && $(MAKE) clean; fi
	if [ -d Alphabet ]; then cd Alphabet && $(MAKE) clean; fi
	if [ -d Graph ]; then cd Graph && $(MAKE) clean; fi
	if [ -d Group ]; then cd Group && $(MAKE) clean; fi
	if [ -d SbgpFG ]; then cd SbgpFG && $(MAKE) clean; fi
	if [ -d FreeGroup ]; then cd FreeGroup && $(MAKE) clean; fi
	if [ -d HigmanGroup ]; then cd HigmanGroup && $(MAKE) clean; fi
	if [ -d BraidGroup ]; then cd BraidGroup && $(MAKE) clean; fi
	if [ -d CryptoKL ]; then cd CryptoKL && $(MAKE) clean; fi
	if [ -d CryptoAAG ]; then cd CryptoAAG && $(MAKE) clean; fi
	if [ -d CryptoAE ]; then cd CryptoAE && $(MAKE) clean; fi
	if [ -d CryptoShftConj ]; then cd CryptoShftConj && $(MAKE) clean; fi
	if [ -d CryptoTripleDecomposition ]; then cd CryptoTripleDecomposition && $(MAKE) clean; fi
	if [ -d Equation ]; then cd Equation && $(MAKE) clean; fi
	if [ -d Examples ]; then cd Examples && $(MAKE) clean; fi
	if [ -d Experiments ]; then cd Experiments && $(MAKE) clean; fi
	if [ -d Maps ]; then cd Maps && $(MAKE) clean; fi
	if [ -d Examples ]; then cd Examples && $(MAKE) clean; fi
	if [ -d Graphics ]; then cd Graphics && $(MAKE) clean; fi
	if [ -d StringSimilarity ]; then cd StringSimilarity && $(MAKE) clean; fi
	if [ -d TheGrigorchukGroup ]; then cd TheGrigorchukGroup && $(MAKE) clean; fi
	if [ -d python ]; then cd python && $(MAKE) clean; fi
