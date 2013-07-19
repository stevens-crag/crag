#! /usr/bin/env python

import sys
from subprocess import call
import glob

pattern = sys.argv[1]
print "matching", pattern 
for filename in glob.glob(pattern):
	print "processing", filename
	call(["./lba_attack", "--create_csv", filename])