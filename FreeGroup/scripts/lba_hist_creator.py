#! /usr/bin/env python

import csv
import numpy
import matplotlib.pyplot as plt
import sys
import math


def process_file(filename):
	morph=[]
	conj=[]
	min_conj=[]
	with open(filename, 'r') as csvfile:
		reader = csv.DictReader(csvfile, delimiter=';')
		#checking file validity
		names = set(reader.fieldnames)
		if 'val' not in names or 'conj_val' not in names or 'min_val' not in names:
			print 'Skipping {0} because it has wrong format'.format(filename)
			return

		n = 0
		for row in reader:
		  morph.append(float(row['val']))
		  conj.append(float(row['conj_val']))
		  min_conj.append(float(row['min_val']))

	max1 = max(morph)
	max2 = max(min_conj)
	line = [float(0)]
	if max1 > max2:
		line.append(max1)
	else:
		line.append(max2)

	log_ratio = map(lambda x,y: math.log(x / y, 2), min_conj, morph)
	# print log_ratio

	plt.hist(log_ratio, bins=20)
	# plt.figure(1, figsize=(6, 3))

	# plt.subplot(121)
	# plt.plot(morph, min_conj, '.')
	# # plt.plot(line, line, 'r-')
	# plt.xlabel("Morphism values")
	# plt.ylabel("Minimized values")

	# plt.subplot(122)
	# plt.plot(conj, min_conj, 'r.')
	# plt.xlabel("Conjugation values")
	# plt.ylabel("Minimized values")



	plt.savefig('histogram_' + filename.split(".")[0] + '.pdf', format='pdf')
	# plt.show()

	plt.close()
		# spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
	 #  for row in spamreader:
	 #      print ', '.join(row)

if __name__ == '__main__':
	import glob
	pattern = sys.argv[1]
	print "matching ", pattern 
	for filename in glob.glob(pattern):
		print "processing ", filename
		process_file(filename)
    