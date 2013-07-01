#! /usr/bin/env python

import csv
import numpy as np
import matplotlib.pyplot as plt
import sys

#change these values
filename = sys.argv[1]

rows = []
indices = {}
with open(filename, 'r') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	first_row = True
	for row in reader:
		row = row[:-1]
		if first_row:
			for i in range(0, len(row)):
				indices[row[i]] = i
			first_row = False
		else:
			rows.append(map(int, row))



u_id = indices['|u|']
v_id = indices['|v|']
c_id = indices['|c|']
time_id = indices['time']



def plot_time_progression(ax, *initial_values):
	u_0, v_0, c_0 = initial_values

	u_time_pairs = [(row[u_id], row[time_id]) for row in rows if row[v_id] == v_0 and row[c_id] == c_0]
	v_time_pairs = [(row[v_id], row[time_id]) for row in rows if row[u_id] == u_0 and row[c_id] == c_0]
	c_time_pairs = [(row[c_id], row[time_id]) for row in rows if row[u_id] == u_0 and row[v_id] == v_0]

	
 
 	def find_average_value(pairs):
		values = {}
		for u in {pair[0] for pair in pairs}:
			values[u] = np.mean([pair[1] for pair in pairs if pair[0] == u]) / 1000

		return values
	
	u_sizes = find_average_value(u_time_pairs)
	v_sizes = find_average_value(v_time_pairs)
	c_sizes = find_average_value(c_time_pairs)


	u_sizes = zip(*u_sizes.items())
	v_sizes = zip(*v_sizes.items())
	c_sizes = zip(*c_sizes.items())

	ax.plot(u_sizes[0], u_sizes[1], label='u')
	ax.plot(v_sizes[0], v_sizes[1], label='v')
	ax.plot(c_sizes[0], c_sizes[1], label='c')
	ax.legend(loc='best')
	ax.set_title("({u},{v},{c})".format(u=u_0, v=v_0, c=c_0))
	

#plotting performances
fig, axes = plt.subplots(1, 2, sharex = True, sharey = True) 

plt.ylabel("time, s")
plt.xlabel("parameter value")

plt.title("Average calculation time for keys")
plot_time_progression(axes[0], 1, 1, 1)
plot_time_progression(axes[1], 2, 2, 2)

plt.show()
plt.close()



# [row in rows if row[indices[]]]

# |u|,|v|,|c|,time,total_length,height,vertices_num