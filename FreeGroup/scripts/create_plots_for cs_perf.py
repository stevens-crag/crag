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
	reader = csv.reader(csvfile)
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
v_num_id = indices['vertices_num']

def y_average_value(pairs):
		values = {}
		for u in {pair[0] for pair in pairs}:
			values[u] = np.mean([pair[1] for pair in pairs if pair[0] == u])

		return values


def plot_progression(ax, y_id, fixed_value):
	u_time_pairs = [(row[u_id], row[y_id]) for row in rows if row[v_id] == fixed_value and row[c_id] == fixed_value]
	v_time_pairs = [(row[v_id], row[y_id]) for row in rows if row[u_id] == fixed_value and row[c_id] == fixed_value]
	c_time_pairs = [(row[c_id], row[y_id]) for row in rows if row[u_id] == fixed_value and row[v_id] == fixed_value]

	u_sizes = y_average_value(u_time_pairs)
	v_sizes = y_average_value(v_time_pairs)
	c_sizes = y_average_value(c_time_pairs)


	u_sizes = zip(*u_sizes.items())
	v_sizes = zip(*v_sizes.items())
	c_sizes = zip(*c_sizes.items())

	ax.plot(u_sizes[0], u_sizes[1], label='u')
	ax.plot(v_sizes[0], v_sizes[1], label='v')
	ax.plot(c_sizes[0], c_sizes[1], label='c')
	ax.legend(loc='best')
	ax.set_title("fixed_value={0}".format(fixed_value))

def plot_equals(ax, y_id):
	values = [(row[u_id], row[y_id]) for row in rows if row[u_id] == row[v_id] == row[c_id]]

	avgs = y_average_value(values)

	avgs = zip(*avgs.items())

	ax.plot(avgs[0], avgs[1])


	

#plotting performances
fig, axes = plt.subplots(2, 3, sharex = True, sharey = True) 
fig.figsize = (4, 6)

plt.ylabel("time, ms")
plt.xlabel("parameter value")

plt.title("Average calculation time for keys")
plot_progression(axes[0][0], time_id, 1)
plot_progression(axes[0][1], time_id, 2)
plot_progression(axes[0][2], time_id, 3)
plot_progression(axes[1][0], time_id, 4)
plot_progression(axes[1][1], time_id, 5)
plot_progression(axes[1][2], time_id, 6)

plt.savefig("progression_times")
plt.close()

fig, axes = plt.subplots(1, 1, sharex = True, sharey = True) 
fig.figsize = (2, 2)

plt.ylabel("time, ms")
plt.xlabel("parameter value")

plt.title("Average calculation time for keys when all values are equal each other")
plot_equals(axes, time_id)

plt.savefig("progression_equal_time")
plt.close()





fig, axes = plt.subplots(2, 3, sharex = True, sharey = True) 
fig.figsize = (4, 6)

plt.ylabel("vertices num")
plt.xlabel("parameter value")

plt.title("Average vertices_num time for keys")
plot_progression(axes[0][0], v_num_id, 1)
plot_progression(axes[0][1], v_num_id, 2)
plot_progression(axes[0][2], v_num_id, 3)
plot_progression(axes[1][0], v_num_id, 4)
plot_progression(axes[1][1], v_num_id, 5)
plot_progression(axes[1][2], v_num_id, 6)

plt.savefig("progression_vertices")
plt.close()

fig, axes = plt.subplots(1, 1, sharex = True, sharey = True) 
fig.figsize = (2, 2)

plt.ylabel("time, ms")
plt.xlabel("parameter value")

plt.title("Average vertices num for keys when all values are equal each other")
plot_equals(axes, v_num_id)

plt.savefig("progression_vertices_equal_time")
plt.close()

#trivial vertices num
count = len(rows)
trivial_count = 0
for row in rows:
	if row[v_num_id] == 0:
		trivial_count += 1
print trivial_count, count


# [row in rows if row[indices[]]]

# |u|,|v|,|c|,time,total_length,height,vertices_num