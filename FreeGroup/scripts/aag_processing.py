#! /usr/bin/env python

import csv
import numpy
import matplotlib.pyplot as plt
import sys
import pprint


__latex_table_line_end__ = '\\\\ \\hline'


def __list_to_latex_table_line__(l):
	str = '{0} '.format(l[0])
	for y in l[1:]:
		if type(y) == numpy.float32 or type(y) == numpy.float64:
			str += '& {0}'.format(y)
		else:
			str += '& {0} '.format(y)
	return str + __latex_table_line_end__


class StatPlot:

	def __init__(self, filename):
		self.rows = []
		self.indices = {}
		with open(filename, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter=';')
			first_row = True
			for row in reader:
				row = row[:-1]
				if first_row:
					for i in range(0, len(row)):
						self.indices[row[i]] = i
					first_row = False
				else:
					self.rows.append(map(int, row))


	def get_aggregates(self, id1, id2, aggr_func, filtered_rows, xs):
		return [(aggr_func([row[id2] for row in filtered_rows if row[id1] == x])) for x in xs]


	def plot_aggregate(self, ax, id1, id2, aggr_func, agg_name, filter = None, map_func = None):
		if filter:
			filtered_rows = [row for row in self.rows if filter(row)]
		else:
			filtered_rows = self.rows
		xs = sorted(list({row[id1] for row in filtered_rows}))
		ys = self.get_aggregates(id1, id2, aggr_func, filtered_rows, xs)
		if map_func:
			ys = map_func(ys)
		ax.plot(xs, ys, label=agg_name)


	def scatter(self, ax, id1, id2, name, filter = None, dots_color='b', map_func = None):
		if filter:
			xs = [row[id1] for row in self.rows if filter(row)]
			ys = [row[id2] for row in self.rows if filter(row)]
		else:
			xs = [row[id1] for row in self.rows]
			ys = [row[id2] for row in self.rows]
		if map_func:
			ys = map_func(ys)
		ax.scatter(xs, ys, color=dots_color, marker='.', label=name)


	def boxplots(self, ax, id1, id2, filter = None, map_func = None):
		if filter:
			filtered_rows = [row for row in self.rows if filter(row)]
		else:
			filtered_rows = self.rows
		xs = sorted(list({row[id1] for row in filtered_rows}))
		ys = [[row[id2] for row in filtered_rows if row[id1] == x] for x in xs]
		if map_func:
			ys = map(map_func, ys)
		ax.boxplot(ys)
		ax.set_xticklabels(xs)


def time_normalizer(time_row):
		return numpy.divide(map(float, time_row), 1000)


def draw_aag_plots(stat, filename):
	# key_length;time;height;vertices_num;
	key_id = stat.indices['key_length']
	
	fig, axes = plt.subplots(3, 1, sharex = True) 
	fig.figsize = (6, 2)
	plt.xlabel("|e|")

	ax = axes[0]
	ax.set_ylabel("time, s")
	# ax.set_xlabel("|e|")
	# stat.plot_aggregate(ax, key_id, stat.indices['time'], numpy.median, 'median', map_func=time_normalizer)
	stat.scatter(ax, key_id, stat.indices['time'], 'samples', map_func=time_normalizer)
	ax.set_xlim([0, ax.get_xlim()[1]])
	ax.set_ylim([0, ax.get_ylim()[1]])
	ax.legend( loc =' best')

	ax = axes[1]
	ax.set_ylabel("vertices num")
	# stat.plot_aggregate(ax, key_id, stat.indices['vertices_num'], numpy.median, 'median')
	stat.scatter(ax, key_id, stat.indices['vertices_num'], 'samples')
	ax.set_xlim([0, ax.get_xlim()[1]])
	ax.set_ylim([0, ax.get_ylim()[1]])
	ax.legend( loc =' best')

	ax = axes[2]
	ax.set_ylabel("height")
	# stat.plot_aggregate(ax, key_id, stat.indices['height'], numpy.median, 'median')
	stat.scatter(ax, key_id, stat.indices['height'], 'samples')
	ax.set_xlim([0, ax.get_xlim()[1]])
	ax.set_ylim([0, ax.get_ylim()[1]])
	ax.legend( loc =' best')

	plt.show()

	plt.savefig("aag_perf")
	plt.close()
		

if __name__ == "__main__":
	filename = sys.argv[1]

	#rank;|e|;vertices_num;height;time;free_red_vn;free_red_h;free_red_time;nf_vn;nf_ht;nf_time;

	stat = StatPlot(filename)	
	draw_aag_plots(stat, filename[:-4])