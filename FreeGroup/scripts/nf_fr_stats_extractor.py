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


	def print_stat(self):
		for k, v in self.indices.iteritems():
			if k != "rank" and k != "|e|" and k != "":
				values = [row[v] for row in self.rows]
				print "{0}: mean={1:.3f}, sd={2:.3f}, max={3:.3f}, min={4:.3f}".format(k, numpy.mean(values), numpy.std(values), numpy.max(values), numpy.min(values))

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


	aggregates = (('max', numpy.max),
		('3/4 quantile', lambda row: numpy.percentile(row, [75])[0]),
		('median', numpy.median),
		('mean', numpy.mean),
		('1/4 quantile', lambda row: numpy.percentile(row, [25])[0]),
		('min', numpy.min))


	def stats_to_list(self, id1, id2, filter, map_func = None):
		if filter:
			filtered_rows = [row for row in self.rows if filter(row)]
		else:
			filtered_rows = rows
		xs = sorted(list({row[id1] for row in filtered_rows}))
		rows = []

		rows.append(xs)

		for name, func in StatPlot.aggregates:
			agg = self.get_aggregates(id1, id2, func, filtered_rows, xs)	
			if map_func:
				agg = map_func(agg)
			rows.append([name] + list(agg))

		return rows

	def to_latex_table(self, name, id1, id2, filter, map_func = None):
		lists = self.stats_to_list(id1, id2, filter, map_func)
		format_string = 'l'
		for x in lists[0]:
			format_string += '|c'

		lines = ['\\begin{center}', name, '\\begin{tabular}{' + format_string + '}']
		lines.append('& ' + __list_to_latex_table_line__(lists[0]))
		for l in lists[1:]:
			lines.append(__list_to_latex_table_line__(l))	

		lines.append('\\end{tabular}')	
		lines.append('\\end{center}')
		return lines


	# def plot_mapped_aggregate(self, ax, id1, id2, aggr_func, name, filter, ):
	# 	filtered_rows = [row for row in self.rows if filter(row)]
	# 	xs = sorted(list({row[id1] for row in filtered_rows}))
	# 	ys = [(aggr_func([row[id2] for row in filtered_rows if row[id1] == x])) for x in xs]
	# 	ax.plot(xs, ys, label = name)

def time_normalizer(time_row):
		return numpy.divide(map(float, time_row), 1000)

def print_tables(stat, ranks):
	#print latex tables
	for rank in ranks:
		def rank_filter(row):
			return row[rank_id] == rank

		print 'rank = ', rank

		for line in stat.to_latex_table('Nielsen composition time, ms', size_id, stat.indices['time'], 
			rank_filter):
			print line

		for line in stat.to_latex_table('Nielsen composition vertices num', size_id, stat.indices['vertices_num'], 
			rank_filter):
			print line

		for line in stat.to_latex_table('Nielsen composition height', size_id, stat.indices['height'], 
			rank_filter):
			print line

		for line in stat.to_latex_table('free reduction time, s', size_id, stat.indices['free_red_time'], 
			rank_filter, time_normalizer):
			print line

		for line in stat.to_latex_table('free reduction vertices num', size_id, stat.indices['free_red_vn'], 
			rank_filter):
			print line

		for line in stat.to_latex_table('free reduction height', size_id, stat.indices['free_red_h'], 
			rank_filter):
			print line

		for line in stat.to_latex_table('normal form time, s', size_id, stat.indices['nf_time'], 
			rank_filter, time_normalizer):
			print line

		for line in stat.to_latex_table('normal form vertices num', size_id, stat.indices['nf_vn'], 
			rank_filter):
			print line

		for line in stat.to_latex_table('normal form height', size_id, stat.indices['nf_ht'], 
			rank_filter):
			print line

def draw_vertices_num(stat, ranks):
	for rank in ranks:
		def rank_filter(row):
			return row[rank_id] == rank

		fig, axes = plt.subplots(3, 1, sharex = True) 
		fig.figsize = (6, 2)
		plt.xlabel("|e|")

		ax = axes[0]
		ax.set_ylabel("vertices num")
		# ax.set_xlabel("|e|")
		stat.plot_aggregate(ax, size_id, stat.indices['vertices_num'], numpy.median, 'composition median', 
			rank_filter)
		stat.scatter(ax, size_id, stat.indices['vertices_num'], 'composition', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])

		ax.legend( loc =' best')

		ax = axes[1]
		ax.set_ylabel("vertices num")
		stat.plot_aggregate(ax, size_id, stat.indices['free_red_vn'], numpy.median, 'free reduction median',  
			rank_filter)
		stat.scatter(ax, size_id, stat.indices['free_red_vn'], 'free reduction', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		ax = axes[2]
		ax.set_ylabel("vertices num")
		stat.plot_aggregate(ax, size_id, stat.indices['nf_vn'], numpy.median, 'norm form median',
			rank_filter)
		stat.scatter(ax, size_id, stat.indices['nf_vn'], 'normal form', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		plt.show()

		plt.savefig("vertices_num_size_plot_rank{0}.pdf".format(rank))
		plt.close()


def draw_height(stat, ranks):
	for rank in ranks:
		def rank_filter(row):
			return row[rank_id] == rank

		fig, axes = plt.subplots(3, 1, sharex = True) 
		fig.figsize = (6, 2)
		plt.xlabel("|e|")

		ax = axes[0]
		ax.set_ylabel("height")
		# ax.set_xlabel("|e|")
		stat.plot_aggregate(ax, size_id, stat.indices['height'], numpy.median, 'composition median', 
			rank_filter)
		stat.scatter(ax, size_id, stat.indices['height'], 'composition', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		ax = axes[1]
		ax.set_ylabel("height")
		stat.plot_aggregate(ax, size_id, stat.indices['free_red_h'], numpy.median, 'free reduction median',  
			rank_filter)
		stat.scatter(ax, size_id, stat.indices['free_red_h'], 'free reduction', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		ax = axes[2]
		ax.set_ylabel("height")
		stat.plot_aggregate(ax, size_id, stat.indices['nf_ht'], numpy.median, 'norm form median',
			rank_filter)
		stat.scatter(ax, size_id, stat.indices['nf_ht'], 'normal form', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		plt.show()

		plt.savefig("height_size_plot_rank{0}.pdf".format(rank))
		plt.close()


def draw_time_on_vertices_num(stat, ranks):
	for rank in ranks:
		def rank_filter(row):
			return row[rank_id] == rank

		vn_id = stat.indices['vertices_num']

		fig, axes = plt.subplots(3, 1, sharex = True) 
		fig.figsize = (6, 2)
		plt.xlabel("vertices num")

		ax = axes[0]
		ax.set_ylabel("time")
		# ax.set_xlabel("|e|")
		# stat.plot_aggregate(ax, vn_id, stat.indices['vertices_num'], numpy.median, 'composition median', 
			# rank_filter)
		stat.scatter(ax, vn_id, stat.indices['vertices_num'], 'composition', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		ax = axes[1]
		ax.set_ylabel("time")
		# stat.plot_aggregate(ax, vn_id, stat.indices['free_red_time'], numpy.median, 'free reduction median',  
			# rank_filter)
		stat.scatter(ax, vn_id, stat.indices['free_red_time'], 'free reduction', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		ax = axes[2]
		ax.set_ylabel("time")
		# stat.plot_aggregate(ax, vn_id, stat.indices['nf_time'], numpy.median, 'norm form median',
			# rank_filter)
		stat.scatter(ax, vn_id, stat.indices['nf_time'], 'normal form', 
			rank_filter)
		ax.set_xlim([0, ax.get_xlim()[1]])
		ax.set_ylim([0, ax.get_ylim()[1]])
		ax.legend( loc =' best')

		plt.show()

		plt.savefig("time_vertices_num_plot_rank{0}.pdf".format(rank))
		plt.close()


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
	# size_id = stat.indices['|e|']
	# rank_id = stat.indices['rank']


	# ranks = {row[rank_id] for row in stat.rows}
	# draw_vertices_num(stat, ranks)
	# draw_time_on_vertices_num(stat, ranks)
	# draw_height(stat,ranks)
	

	




		


