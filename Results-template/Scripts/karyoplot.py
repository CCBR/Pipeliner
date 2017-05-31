import os
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import math
import numpy
plt.switch_backend('agg')

def karyoplot(karyo_filename, metadata={}, part=1):
	'''
	To create a karyo_filename go to: http://genome.ucsc.edu/cgi-bin/hgTables 
	group: Mapping and Sequencing
	track: Chromosome Band 
	An example of an output (hg19, Human) is here: http://pastebin.com/6nBX6sdE 
	The script will plot dots next to loci defined in metadata as:
	metadata = {
		'1' : [2300000, 125000000, 249250621],
	}
	'''


	karyo_dict={}
	with open(karyo_filename) as karyo_f:
		lines = [x.replace(os.linesep, '').split() for x in karyo_f.readlines()]
		for chromosome in [str(x) for x in range(1,23)] + ['X', 'Y']:
			karyo_dict[chromosome] = [[y[0], int(y[1]), int(y[2]), y[3], y[4]] for y in [x for x in lines if x[0] == 'chr' + chromosome]]
	scores=[]
	for c in karyo_dict.values():
		for y in c:
			scores.append(float(y[4]))
	maxscore=max([abs(min(scores)),abs(max(scores))])
	minscore=-1*maxscore
	print(minscore,maxscore)
	
	fig, ax = plt.subplots()

	DIM = 1.0

	ax.set_xlim([0.0, DIM * (1.3)])
	ax.set_ylim([0.0, DIM])

	def get_chromosome_length(chromosome):
		chromosome_start = float(min([x[1] for x in karyo_dict[chromosome]]))
		chromosome_end = float(max(x[2] for x in karyo_dict[chromosome]))
		chromosome_length = chromosome_end - chromosome_start

		#return int(chromosome_length*1.05)
		return chromosome_length

	# def convert_to_rgb(minval, maxval, val, colors):
	# 	width=maxval-minval
	# 	val-=minval # shift origin to minval
	# 	val/=width # change to a 0-100 scale
	# 	val=int(val*100)
	# 	if val > 50:
	# 		r=1-2*(val-50)/100.0
	# 		g=1.0
	# 	else:
	# 		r=1.0
	# 		g=2*val/100.0
	# 	b=0.0
	# 	if r<0:
	# 		r=0.0
	# 	if g<0:
	# 		g=0.0
	# 	if r>1.0:
	# 		r=1.0
	# 	if g>1.0:
	# 		g=1.0
	# 	return r,g,b

	def convert_to_rgb(minval,maxval,val,colors):
		step=(maxval*1.05-minval*1.05)/13
		bins=numpy.arange(minval*1.05,maxval*1.05+step,step)
		return colors[numpy.digitize(numpy.array([val]),bins=bins)[0]]

	def plot_chromosome(chromosome, order):

		chromosome_length = get_chromosome_length(chromosome)
		chromosome_length_1 = get_chromosome_length('1')

		x_start = order * DIM * 0.1 
		x_end = x_start + (DIM * 0.04)
		y_start = DIM * 0.8 * (chromosome_length/chromosome_length_1)
		y_end = DIM * 0.1


		# We use the same colors as: http://circos.ca/tutorials/lessons/2d_tracks/connectors/configuration 
		colors = {
			1 : (0/255.0,102/255.0,0/255.0),
			2 : (0/255.0,204/255.0,0/255.0),
			3 : (51/255.0,204/255.0,51/255.0),
			4 : (102/255.0,255/255.0,51/255.0),
			5 : (153/255.0,255/255.0,51/255.0),
			6 : (204/255.0,255/255.0,51/255.0),
			7 : (255/255.0,255/255.0,0/255.0),
			8 : (255/255.0,204/255.0,0/255.0),
			9 : (255/255.0,153/255.0,51/255.0),
			10 : (255/255.0,102/255.0,51/255.0),
			11 : (204/255.0,51/255.0,51/255.0),
			12 : (204/255.0,0/255.0,0/255.0),
			13 : (255/255.0,0/255.0,0/255.0),
		}

		for index, piece in enumerate(karyo_dict[chromosome]):

			current_height = piece[2] - piece[1]
			current_height_sc = ((y_end - y_start) / chromosome_length) * current_height
			if index == 0:
				y_previous = y_start

			y_next = y_previous + current_height_sc

			#color = colors[piece[4]]
			color = convert_to_rgb(minscore,maxscore,float(piece[4]), colors)  # [BLUE, GREEN, RED]
			#print(piece[4],color)

			#plot the caryotypes
			r = Rectangle((x_start, y_previous), x_end-x_start, current_height_sc, color = color)
			ax.add_patch(r)

			y_previous = y_next

		#Plot semicircles at the beginning and end of the chromosomes
		center_x = x_start + (x_end-x_start)/2.0
		radius = (x_end-x_start)/2.0
		theta1 = 0.0
		theta2 = 180.0
		w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
		w2 = Wedge((center_x, y_end), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
		ax.add_patch(w1)
		ax.add_patch(w2)
		ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
		ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')

		#Plot metadata
		if chromosome in metadata:
			for md in metadata[chromosome]:
				ax.plot([x_end + (DIM*0.015)], [y_start + (y_end-y_start) * (md/chromosome_length)], '.', color='black')

		ax.text(x_start, y_end - (DIM * 0.07), chromosome)

	if part==1:
		plot_chromosome('1', 1)
		plot_chromosome('2', 2)
		plot_chromosome('3', 3)
		plot_chromosome('4', 4)
		plot_chromosome('5', 5)
		plot_chromosome('6', 6)
		plot_chromosome('7', 7)
		plot_chromosome('8', 8)
		plot_chromosome('9', 9)
		plot_chromosome('10', 10)
		plot_chromosome('11', 11)
		plot_chromosome('12', 12)


	elif part==2:
		plot_chromosome('13', 1)
		plot_chromosome('14', 2)
		plot_chromosome('15', 3)
		plot_chromosome('16', 4)
		plot_chromosome('17', 5)
		plot_chromosome('18', 6)
		plot_chromosome('19', 7)
		
		if "hg" in species:
			plot_chromosome('20', 8)
			plot_chromosome('21', 9)
			plot_chromosome('22', 10)

		plot_chromosome('X', 11)
		plot_chromosome('Y', 12)

	elif part==0:
		#print("HERE1")
		plot_chromosome('1', 1)
		#print("HERE3")
	
	else:
		raise Exception('plot argument should be either "1" or "2"')

	plt.axis('off')
	return

import sys
fn = sys.argv[1]
species = sys.argv[2]
#print('plotting..')
karyoplot(fn, part=1)
#print('HERE4')
plt.savefig(fn+'1.png')
karyoplot(fn, part=2)
#print('HERE4')
plt.savefig(fn+'2.png')
#print('HERE5')
#karyoplot(fn, part=1)
#karyoplot(fn, part=2)
