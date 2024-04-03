import json
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cbs
import os
import pysam
import sys


def plot_bins(chr_log2, threshold, segment_regions, output_dir, c, ch):

	# Calculate number of rows and columns for subplot grid
	num_chromosomes = len(chr_log2.keys())
	num_cols = 5
	num_rows = -(-num_chromosomes // num_cols)  # Ceiling division to ensure enough rows

	# Create a grid of subplots
	fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 20))

	# Plot copy number segments for each chromosome
	for i, chromosome in enumerate(chr_log2.keys()):
		row_index = i // num_cols
		col_index = i % num_cols
		ax = axs[row_index, col_index] if num_rows > 1 else axs[col_index]
	    # chromosome_data = cnr_data[cnr_data['chromosome'] == chromosome]
	    # plot_copy_number_segments(ax, chromosome_data, chromosome)
		x_axis = [i*threshold for i in range(len(chr_log2[chromosome]))]
		ax.scatter(x_axis, chr_log2[chromosome], s=10)
		ax.set_ylim((-2,2))
		ax.set_title(f"sample {str(chromosome.replace('_log2_sample', '').replace('.txt', ''))}")
		ax.set_xlabel('Genomic Position')
		ax.set_ylabel('Log2 Copy Number Ratio')
		ax.grid(True)
		# ax.axhline(y = 0,  color='orange', linestyle='-', linewidth=2) ## mean line at 0
		for l in segment_regions[chromosome.replace('log2', 'segments')]:
			ax.hlines(l[2], xmin=l[0], xmax=l[1], color='orange', linestyle='-', linewidth=2)


	# Hide empty subplots
	for i in range(num_chromosomes, num_rows * num_cols):
		row_index = i // num_cols
		col_index = i % num_cols
		axs[row_index, col_index].axis('off')

	# Adjust layout
	plt.tight_layout()

	# Save or show the plot
	plt.savefig( output_dir + 'CNV_plots_chr_' + str(ch) + '_bins_norm_800K_segements' + str(c) + '.png') ##variable for the number of chs
	# plt.show()

if __name__ == "__main__":

	chs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

	iteration = [['ALU010_segments_sample_1-2.txt', 'ALU010_segments_sample_4-3.txt', 'ALU010_segments_sample_6-5.txt', 'ALU010_segments_sample_8-9.txt', 'ALU010_segments_sample_10-11.txt', 'ALU010_segments_sample_12-13.txt', 'ALU011_segments_sample_2-1.txt', 'ALU011_segments_sample_4-3.txt', 'ALU011_segments_sample_6-5.txt', 'ALU011_segments_sample_8-9.txt', 'ALU011_segments_sample_10-11.txt', 'ALU011_segments_sample_12-13.txt', 'ALU012_segments_sample_1-2.txt', 'ALU012_segments_sample_3-4.txt', 'ALU012_segments_sample_5-6.txt', 'ALU012_segments_sample_8-9.txt', 'ALU012_segments_sample_10-11.txt', 'ALU012_segments_sample_12-13.txt', 'ALU013_segments_sample_8-9.txt', 'ALU013_segments_sample_10-11.txt', 'ALU013_segments_sample_12-13.txt', 'ALU014_segments_sample_1-2.txt', 'ALU014_segments_sample_3-4.txt', 'ALU014_segments_sample_5-6.txt', 'ALU014_segments_sample_8-9.txt'], ['ALU014_segments_sample_10-11.txt', 'ALU014_segments_sample_12-13.txt', 'ALU015_segments_sample_1-2.txt', 'ALU015_segments_sample_3-4.txt', 'ALU015_segments_sample_5-6.txt', 'ALU015_segments_sample_8-9.txt', 'ALU015_segments_sample_10-11.txt', 'ALU015_segments_sample_12-13.txt', 'ALU016_segments_sample_1-2.txt', 'ALU016_segments_sample_3-4.txt', 'ALU016_segments_sample_5-6.txt', 'ALU016_segments_sample_8-9.txt', 'ALU016_segments_sample_10-11.txt', 'ALU016_segments_sample_12-13.txt', 'ALU017_segments_sample_1-2.txt', 'ALU017_segments_sample_3-4.txt', 'ALU017_segments_sample_5-6.txt', 'ALU017_segments_sample_8-9.txt', 'ALU017_segments_sample_10-11.txt', 'ALU017_segments_sample_12-13.txt', 'ALU018_segments_sample_1-2.txt', 'ALU018_segments_sample_3-4.txt', 'ALU018_segments_sample_5-6.txt', 'ALU018_segments_sample_8-9.txt', 'ALU018_segments_sample_10-11.txt']]
	output_dir = '/Users/Victoria/Downloads/'

	for ch in chs:

		for i, files in zip(['1','2'], iteration):
			sample_chs = {}
			sample_segments = {}
			for file in files:
				print(i)
				d = json.load(open('segments+log2/' + file.replace('segments', 'log2')))
				p = json.load(open('segments+log2/' + file))
				sample_chs[file.replace('segments', 'log2')] = d[ch]
				sample_segments[file] = p[ch]
			plot_bins(sample_chs, 1000000, sample_segments, output_dir, i, ch)


