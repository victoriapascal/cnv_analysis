import json
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
import statsmodels.api as sm
import matplotlib.patches as mpatches

def plot_scatter_chs(bins1, bins2, chromosome_lengths, output_dir, ch, c):
	# print('Plotting results')
	num_chromosomes = len(bins1.keys())
	num_cols = 5
	num_rows = -(-num_chromosomes // num_cols)  # Ceiling division to ensure enough rows

	# Create a grid of subplots
	fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 20))

	# Plot copy number segments for each chromosome
	for i, chromosome in enumerate(bins1.keys()):
		row_index = i // num_cols
		col_index = i % num_cols

		ax = axs[row_index, col_index] if num_rows > 1 else axs[col_index]
		

		x_axis = [int(i) for i in bins1[chromosome]]
		mean = [bins1[chromosome][i] for i in bins1[chromosome]]
		mean_filtered = [0.8 if value is None else value for value in mean]

		result = None
		if mean[0] == None:
			first_non_none_index = next((i for i, x in enumerate(mean) if x is not None), None)
			keys = list(bins1[chromosome].keys())[first_non_none_index:]
			x_axis = [int(i) for i in keys]
			result = mean[first_non_none_index:] if first_non_none_index is not None else []
			mean_filtered = [0.8 if value is None else value for value in result]

		x_axis2 = [int(i) for i in bins2[chromosome]]
		mean2 = [bins2[chromosome][i] for i in bins2[chromosome]]
		mean2_filtered = [0.8 if value is None else value for value in mean2]

		if mean2[0] == None:
			first_non_none_index2 = next((i for i, x in enumerate(mean2) if x is not None), None)
			keys = list(bins1[chromosome].keys())[first_non_none_index2:]
			x_axis2 = [int(i) for i in keys]
			result2 = mean2[first_non_none_index2:] if first_non_none_index2 is not None else []
			mean2_filtered = [0.8 if value is None else value for value in result2]
		

		lowess1 = sm.nonparametric.lowess(mean_filtered, x_axis, frac = 0.1) #aprox. 10% of the data points are included in the local regression
		lowess2 = sm.nonparametric.lowess(mean2_filtered, x_axis2, frac = 0.1)

		# Plot the first scatter plot
		if result is not None: ## check if the variable results exist
			ax.scatter(x_axis, result, color='#FF204E', label='Tumor', s = 0.5)
			ax.scatter(x_axis2, result2, color='#135D66', label='Normal', s = 0.5)
		else:
			ax.scatter(x_axis, mean, color='#FF204E', label='Tumor', s = 0.5)
			ax.scatter(x_axis2, mean2, color='#135D66', label='Normal', s = 0.5)

		##add lowess lines
		shading_multiplier = 0.2
		ax.plot(lowess1[:, 0], lowess1[:, 1], color='#FF204E', linewidth=1)
		ax.plot(lowess2[:, 0], lowess2[:, 1], color='#135D66', linewidth=1)
		
		# Add labels and title
		ax.set_ylim((0,1.1))
		ax.set_xlim((0,chromosome_lengths[ch]))
		ax.set_ylabel('Degree of methylation')
		ax.set_title(f"Sample {str(chromosome.replace('counts_', '').replace('.txt', ''))}")
		ax.grid(True)


	# Hide empty subplots
	for i in range(num_chromosomes, num_rows * num_cols):
		row_index = i // num_cols
		col_index = i % num_cols
		axs[row_index, col_index].axis('off')

	# Adjust layout
	plt.tight_layout()

	##include legend
	labels = ['Tumor', 'Normal']
	colors = []
	colors.append(mpatches.Patch(color='#FF204E'))
	colors.append(mpatches.Patch(color='#135D66'))
	# fig.legend(colors, labels, loc='lower right',bbox_to_anchor=(0.95, 0.1), shadow=True, ncol=2)
	
	# Save or show the plot
	plt.savefig( output_dir + 'MET_plots_chr' + str(ch) + '_bins_norm_800K_' + str(c) + '.png') ##variable for the number of chs
	# plt.savefig('MET_plots_sample_' + name1 + '_' + name2 + '.png')

	# plt.show()

if __name__ == "__main__":

	##chromosome length dict
	chromosome_lengths = {
	    '1': 249237862,
	    '2': 243188438,
	    '3': 197959720,
	    '4': 191039922,
	    '5': 180903375,
	    '6': 171046147,
	    '7': 159121704,
	    '8': 146298074,
	    '9': 141149219,
	    '10': 135520070,
	    '11': 134941579,
	    '12': 133839031,
	    '13': 115108862,
	    '14': 107287443,
	    '15': 102491545,
	    '16': 90282411,
	    '17': 81192741,
	    '18': 78009649,
	    '19': 59114562,
	    '20': 62963629,
	    '21': 48117466,
	    '22': 51244339,
	    'X': 155246021,
	    'Y': 59033256}

	# chs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
	iteration = [['counts_run75_01.txt,counts_run75_02.txt', 'counts_run75_04.txt,counts_run75_03.txt', 'counts_run75_06.txt,counts_run75_05.txt', 'counts_run75_08.txt,counts_run75_09.txt', 'counts_run75_10.txt,counts_run75_11.txt', 'counts_run75_12.txt,counts_run75_13.txt', 'counts_run77_02.txt,counts_run77_01.txt', 'counts_run77_04.txt,counts_run77_03.txt', 'counts_run77_06.txt,counts_run77_05.txt', 'counts_run77_08.txt,counts_run77_09.txt', 'counts_run77_10.txt,counts_run77_11.txt', 'counts_run77_12.txt,counts_run77_13.txt', 'counts_run80_01.txt,counts_run80_02.txt', 'counts_run80_03.txt,counts_run80_04.txt', 'counts_run80_05.txt,counts_run80_06.txt', 'counts_run80_08.txt,counts_run80_09.txt', 'counts_run80_10.txt,counts_run80_11.txt', 'counts_run80_12.txt,counts_run80_13.txt', 'counts_Run86_08.txt,counts_Run86_09.txt', 'counts_Run86_10.txt,counts_Run86_11.txt', 'counts_Run86_12.txt,counts_Run86_13.txt', 'counts_run89_01.txt,counts_run89_02.txt', 'counts_run89_03.txt,counts_run89_04.txt', 'counts_run89_05.txt,counts_run89_06.txt','counts_run89_08.txt,counts_run89_09.txt'],['counts_run89_10.txt,counts_run89_11.txt', 'counts_run89_12.txt,counts_run89_13.txt', 'counts_run92_01.txt,counts_run92_02.txt', 'counts_run92_03.txt,counts_run92_04.txt', 'counts_run92_05.txt,counts_run92_06.txt', 'counts_run92_08.txt,counts_run92_09.txt', 'counts_run92_10.txt,counts_run92_11.txt', 'counts_run92_12.txt,counts_run92_13.txt', 'counts_run94_01.txt,counts_run94_02.txt', 'counts_run94_03.txt,counts_run94_04.txt', 'counts_run94_05.txt,counts_run94_06.txt', 'counts_run94_08.txt,counts_run94_09.txt', 'counts_run94_10.txt,counts_run94_11.txt', 'counts_run94_12.txt,counts_run94_13.txt', 'counts_run97_01.txt,counts_run97_02.txt', 'counts_run97_03.txt,counts_run97_04.txt', 'counts_run97_05.txt,counts_run97_06.txt', 'counts_run97_08.txt,counts_run97_09.txt', 'counts_run97_10.txt,counts_run97_11.txt', 'counts_run97_12.txt,counts_run97_13.txt', 'counts_run99_01.txt,counts_run99_02.txt', 'counts_run99_03.txt,counts_run99_04.txt', 'counts_run99_05.txt,counts_run99_06.txt', 'counts_run99_08.txt,counts_run99_09.txt', 'counts_run99_10.txt,counts_run99_11.txt']] 

	# iteration = [['ALU010_segments_sample_1-2.txt', 'ALU010_segments_sample_4-3.txt', 'ALU010_segments_sample_6-5.txt', 'ALU010_segments_sample_8-9.txt', 'ALU010_segments_sample_10-11.txt', 'ALU010_segments_sample_12-13.txt', 'ALU011_segments_sample_2-1.txt', 'ALU011_segments_sample_4-3.txt', 'ALU011_segments_sample_6-5.txt', 'ALU011_segments_sample_8-9.txt', 'ALU011_segments_sample_10-11.txt', 'ALU011_segments_sample_12-13.txt', 'ALU012_segments_sample_1-2.txt', 'ALU012_segments_sample_3-4.txt', 'ALU012_segments_sample_5-6.txt', 'ALU012_segments_sample_8-9.txt', 'ALU012_segments_sample_10-11.txt', 'ALU012_segments_sample_12-13.txt', 'ALU013_segments_sample_8-9.txt', 'ALU013_segments_sample_10-11.txt', 'ALU013_segments_sample_12-13.txt', 'ALU014_segments_sample_1-2.txt', 'ALU014_segments_sample_3-4.txt', 'ALU014_segments_sample_5-6.txt', 'ALU014_segments_sample_8-9.txt'], ['ALU014_segments_sample_10-11.txt', 'ALU014_segments_sample_12-13.txt', 'ALU015_segments_sample_1-2.txt', 'ALU015_segments_sample_3-4.txt', 'ALU015_segments_sample_5-6.txt', 'ALU015_segments_sample_8-9.txt', 'ALU015_segments_sample_10-11.txt', 'ALU015_segments_sample_12-13.txt', 'ALU016_segments_sample_1-2.txt', 'ALU016_segments_sample_3-4.txt', 'ALU016_segments_sample_5-6.txt', 'ALU016_segments_sample_8-9.txt', 'ALU016_segments_sample_10-11.txt', 'ALU016_segments_sample_12-13.txt', 'ALU017_segments_sample_1-2.txt', 'ALU017_segments_sample_3-4.txt', 'ALU017_segments_sample_5-6.txt', 'ALU017_segments_sample_8-9.txt', 'ALU017_segments_sample_10-11.txt', 'ALU017_segments_sample_12-13.txt', 'ALU018_segments_sample_1-2.txt', 'ALU018_segments_sample_3-4.txt', 'ALU018_segments_sample_5-6.txt', 'ALU018_segments_sample_8-9.txt', 'ALU018_segments_sample_10-11.txt']]
	output_dir = '/Users/Victoria/Downloads/'

	for ch in chromosome_lengths.keys():

		for i, files in zip(['1','2'], iteration):
			print(ch, i)
			sample1_chs = {}
			sample2_chs = {}
			for file in files:
				# d = json.load(open('segments+log2/' + file.replace('segments', 'log2')))
				p1 = json.load(open('counts_bins_norm_mean/' + file.split(',')[0]))
				p2 = json.load(open('counts_bins_norm_mean/' + file.split(',')[1]))

				sample1_chs[file.split(',')[0]] = p1[ch]
				sample2_chs[file.split(',')[0]] = p2[ch]
			# chs)
			plot_scatter_chs(sample1_chs, sample2_chs, chromosome_lengths, output_dir, ch, i)


