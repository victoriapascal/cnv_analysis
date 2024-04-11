import matplotlib.pyplot as plt
import numpy as np
import json
import statsmodels.api as sm
import matplotlib.patches as mpatches
import os
import sys

def alus_to_discard(txt):
	remove = {}
	with open(txt, 'r') as f:
		for line in f:
			l = line.strip().split('_')
			ch = l[0]
			if ch in remove.keys():
				remove[ch].append(l[1])
			else:
				remove[ch] = [l[1]]

	return remove


def get_met_counts(inf, remove, outf):
	'''
	Using the methilation files as inut, assess the methylation at CpG and ALU level
	'''
	print('Assessing methylation')
	total_reads = 0
	cpgs_counts = {}
	with open(inf, 'r') as f:
		for line in f:
			l = line.strip().split('\t')
			com = l[8] ## the different compartments
			
			ch = l[0]
			cpg = l[3]
			met = l[5]
			m = int(l[6])
			u = int(l[7])
			
			if com == 'Alu':
				total_reads += m + u
				if len(l) == 10:
					alu = l[9]
					if not alu in remove[ch]:
						if ch in cpgs_counts.keys():
							if alu in cpgs_counts[ch].keys():
								rec = [alu, met, m, u]
								cpgs_counts[ch][alu].append(rec)

							else:
								rec = [alu, met, m, u]
								cpgs_counts[ch][alu] = [rec]

						else:
							cpgs_counts[ch] = {}
							rec = [alu, met, m, u]
							cpgs_counts[ch][alu] = [rec]

	
	alu_mean = {}
	alu_anot = {}
	for ch in cpgs_counts.keys():
		alu_mean[ch] = {}
		alu_anot[ch] = {}
		for k, v in cpgs_counts[ch].items():
			
			alu_anot[ch][k] = []
			for a in v:
				m = int(a[2])
				u = int(a[3])
				met = int(a[1])
				div = float(m/(m+u))
				reads = (m, u)
				alu_anot[ch][k].append(reads)
				if k in alu_mean[ch].keys():
					sum1 = alu_mean[ch][k][0] + m
					sum2 = alu_mean[ch][k][1] + u
					alu_mean[ch][k] = [sum1, sum2]
				else:
					alu_mean[ch][k] = [m, u]


	with open(outf, 'w') as out:
		out.write('Ch\tAlu\tlength\tNum. CpGs\tMet\tTotal M\t Total U \tCpGs (m,u)\n')
		for ch in alu_mean.keys():
			for alu in alu_mean[ch].keys():
				coords = alu.split(':')[1].split(' ')[0].split('-')
				lon = int(coords[1]) - int(coords[0])
				cpgs = '\t'.join(map(str, alu_anot[ch][alu]))
				mean = float(alu_mean[ch][alu][0]/(alu_mean[ch][alu][0] + alu_mean[ch][alu][1]))
				out.write(ch + '\t'+ alu + '\t'+ str(lon) + '\t'+ str(len(alu_anot[ch][alu])) + '\t'+ str(mean) + '\t' + str(alu_mean[ch][alu][0]) + '\t'  + str(alu_mean[ch][alu][1]) + '\t' + str(cpgs) + '\n')
	
	return alu_mean, total_reads

def get_bins(counts, chs, threshold, out):
	print('Getting the counts binned')
	counts_bins = {}
	for ch in chs:
		counts_bins[ch] = {}
		for i in range(threshold, chs[ch] + threshold, threshold):
			# counts_bins[ch][i] = [0, 0]
			for k in counts[ch].keys():
				end = int(k.split(':')[1].split(' ')[0].split('-')[1])
				# if end >= i:
				# 	counts_bins[ch][i + threshold] = [counts[ch][k][0], counts[ch][k][1]]
				# 	# print(i, end, counts[ch][k], counts_bins[ch][i])
				# 	break

				# else:
				if i == threshold:
					if 0 < end <  i: ######TO CHECK
						if i in counts_bins[ch].keys():
							sum1 = counts_bins[ch][i][0] + counts[ch][k][0]
							sum2 = counts_bins[ch][i][1] + counts[ch][k][1]
							counts_bins[ch][i] = [sum1, sum2]
							# print(i, end, counts[ch][k], counts_bins[ch][i])
						
						else:
							counts_bins[ch][i] = [counts[ch][k][0], counts[ch][k][1]]
						# print(i, end, counts[ch][k], counts_bins[ch][i])

				else:
					if (i - threshold) < end <  i: ## to not append values from the very beginning
						if i in counts_bins[ch].keys():
							sum1 = counts_bins[ch][i][0] + counts[ch][k][0]
							sum2 = counts_bins[ch][i][1] + counts[ch][k][1]
							counts_bins[ch][i] = [sum1, sum2]
							# print(i, end, counts[ch][k], counts_bins[ch][i])
						
						else:
							counts_bins[ch][i] = [counts[ch][k][0], counts[ch][k][1]]
								# print(i, end, counts[ch][k], counts_bins[ch][i])


	with open(out, 'w') as f1:
		for a in counts_bins.keys():
			for b in counts_bins[a].keys():
				f1.write(str(a) + '\t' + str(b) + '\t' + str(counts_bins[a][b][0]) + '\t' + str(counts_bins[a][b][1]) + '\n')

	return counts_bins

def normalize_bins_total_counts(counts_bins, total_reads):
	print('Normalizing the counts')
	counts_bin_norm = {}

	ratio = float(800000/total_reads)
	for ch in counts_bins.keys():
		counts_bin_norm[ch] = {}
		for a in counts_bins[ch]:
			v1 = counts_bins[ch][a][0] * ratio
			v2 = counts_bins[ch][a][1] * ratio
			counts_bin_norm[ch][a] = [v1, v2]

	return counts_bin_norm

def compute_mean(norm1, chromosome_lengths, threshold, outf):
	norm_means = {}
	with open(outf, 'w') as out:
		for ch in chromosome_lengths.keys():
			norm_means[ch] = {}
			x_axis = [str(i) for i in range(threshold, chromosome_lengths[ch] + threshold, threshold)]
			for x in x_axis:
				norm_means[ch][x] = None
				# print(ch, x, norm1[ch][x])
				if x in norm1[ch].keys():
					mean = float(norm1[ch][x][0]/ (norm1[ch][x][0] + norm1[ch][x][1]))
					norm_means[ch][x] = mean
					out.write(str(ch) + '\t' + str(x) + '\t' + str(mean) + '\n')


	return norm_means


def plot_scatter_chs(bins1, bins2, chromosome_lengths, name1, name2):
	print('Plotting results')
	num_chromosomes = len(chromosome_lengths.keys())
	num_cols = 5
	num_rows = -(-num_chromosomes // num_cols)  # Ceiling division to ensure enough rows

	# Create a grid of subplots
	fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 20))

	# Plot copy number segments for each chromosome
	for i, chromosome in enumerate(bins1.keys()):
		row_index = i // num_cols
		col_index = i % num_cols
		ax = axs[row_index, col_index] if num_rows > 1 else axs[col_index]
		
		# x_axis_all = [int(x.split(':')[1].split(' ')[0].split('-')[1]) for x in c1[chromosome]]
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

		# mean_all = [float(c1[chromosome][x][0]/ (c1[chromosome][x][0] + c1[chromosome][x][1])) for x in c1[chromosome]]

		# x_axis2_all = [int(x.split(':')[1].split(' ')[0].split('-')[1]) for x in c2[chromosome]] ##start coordinates leads
		x_axis2 = [int(i) for i in bins2[chromosome]]
		mean2 = [bins2[chromosome][i] for i in bins2[chromosome]]
		mean2_filtered = [0.8 if value is None else value for value in mean2]
		# mean2_all = [float(c2[chromosome][x][0]/ (c2[chromosome][x][0] + c2[chromosome][x][1])) for x in c2[chromosome]]

		if mean2[0] == None:
			first_non_none_index2 = next((i for i, x in enumerate(mean2) if x is not None), None)
			keys = list(bins1[chromosome].keys())[first_non_none_index2:]
			x_axis2 = [int(i) for i in keys]
			result2 = mean2[first_non_none_index2:] if first_non_none_index2 is not None else []
			mean2_filtered = [0.8 if value is None else value for value in result2]
		

		lowess1 = sm.nonparametric.lowess(mean_filtered, x_axis, frac = 0.1) #aprox. 10% of the data points are included in the local regression
		lowess2 = sm.nonparametric.lowess(mean2_filtered, x_axis2, frac = 0.1)

		# Calculate standard errors
		# smoothed_values1 = lowess1[:, 1]
		# residuals1 = mean - np.interp(x_axis, lowess1[:, 0], smoothed_values1)
		# standard_errors1 = np.std(residuals1)

		# smoothed_values2 = lowess2[:, 1]
		# residuals2 = mean2 - np.interp(x_axis2, lowess2[:, 0], smoothed_values2)
		# standard_errors2 = np.std(residuals2)


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
		# ax.fill_between(lowess1[:, 0], lowess1[:, 1] - shading_multiplier * standard_errors1, lowess1[:, 1] + shading_multiplier * standard_errors1, color='#FF204E', alpha=0.2, label='95% CI')
		# ax.fill_between(lowess2[:, 0], lowess2[:, 1] - shading_multiplier * standard_errors2, lowess2[:, 1] + shading_multiplier * standard_errors2, color='#135D66', alpha=0.2, label='95% CI')

		# Add labels and title
		ax.set_ylim((0,1.1))
		ax.set_xlim((0,chromosome_lengths[chromosome]))
		ax.set_ylabel('Degree of methylation')
		ax.set_title(f"Chromosome {str(chromosome)}")
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
	fig.legend(colors, labels, loc='lower right',bbox_to_anchor=(0.95, 0.1), shadow=True, ncol=2)
	
	# Save or show the plot
	plt.savefig('MET_plots_sample_' + name1 + '_' + name2 + '.png')

	# plt.show()

if __name__ == "__main__":
	
	##set bin thresholds
	threshold = 1000000

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
	
	alu = os.listdir(sys.argv[1]) 
	pair = sys.argv[2] #1-2

	d1 = '{:02d}'.format(int(pair.split('-')[0]))
	d2 = '{:02d}'.format(int(pair.split('-')[1]))

	name1 = alu[0].split('_')[0] + '_' + str(d1)
	name2 = alu[0].split('_')[0] + '_' + str(d2)

	for a in alu:
		if a.startswith(name1):
			f = sys.argv[1] + os.sep + a	
		if a.startswith(name2):
			f2 = sys.argv[1] + os.sep + a

	# name1 = 'run75_01'
	# name2 = 'run75_02'

	# f = 'run75_01_001_GGTAGC.RGok.MAPQ.sorted.markdupAniling.merged.mapped'
	# f2 = 'run75_02_001_ATGAGC.RGok.MAPQ.sorted.markdupAniling.merged.mapped'
	
	outf =  name1 + '_meth_per_alu.tsv'
	outf2 =  name2 + '_meth_per_alu.tsv'

	remove = alus_to_discard('Alus_to_discard.txt')
	
	counts1, total_reads1 = get_met_counts(f, remove, outf)
	counts2, total_reads2 = get_met_counts(f2, remove, outf2)
	
	outb1 =  name1 + '_meth_per_alu_bins.tsv'
	outb2 =  name2 + '_meth_per_alu_bins.tsv'

	
	bins1 = get_bins(counts1, chromosome_lengths, threshold, outb1)
	bins2 = get_bins(counts2, chromosome_lengths, threshold, outb2)


	bins1 = json.load(open('counts_' + name1 + '.txt'))
	bins2 = json.load(open('counts_' + name2 + '.txt'))

	total_reads1 = 1528645
	total_reads2 = 1118578

	
	norm1 = normalize_bins_total_counts(bins1, total_reads1)
	norm2 = normalize_bins_total_counts(bins2, total_reads2)

	outm1 = name1 + '_meth_per_alu_bins_norm.tsv'
	outm2 = name2 + '_meth_per_alu_bins_norm.tsv'
	
	mean1 = compute_mean(norm1, chromosome_lengths, threshold, outm1)
	mean2 = compute_mean(norm2, chromosome_lengths, threshold, outm2)

	json.dump(mean1, open('counts_' + name1 + '.txt','w'))
	json.dump(mean2, open('counts_' + name2 +'.txt','w'))


	# mean1 = json.load(open('counts_' + name1 + '.txt'))
	# mean2 = json.load(open('counts_' + name2 + '.txt'))


	plot_scatter_chs(mean1, mean2, chromosome_lengths, name1, name2)


