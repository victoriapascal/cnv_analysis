import json
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cbs
import os
import pysam
import sys

def sort_tsv_file(tsv, out_tsv):
	df = pd.read_csv(tsv, sep='\t', skiprows=1, dtype={'Chr': str})
	
	# Filter out lines where the 'Chr' column contains an underscore ('_')
	df = df[~df['Chr'].str.contains('_')]

	# Define a custom sorting order for the 'Chr' column
	custom_order = {'X': 23, 'Y': 24}
	df['Chr'] = df['Chr'].replace(custom_order)
	df['Chr'] = pd.to_numeric(df['Chr'], errors='coerce')

	# Sort the DataFrame based on the 'Chr' column
	sorted_df = df.sort_values(by=['Chr', 'Start'])
	reverse_order = {23: 'X', 24: 'Y'}
	sorted_df['Chr'] = sorted_df['Chr'].replace(reverse_order)	
	sorted_df.to_csv(out_tsv, sep='\t', index=False)	

def bin_fc_file(chromosome_lengths, inf, threshold, pair, output_dir, run):
	with open(inf, 'r') as f:
		next(f)
		lines = f.readlines()
	
	##split each chromosome in regions of aprox. 1MB
	print('Binning the chromosomes')
	chr_bins_normal = {} ##break points of the bins
	chr_bins_tumor = {}
	chr_bin_alu_counts = {}

	sample1 = int(pair.split('-')[0]) + 5 
	sample2 = int(pair.split('-')[1]) + 5 
	for ch in chromosome_lengths.keys():
		print(ch)
		chr_bins_tumor[ch] = {}
		chr_bins_normal[ch] = {}
		chr_bin_alu_counts[ch] = {}


		for i in range(threshold, chromosome_lengths[ch]+ threshold, threshold): ##range over the chromosome lenght in steps of thr size
			# print(ch, i)
			chr_bins_tumor[ch][i] = 0
			chr_bins_normal[ch][i] = 0
			chr_bin_alu_counts[ch][i] = 0
			if ch == '1' and i == threshold:
				start_pos = 0
			for l in lines[start_pos:]:
				line = l.strip().split('\t')
				chs = line[1]
				#end = int(line[3])
				#reads_tumor = int(line[6])
				#reads_normal = int(line[7])
				if ch == chs:
					end = int(line[3])
					reads_tumor = int(line[sample1]) ## depending on the sample, line[6] for sample 1
					reads_normal = int(line[sample2]) ## depending on the sample, line[7] for sample 2

					if end >= i:
						# print(ch, end, i)
						start_pos = lines.index(l)
						break
						# if i in chr_bins_tumor[ch].keys():
						# 	chr_bins_tumor[ch][i] = chr_bins_tumor[ch][i] + reads_tumor
						# 	chr_bins_normal[ch][i] = chr_bins_normal[ch][i] + reads_normal
						# 	chr_bin_alu_counts[ch][i] = chr_bin_alu_counts[ch][i] + 1
						
						# else:
						# 	chr_bins_tumor[ch][i] = reads_tumor
						# 	chr_bins_normal[ch][i] = reads_normal
						# 	chr_bin_alu_counts[ch][i] = 1

						# start_pos = lines.index(l)
						
					
					else:
						if i in chr_bins_tumor[ch].keys():
							chr_bins_tumor[ch][i] = chr_bins_tumor[ch][i] + reads_tumor
							chr_bins_normal[ch][i] = chr_bins_normal[ch][i] + reads_normal
							chr_bin_alu_counts[ch][i] = chr_bin_alu_counts[ch][i] + 1
						else:
							chr_bins_tumor[ch][i] = reads_tumor
							chr_bins_normal[ch][i] = reads_normal
							chr_bin_alu_counts[ch][i] = 1


	total_reads_tumor = 0
	total_reads_normal = 0

	with open(output_dir + 'report_' + run + '_read_counts_bins_T-N' + pair + '.tsv', 'w') as out:
		out.write('Chr_bin\tReads_Tumor\tReads_Normal\tNum_Alus\n')
		for k in chr_bins_tumor.keys():
			for v in chr_bins_tumor[k].keys():
				total_reads_tumor = total_reads_tumor + chr_bins_tumor[k][v]
				total_reads_normal = total_reads_normal + chr_bins_normal[k][v]
				out.write(str(k) + '_' + str(v) + '\t' + str(chr_bins_tumor[k][v]) + '\t' + str(chr_bins_normal[k][v]) +  '\t' + str(chr_bin_alu_counts[k][v]) + '\n')
	
	return chr_bins_tumor, chr_bins_normal, total_reads_tumor, total_reads_normal

def get_total_read_counts(bam_file):
	bam = pysam.AlignmentFile(bam_file, "rb")
	read_count = 0

	# Initialize a counter for the number of reads
	read_count = 0

	# Iterate through each read in the BAM file and count them
	for read in bam.fetch():
	    read_count += 1

	# Close the BAM file
	bam.close()

	return read_count

def compute_log2_tumor_normal(chr_bins_tumor, chr_bins_normal, ratioT, ratioN, output_dir, run, pair):
	out = open(output_dir + 'log2_values_' + run + '_' + pair + '_for_bin_and_chs.tsv', 'w')
	out.write('Clone\tChromosome\tPosition\tlog2\n') ## for the DNAcopy package
	out2 = open(output_dir + 'log2_counts_'+ run + '_' + pair  +'_for_bin_and_chs.tsv', 'w')
	out2.write('Chromosome_bin\treads_tumor\treads_normal\tLog2(T) - Log2(N)\n')

	chr_log2 = {}
	for ch in chr_bins_tumor.keys():
		x_axis = [i for i in range(threshold, chromosome_lengths[ch] + threshold, threshold)]
		chr_log2[ch] = []
		for b in chr_bins_tumor[ch].keys():
			
			# if the Tumor does not have any read in this specific bin but normal does
			if chr_bins_tumor[ch][b] == 0 and chr_bins_normal[ch][b] > 0:
				diff = math.log2(0.00001) - math.log2(chr_bins_normal[ch][b]*ratioN)
				chr_log2[ch].append(2)
				out.write('bin_' + str(b) + '\t' + str(ch) + '\t' + str(b) + '\t' + '2' + '\n')
				out2.write(str(ch) + '_' + str(b) + '\t' + str(chr_bins_tumor[ch][b]*ratioT) + '\t' + str(chr_bins_normal[ch][b]*ratioN) + '\t' + '2' + '\n')
			
			# if the Normal does not have any read in this specific bin but Tumor does
			if chr_bins_normal[ch][b] == 0 and chr_bins_tumor[ch][b] > 0:
				diff = math.log2(chr_bins_tumor[ch][b]*ratioT) - math.log2(0.00001)
				chr_log2[ch].append(-2)
				out.write('bin_' + str(b) + '\t' + str(ch) + '\t' + str(b) + '\t' + '-2' + '\n')
				out2.write(str(ch) + '_' + str(b) + '\t' + str(chr_bins_tumor[ch][b]*ratioT) + '\t' + str(chr_bins_normal[ch][b]*ratioN) + '\t' + '-2' + '\n')
			
			# if none Tumor or Normal have reads 
			if chr_bins_tumor[ch][b] == 0 and chr_bins_normal[ch][b] == 0:
				chr_log2[ch].append(None)
				out2.write(str(ch) + '_' + str(b) + '\t' + str(chr_bins_tumor[ch][b]*ratioT) + '\t' + str(chr_bins_normal[ch][b]*ratioN) + '\t' + 'None' + '\n')

			# Both Tumor and normal have reads (most likely situation)
			if chr_bins_tumor[ch][b] > 0 and chr_bins_normal[ch][b] > 0:
				norm_n = float(chr_bins_normal[ch][b]*ratioN)
				diff = math.log2(chr_bins_tumor[ch][b]*ratioT) - math.log2(norm_n)
				chr_log2[ch].append(diff)
				# print(chr_bins_tumor[ch][b],norm_t, chr_bins_normal[ch][b], norm_n, diff)
				out.write('bin_' + str(b) + '\t' + str(ch) + '\t' + str(b) + '\t' + str(diff) + '\n')
				out2.write(str(ch) + '_' + str(b) + '\t' + str(chr_bins_tumor[ch][b]*ratioT) + '\t' + str(norm_n) + '\t' + str(diff) + '\n')

	out.close()
	out2.close()
	return chr_log2


def smooth_log2values_bin(chr_log2):
	"""Apply moving average smoothing to data."""
	
	value_conversion = {}
	chr_bin_counts_smoothed = {}
	for a in chr_log2.keys():
		value_conversion[a] = {}
		chr_bin_counts_smoothed[a] = []
		
		cpr_filtered = [a for a in chr_log2[a] if not a == None]
		window_size = 3
		mask = np.array([val is not None for val in chr_log2[a]])

		cpr = np.array(cpr_filtered)

		smoothed_data = np.convolve(cpr, np.ones(window_size)/window_size, mode='same')
		smoothed_list = smoothed_data.tolist()
		for k, v in zip(cpr_filtered, smoothed_list):
			value_conversion[a][k] = v


		for value in chr_log2[a]:
			if value is None:
				chr_bin_counts_smoothed[a].append(None)
			else:
				chr_bin_counts_smoothed[a].append(value_conversion[a][value])

	return chr_bin_counts_smoothed, value_conversion

def get_segments(chr_log2):
	chr_segments = {}
	for a in chr_log2.keys():
		values = [v for v in chr_log2[a] if v is not None] ## we compute the segments deicarting the None values
		L = cbs.segment(values)
		chr_segments[a] = L

	return chr_segments

def get_segment_regions(chromosome_lengths, chr_segments, chr_log2):
	"""Link the smoothed log2 values to the non-smoothed value so that we can get the actual segment length from the bin coordinates"""
	
	segment_regions = {}
	for ch in chr_segments.keys():
		segment_regions[ch] = []
		x_axis = [i for i in range(threshold, chromosome_lengths[ch] + threshold, threshold)]
		values = [v for v in chr_log2[ch] if v is not None]
		coords_end = []
		for v in chr_segments[ch]:
			v1 = values[v[0]]
			if v[1] == len(values):
				v2 = values[-1]
			else:
				v2 = values[v[1]]
			
			start_segment = chr_log2[ch].index(v1)
			end_segment = chr_log2[ch].index(v2)
			counts_start = chr_log2[ch].count(v1)
			counts_end = chr_log2[ch].count(v2)
			indices_s = [i for i in range(len(chr_log2[ch])) if chr_log2[ch][i] == v1]
			indices_e = [i for i in range(len(chr_log2[ch])) if chr_log2[ch][i] == v2]
			
			print(ch, len(chr_segments[ch]), start_segment, end_segment, counts_start, counts_end)
			if counts_start == 1 and counts_end == 1:
			# if counts_start == 1 and counts_end == 1 or not coords_end:
				log2_segment = chr_log2[ch][start_segment:end_segment]
				log2_segment_wo_nones = [a for a in log2_segment if a is not None]
				if log2_segment_wo_nones: 
					segment_mean = np.mean(np.array(log2_segment_wo_nones))
					record = [x_axis[start_segment], x_axis[end_segment], segment_mean]
					segment_regions[ch].append(record)
					coords_end.append(end_segment)
			else:
				if not coords_end:
					log2_segment = chr_log2[ch][start_segment:end_segment]
					log2_segment_wo_nones = [a for a in log2_segment if a is not None]
					if log2_segment_wo_nones: 
						segment_mean = np.mean(np.array(log2_segment_wo_nones))
						record = [x_axis[start_segment], x_axis[end_segment], segment_mean]
						segment_regions[ch].append(record)
						coords_end.append(end_segment)

				if counts_start > 1:
					log2_segment = chr_log2[ch][coords_end[-1]:end_segment]
					log2_segment_wo_nones = [a for a in log2_segment if a is not None]
					segment_mean = np.mean(np.array(log2_segment_wo_nones))
					record = [x_axis[coords_end[-1]], x_axis[end_segment], segment_mean]
					segment_regions[ch].append(record)
					coords_end.append(end_segment)

				elif counts_end > 1 and len(chr_segments[ch]) == 1:
					log2_segment = chr_log2[ch][start_segment:indices_e[-1]]
					log2_segment_wo_nones = [a for a in log2_segment if a is not None]
					segment_mean = np.mean(np.array(log2_segment_wo_nones))
					record = [x_axis[start_segment], x_axis[indices_e[-1]], segment_mean]
					segment_regions[ch] = [record]
	
					# coords_end.append(end_segment)

				else:
					for num in indices_e: ## use the next index larger than v[1]
						if num >= v[1]:
							win = indices_e.index(num)
							break
					log2_segment = chr_log2[ch][start_segment:indices_e[win]]
					log2_segment_wo_nones = [a for a in log2_segment if a is not None]
					segment_mean = np.mean(np.array(log2_segment_wo_nones))

					record = [x_axis[start_segment], x_axis[indices_e[win]], segment_mean]
					segment_regions[ch].append(record)
					coords_end.append(indices_e[win])

	return segment_regions


def plot_bins(chr_log2, chromosome_lengths, threshold, segment_regions, output_dir, run, pair):

	# Calculate number of rows and columns for subplot grid
	num_chromosomes = len(chromosome_lengths.keys())
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
		x_axis = [i for i in range(threshold, chromosome_lengths[chromosome] + threshold, threshold)]
		#print(chromosome, chr_log2[chromosome])
		ax.scatter(x_axis, chr_log2[chromosome], s=10)
		ax.set_ylim((-2,2))
		ax.set_title(f"Chromosome {str(chromosome)}")
		ax.set_xlabel('Genomic Position')
		ax.set_ylabel('Log2 Copy Number Ratio')
		ax.grid(True)
		# ax.axhline(y = 0,  color='orange', linestyle='-', linewidth=2) ## mean line at 0
		for l in segment_regions[chromosome]:
			ax.hlines(l[2], xmin=l[0], xmax=l[1], color='orange', linestyle='-', linewidth=2)


	# Hide empty subplots
	for i in range(num_chromosomes, num_rows * num_cols):
		row_index = i // num_cols
		col_index = i % num_cols
		axs[row_index, col_index].axis('off')

	# Adjust layout
	plt.tight_layout()

	# Save or show the plot
	plt.savefig( output_dir + 'CNV_plots_' + run + '_' + pair + '_bins_norm_800K_segements.png')
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

	tsv = sys.argv[1] ##featureCounts file e.g. fc_ALU010.tsv
	pair = sys.argv[2].strip() ## e.g. 1-2 (sample 1-2)
	# output_dir = '/home/labs/maplab/vpascal/Documents/ALU_runs/'
	output_dir = '/Users/Victoria/Downloads/'
	run = tsv.replace('fc_', '').replace('.tsv', '')

	# loc_bams = '/media/vpascal/Expansion/fmoron/projects/Alu-GEUS/data/'
	 
	# files = os.list(loc_bams + run)

	# ##based on the run and pair desired to analyze get the corresponding file in case we want to normalize data based on total number of reads
	# bam_tumor = [f for f in files if f.endswith('.bam') and '_0' + pair[0] +'_001_' in f][0]
	# bam_normal = [f for f in files if f.endswith('.bam') and '_0' + pair[1] +'_001_' in f][0]
	# bamT = loc_bams + run + bam_tumor
	# bamN = loc_bams + run + bam_normal

	# bamT = '/Users/Victoria/Documents/GEUS-Alu/runs/ALU012/bams/run80_03_001_TAGCTT.markdupAnilng.notags.bam'
	# bamN = '/Users/Victoria/Documents/GEUS-Alu/runs/ALU012/bams/run80_04_001_GGCTAC.markdupAnilng.notags.bam'
	# total_N = get_total_read_counts(bamN)
	# total_T = get_total_read_counts(bamT)
	# ratio = float(total_T/total_N)

	## sort the featureCounts file to list the chrs in order (1:22 + X,Y)
	# tsv = '/Users/Victoria/Downloads/duplicate/count_Alu.unique.tsv'
	
	# loc_fc = '/home/labs/maplab/vpascal/Documents/ALU_analysis/cnv_analysis/featureCounts/'
	# num_run = run.replace('ALU0', '')
	# tsv_file = [f for f in os.list(loc_fc) if f == 'fc_RUN' + str(num_run) + '.tsv'][0]
	# tsv = loc_fc + tsv_file
	
	sorted_tsv = tsv.replace('.tsv', '_sorted.tsv')
	if not sorted_tsv in os.listdir(output_dir):
		sort_tsv_file(tsv, sorted_tsv)
	

	##set bin thresholds
	threshold = 1000000

	##Bin the chromosome
	chr_bin_counts_tumor, chr_bin_counts_normal, total_T, total_N = bin_fc_file(chromosome_lengths, sorted_tsv, threshold, pair, output_dir, run )
	# ratio = float(total_T/total_N)
	ratioT = float(800000/total_T) ## using the treshold of 800K
	ratioN = float(800000/total_N) ## using the treshold of 800K
	print(total_T, total_N)

	##save the files into a dict
	json.dump(chr_bin_counts_tumor, open(output_dir + run + '_read_counts_T_dict_sample_' + pair + '.txt','w'))
	json.dump(chr_bin_counts_normal, open(output_dir + run + '_read_counts_N_dict_sample_' + pair + '.txt','w'))

	# chr_bin_counts_tumor = json.load(open('/Users/Victoria/Downloads/duplicate/run12_chr_read_counts_tumor_dict_sample3-4.txt'))
	# chr_bin_counts_normal = json.load(open('/Users/Victoria/Downloads/duplicate/run12_chr_read_counts_normal_dict_sample3-4.txt'))

	#Compute log2(T) - Log2(N), considering normalized N counts using the ratio
	chr_log_values = compute_log2_tumor_normal(chr_bin_counts_tumor, chr_bin_counts_normal, ratioT, ratioN, output_dir, run, pair)

	##Smooth based on log2 values
	# chr_log_values_smoothed, value_conversion = smooth_log2values_bin(chr_log_values)

	##Segment regions based on log2 values
	tuple_index = get_segments(chr_log_values)
	segment_coords = get_segment_regions(chromosome_lengths, tuple_index, chr_log_values) ##alternate between smoothed and non-smoothed

	##plot the log2 values per bin and segments
	plot_bins(chr_log_values, chromosome_lengths, threshold, segment_coords, output_dir, run, pair)

	##save log2 counts and segments in a dict
	json.dump(chr_log_values, open(output_dir + run + '_log2_sample_' + pair + '.txt','w'))
	json.dump(segment_coords, open(output_dir + run + '_segments_sample_' + pair + '.txt','w'))
