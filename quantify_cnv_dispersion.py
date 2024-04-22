import json
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def quantify_global_cn_changes(inf):
	cn_ch = {}
	for ch in inf.keys():
		filtered = [num for num in inf[ch] if not num == None]
		glo = "{:.2f}".format(float(sum(filtered) / len(filtered))) ## mean log2
		cn_ch[ch] = glo

	return cn_ch

def quantify_local_cn_changes(inf, seg, thr, chr_len):
	ch_local_mean = {}
	for ch in seg.keys():
		for s in seg[ch]:
			size = float((s[1] - s[0]) / thr)
			prop =  "{:.2f}".format(float(size / float(chr_len[ch]/thr)) * 100)

			if size + 1 == len(inf[ch]):
				filtered = [num for num in inf[ch] if not num == None]
				mean = "{:.2f}".format(float(sum(filtered) / len(filtered)))
				ch_local_mean[ch] = [str(mean), '100']
			else:
				
				start = int((s[0] / thr) -1)
				end = int((s[1] / thr) -1)
				# print(ch, prop, size)
				filtered = [num for num in inf[ch][start:end] if not num == None]
				lo = "{:.2f}".format(float(sum(filtered) / len(filtered)))
				rec = [str(lo), str(prop)]
				if ch in ch_local_mean.keys():
					# ch_local_mean[ch].append(rec)
					ch_local_mean[ch].append(str(lo))
					ch_local_mean[ch].append(str(prop))
				else:
					ch_local_mean[ch] = rec


	return ch_local_mean

def plot_segments(chr_segments, chromosome_lengths, thr):

	num_chromosomes = len(chr_segments.keys())
	num_cols = 5
	num_rows = -(-num_chromosomes // num_cols)  # Ceiling division to ensure enough rows

	# # Create a grid of subplots
	fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 20))
	for i, ch in enumerate(chr_segments.keys()):
		row_index = i // num_cols
		col_index = i % num_cols
		ax = axs[row_index, col_index] if num_rows > 1 else axs[col_index]
	   
		for segment in chr_segments[ch]:
			start_x, end_x, value = segment
			rounded_value = (chromosome_lengths[ch] // thr) * thr
			print(ch, rounded_value, end_x-start_x)
			if (end_x- start_x) == rounded_value:
				ax.plot([start_x, end_x], [value, value], color='#3C5B6F', linewidth=0.5)
			else:
				ax.plot([start_x, end_x], [value, value], color='#FA7070', linewidth=0.5)

		# Customize the plot
		ax.set_ylim((-1.5,1.5))
		ax.set_xlabel('Length')
		ax.set_ylabel('Value')
		ax.set_title(f"Segments chr {str(ch)}")
		# ax.grid(True)

	#Hide empty subplots
	for i in range(num_chromosomes, num_rows * num_cols):
		row_index = i // num_cols
		col_index = i % num_cols
		axs[row_index, col_index].axis('off')

	# Adjust layout
	plt.tight_layout()

	# Show the plot
	# plt.show()
	labels = ['Whole', 'Fragmented']
	colors = []
	colors.append(mpatches.Patch(color='#3C5B6F'))
	colors.append(mpatches.Patch(color='#FA7070'))
	fig.legend(colors, labels, loc='lower right',bbox_to_anchor=(0.95, 0.1), shadow=True, ncol=2)
	plt.savefig('CNV_segments_plot.png') ##variable for the number of chs



if __name__ == "__main__":
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
	    '13': 96108862, ##custom
	    '14': 88287443, ##custom
	    '15': 82491545,##custom
	    '16': 90282411,
	    '17': 81192741,
	    '18': 78009649,
	    '19': 59114562,
	    '20': 62963629,
	    '21': 39117466, ##custom
	    '22': 35244339, ##custom
	    'X': 155246021,
	    'Y': 59033256}

	folder = [file for file in os.listdir('/Users/Victoria/Downloads/segments+log2') if 'log2' in file and file.endswith('txt')]
	out_dir = '/Users/Victoria/Downloads/CNV_quantification/'
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)

	chr_segments = {}
	for file in folder:
		infile = '/Users/Victoria/Downloads/segments+log2/' + file
		segments = '/Users/Victoria/Downloads/segments+log2/' + file.replace('log2', 'segments')
		inf = json.load(open(infile))
		seg = json.load(open(segments))

		threshold = 1000000
		# glo = quantify_global_cn_changes(inf)
		# loc = quantify_local_cn_changes(inf, seg, threshold, chromosome_lengths)
		# with open(out_dir + 'Quantification_cnv_changes_' + file.replace('_log2', ''), 'w') as out:
		# 	out.write("Chromosome\tGlobal CNV\tLocal CNV mean \t %Sequence\n")
		# 	for ch, val in loc.items():
		# 		out.write(str(ch) + '\t' + str(glo[ch])  + '\t' + '\t'.join(val) + '\n')

	##plot segments all samples in one plot
		for ch in seg.keys():
			if ch in chr_segments.keys():
				for a in seg[ch]:
					chr_segments[ch].append(a)
			else:
				chr_segments[ch] = seg[ch]

	plot_segments(chr_segments, chromosome_lengths, threshold)


	







