import pandas as pd 
import numpy.ma as ma
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
from string import ascii_letters
import seaborn as sns
import os.path
import itertools    
import re
import os

from string import ascii_letters

def get_NORAD_score(target,position,temp):
	file_open = directory + "/" +  target +  "_" + str(temp) +  ".csv"

	NORAD_score  = pd.read_csv(file_open,sep=';',header=None,index_col = 0) 

	# Remove last position
	return NORAD_score.iloc[:, position:position + 179]

def get_UNIT_score(MOLECULE,temp):

	file_open = directory + "/" +  MOLECULE +  "." + str(temp) +  ".csv"
	
	if (os.path.isfile(file_open)):
		UNIT_score  = pd.read_csv(file_open,sep=';',header=None,index_col = 0) 
		# Remove last position
		return UNIT_score.iloc[:, :-1]
	else:
		print("No encuentro",file_open)
		return pd.DataFrame({'P' : []})



def convertphastConsToTab(name,reverse, chrName, chrSt, chrEnd,BigWig_file):
	# Check chr start and end with http://genome.ucsc.edu/cgi-bin/hgBlat
	
	import pyBigWig
	bw = pyBigWig.open(BigWig_file)

	# Wiggle, bigWig, and bigBed files use 0-based half-open coordinates, which are also used by this extension.
	chrSt = chrSt - 1
	df = pd.DataFrame(bw.values(chrName, chrSt, chrEnd),columns=['Column_Name'])
	
	if (reverse):
		df_save = df.iloc[::-1,]
	else:
		df_save = df

	save = 'conservation_NORAD/' + name + '_phastCons.csv'
	
	df_save.T.to_csv(save,index=False)

	return df_save


def plot_correlation_matrix(dataframe,MOLECULE):
	# ~ sns.set_theme(style="white")
	sns.set(style="white")
	
	# Compute the correlation matrix
	corr = dataframe.T.corr(method='pearson')

	# Generate a mask for the upper triangle
	mask = np.triu(np.ones_like(corr, dtype=bool))

	# Set up the matplotlib figure
	f, ax = plt.subplots(figsize=(18, 13))

	# Generate a custom diverging colormap
	# ~ cmap = sns.diverging_palette(230, 20, as_cmap=True)
	cmap = sns.diverging_palette(220, 20, as_cmap=True)
	
	
	# Draw the heatmap with the mask and correct aspect ratio
	sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1,vmin=-1, center=0,
				square=True, linewidths=.5, cbar_kws={"shrink": .5},annot=True)
	
	output_file = 'plot_NORAD/' + MOLECULE + '.pdf'
	figure = ax.get_figure()    
	figure.savefig(output_file, dpi=400)
	figure.clf()
	figure.clf()
	
	plt.close(f)

def plot_box_plot(df_pos,df_neg,MOLECULE,BigWigFile_list,df_phastCons_all):
	
	# Figure
	sns.set_style("whitegrid") 
	
	df_pos['type'] = 'NRU_patt'
	
	df_neg = df_neg.T
	df_neg['type'] = 'NRU_rest'


	df_phastCons_all = df_phastCons_all.T
	df_phastCons_all.columns=list(BigWigFile_list)
	df_phastCons_all['type'] = 'all'
	# Drop first row
	df_phastCons_all = df_phastCons_all.iloc[1:]
	


	frames = [df_pos, df_neg,df_phastCons_all]
	df = pd.concat(frames,sort=False)


	# Default regular plot
	df.boxplot(by='type',column=list(BigWigFile_list), grid = False,figsize=(12,9))
	pdf_out = 'plot_NORAD/' + MOLECULE + '_boxplot.pdf'
	plt.savefig(pdf_out, dpi=800)
	plt.clf()
	
	
	# TODO: Check this
	# Plot with statistics
	frames = [df_pos, df_phastCons_all]
	df = pd.concat(frames,sort=False)

	# Plot with Statistics
	# ~ from scipy import stats

	# ~ df_long = pd.melt(df, 'type', var_name='Feature', value_name='Value')
	# ~ all_t = list()
	# ~ all_p = list()

	# ~ for case in range(len(BigWigFile_list)):
		# ~ sub_df = df_long[df_long.Feature == BigWigFile_list[case]]
		# ~ g1 = sub_df[sub_df['type'] == 'NRU_patt']['Value'].values
		# ~ g2 = sub_df[sub_df['type'] == 'all']['Value'].values
		
		# ~ # Calculate the T-test for the means of two independent samples of scores.
		# ~ t, p = stats.ttest_ind(g1, g2)
		# ~ all_t.append(t)
		# ~ all_p.append(p)
	# ~ # If the t-value is positive (>0) then the mean of g1 (group 1: NRU_patt) was significantly greater than the mean of g2 (group 2: all).
	# ~ # If the t-value is negative (<0) then the mean of g1 (group 1: NRU_patt) was significantly smaller than the mean of g2 (group 2: all).
	# ~ print(all_t)
	# ~ print(all_p)
	# ~ # We can see that there is a statistically significant difference in all 4 features between NRU_patt and all.
	# ~ x = np.count_nonzero(np.array(BigWigFile_list)[np.array(all_p) < 0.05])
	# ~ print("Therer are ", x , " test that are statistically significant different" )
	
	# ~ # And now the plot
	# ~ # renaming so that class 0 will appear as ...
	# ~ df_long.loc[df_long.type==0, 'type'] = BigWigFile_list[0]
	# ~ df_long.loc[df_long.type==1, 'type'] = BigWigFile_list[1]
	# ~ df_long.loc[df_long.type==2, 'type'] = BigWigFile_list[2]
	# ~ df_long.loc[df_long.type==3, 'type'] = BigWigFile_list[3]
	# ~ df_long.loc[df_long.type==4, 'type'] = BigWigFile_list[4]

	# ~ # Boxplots
	# ~ fig, axes = plt.subplots(2,3, figsize=(14,10), dpi=100)
	# ~ axes = axes.flatten()

	# ~ for idx, feature in enumerate(BigWigFile_list):
		# ~ ax = sns.boxplot(x="Feature", hue="type", y="Value", data = df_long[df_long.Feature == feature], linewidth=2, showmeans=True, meanprops={"marker":"*","markerfacecolor":"white", "markeredgecolor":"black"}, ax=axes[idx])
		# ~ #* tick params
		# ~ axes[idx].set_xticklabels([str(feature)], rotation=0)
		# ~ axes[idx].set(xlabel=None)
		# ~ axes[idx].set(ylabel=None)
		# ~ axes[idx].grid(alpha=0.5)
		# ~ axes[idx].legend(loc="lower right", prop={'size': 11})
	 
		# ~ #*set edge color = black
		# ~ for b in range(len(ax.artists)):
			# ~ ax.artists[b].set_edgecolor('black')
			# ~ ax.artists[b].set_alpha(0.8)
	   
		# ~ #* statistical tests
		# ~ x1, x2 = -0.20, 0.20
		# ~ y, h, col = df_long[df_long.Feature == feature]["Value"].max()+1, 2, 'k'
		# ~ axes[idx].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
		# ~ axes[idx].text((x1+x2)*.5, y+h, "statistically significant", ha='center', va='bottom', color=col)
	# ~ fig.suptitle("Significant feature differences between setosa and versicolor classes/groups", size=14, y=0.93)
	# ~ plt.show()
		

def plot_box_plot_final(df_pos,df_phastCons_all,name,BigWig_no_path):
	
	# Figure
	sns.set_style("whitegrid") 
	
	df_pos['type'] = 'NRU_patt'

	df_phastCons_all = df_phastCons_all.T
	df_phastCons_all.columns=list(BigWig_no_path)
	df_phastCons_all['type'] = 'all'

	# Drop first row
	df_phastCons_all = df_phastCons_all.iloc[1:]


	frames = [df_pos, df_phastCons_all]
	df = pd.concat(frames,sort=False)
	
	df.boxplot(by='type',column=list(BigWig_no_path), grid = False,figsize=(12,9))
	
	pdf_out = 'plot_NORAD/' + name + '_boxplot.pdf'
	plt.savefig(pdf_out, dpi=800)
	
	plt.clf()


def stability_pattern(input_fasta,temps,pattern,st_pos,end_pos,chrName, name, reverse, chrSt, chrEnd,BigWigFile_list):
		
	####################
	# Get phastcons data
	####################
	
	# list of dataframes
	df_phastCons_original = {}
	
	for BigWigFile in BigWigFile_list:
		df_phastCons_original[BigWigFile] = convertphastConsToTab(name,reverse,chrName, chrSt, chrEnd,BigWigFile).T

		# And insert one column
		df_phastCons_original[BigWigFile].insert(0, 'x', os.path.basename(BigWigFile))
		
		# define N in function of number of columns
		N = df_phastCons_original[BigWigFile].shape[1]
		df_phastCons_original[BigWigFile].columns = range(0, N)
		
	####################
	# Get NORAD units data
	####################
	
	data_all_units_positive = []
	
	with open(input_fasta) as fp:
		line = fp.readline()
		while line:
			MOLECULE = line.strip()[1:]
			sequence = fp.readline().strip()
			position = NORADS[0].find(sequence.rstrip())
			target = NORADS_NAME[0]
			
			i=0
			while(position == -1):
				i=i+1
				position = NORADS[i].find(sequence.rstrip())
				target = NORADS_NAME[i]
			
			comb = itertools.combinations(temps, 2)
							
			temp1 = temps[0]
			temp2 = temps[1]
			temp3 = temps[2]
				
			# First compare the units at temperatures
			UNIT_score_first = get_UNIT_score(MOLECULE,temp1)
			UNIT_score_first.index = UNIT_score_first.index + "_" + str(temp1)
			
			UNIT_score_second = get_UNIT_score(MOLECULE,temp2)
			UNIT_score_second.index = UNIT_score_second.index + "_" + str(temp2)
			
			UNIT_score_third = get_UNIT_score(MOLECULE,temp3)
			
			if(UNIT_score_third.index.empty):
				UNIT_score_third.index = {}
			else:
				UNIT_score_third.index = UNIT_score_third.index + "_" + str(temp3)
				
											
			frames = [UNIT_score_first, UNIT_score_second, UNIT_score_third]
			df = pd.concat(frames)

		
			# Crop dataframe
			df = df.iloc[:, st_pos:end_pos]
			df=df.T

			UNIT_score_df = df
			

			# grab the first row for the header
			new_header = UNIT_score_df.iloc[0] 

			# Get dataframe and calculate mean and stdev
			df_temp = UNIT_score_df
			UNIT_score_df['NRU_avg'] = df_temp.mean(axis=1)
			UNIT_score_df['NRU_std'] = df_temp.std(axis=1)



			####################
			# Whole NORAD (5 6 and 7)
			####################
			
			# Then compare NORAD
			NORAD_score_df_first = get_NORAD_score(target,position,temp1)
			# Rename columns
			NORAD_score_df_first.columns = np.arange(1,len(NORAD_score_df_first.columns)+1)
			# Rename index
			NORAD_score_df_first.index = NORAD_score_df_first.index + "_" + str(temp1)
			
			NORAD_score_df_second = get_NORAD_score(target,position,temp2)
			NORAD_score_df_second.columns = np.arange(1,len(NORAD_score_df_second.columns)+1)
			NORAD_score_df_second.index = NORAD_score_df_second.index + "_" + str(temp2)
			
			NORAD_score_df_third = get_NORAD_score(target,position,temp3)
			NORAD_score_df_third.columns = np.arange(1,len(NORAD_score_df_third.columns)+1)
			NORAD_score_df_third.index = NORAD_score_df_third.index + "_" + str(temp3)
			
			frames = [NORAD_score_df_first, NORAD_score_df_second,NORAD_score_df_third ]
			df = pd.concat(frames)

			df = df.iloc[:, st_pos:end_pos]
			df=df.T

			NORAD_score_df = df
			
			new_header = NORAD_score_df.iloc[0] #grab the first row for the header

			# Get dataframe without phastCons and calculate mean and stdev
			df_temp = NORAD_score_df
			NORAD_score_df['NORAD_avg'] = df_temp.mean(axis=1)
			NORAD_score_df['NORAD_std'] = df_temp.std(axis=1)
			
			######################
			# Processing phastCons
			######################
			
			# TODO I've checked and I have to move 1 position (to be the same as the NRU paper) -> TODO, why? 
			offset = 1
			position_in_whole_NORAD = WHOLE_NORAD.find(sequence.rstrip()) + offset

			# Fragment of phastCons NORAD unit
			end_position_in_whole_NORAD = position_in_whole_NORAD + len(sequence.rstrip())
			
			df_phastCons = {}
			
			for BigWigFile in BigWigFile_list:
				# Get a copy of the df_phastCons dataframe
				df_phastCons[BigWigFile] = df_phastCons_original[BigWigFile]
					
				# Slice df_phastCons to UNIT position
				idx=df_phastCons[BigWigFile][0]

				df_phastCons[BigWigFile] = df_phastCons[BigWigFile].iloc[:, position_in_whole_NORAD:end_position_in_whole_NORAD].T
				# Take only the central part
				df_phastCons[BigWigFile] = df_phastCons[BigWigFile].T.iloc[:, st_pos:end_pos].T
				# Rename NORAD score indexes
				df_phastCons[BigWigFile].index = list(NORAD_score_df.index)
				df_phastCons[BigWigFile].columns=idx
			
			# Traspose all dataframes and save them to df_conservation
			for BigWigFile in BigWigFile_list:
				df_phastCons[BigWigFile] = df_phastCons[BigWigFile].T
							
			df_conservation = pd.concat(df_phastCons)
			
			# Reset multiIndex
			BigWig_no_path = list()
			# Remove path and remove extension
			for BigWig in BigWigFile_list:
				BigWig_no_path.append(os.path.basename(BigWig).replace('.bw','').replace('hg38.',''))
				
			df_conservation.index=BigWig_no_path
			
			######################
			# Merge all
			######################
				
			frames = [NORAD_score_df.T, UNIT_score_df.T,df_conservation]

			final_df = pd.concat(frames)
			print (final_df.T.corr())
			
			# TODO CHECK OUTPUT
			csv_file = "out_NORAD" + "/"  +  MOLECULE + '_' +  str(st_pos) + "_" + str(end_pos) +  "_UNIT.csv" 
			final_df.to_csv(csv_file)
						
			csv_file = "out_NORAD" + "/"  +  MOLECULE + '_' +  str(st_pos) + "_" + str(end_pos) +  "_UNIT_corr.csv" 
			final_df.T.corr().to_csv(csv_file)
						
			######################
			# Processing phastCons
			######################
			
			# Add random column
			final_df.loc['random'] = np.random.randint(1, 100, final_df.shape[1])


			plot_correlation_matrix(final_df,MOLECULE)
			
			# Correlation by window
			
			WINDOWS = (3,10,20,30)
			for WIND in WINDOWS:
				new_df = final_df.T.rolling(WIND).mean().T
				
				# Drop the rows where all elements are missing.
				# ~ new_df = new_df.T.dropna(how='all').T
				# Drops all rows containing at least one field with missing data
				# ~ new_df = new_df.T.dropna().T

				new_MOLECULE = MOLECULE + "_" + str(WIND)
				plot_correlation_matrix(new_df,new_MOLECULE)
			
			

			####################
			# Pattern sequence
			####################
			
			positive_col = list()
			negative_col = list()
			all_col = list(final_df.columns)
			if pattern:
				positions_pattern = list()
				if (not re.findall(pattern, sequence)):
					print("No pattern was found")
				
				else:
					# Process all patterns
					for m in re.finditer(pattern, sequence):
						offset = 1
						start = m.start() + offset
						end = m.end() + offset
						print(m,start,end)
						for i in  range(start,end):
							positive_col.append(i)
					
					negative_col = list(set(all_col) - set(positive_col))
					
					df_positive = final_df[final_df.columns.intersection(positive_col)]
					df_negative = final_df[final_df.columns.intersection(negative_col)]
					
					
					MOLECULE_pos = MOLECULE + "_" + pattern
					MOLECULE_neg = MOLECULE + "_NO_" + pattern
					plot_correlation_matrix(df_positive,MOLECULE_pos)
					plot_correlation_matrix(df_negative,MOLECULE_neg)
					
					df_phastCons_all = pd.concat(df_phastCons_original)
					
					# Save to csv # Only pattern positions
					csv_file = "out_NORAD" + "/"  +  MOLECULE_pos + ".csv" 
					df_positive.to_csv(csv_file)
					
					df_positive = df_positive.T
					


					df_positive = df_positive[BigWig_no_path]

					data_all_units_positive.append(df_positive)
					

					plot_box_plot(df_positive,df_negative,MOLECULE_pos,BigWig_no_path,df_phastCons_all)
			else:
				print("No pattern was given")
			
			# Todo: remove this
			# ~ return
			
			line = fp.readline()
		df_all_units_positive = pd.concat(data_all_units_positive)
		name = "all_" + pattern 

		plot_box_plot_final(df_all_units_positive,df_phastCons_all,name,BigWig_no_path)


temperatures = (23,37,55)
directory = 'data_NORAD'


WHOLE_NORAD='AGTTCCGGTCCGGCAGAGATCGCGGAGAGACGCAGAACGCAGCCCGCTCCTCCAGGGCCCTCCAGGCCCTCCGGCCCCGGGCCGGCGGGTGAACTGGGGGGCCCCGGGACAGGCCGAGCCCTCTGCCCTGCAGATAACGGAGGCCTCTGCTGTGGCTGCCCACTGGCTGTGCCCGCCCACTGGCTGTGCCCAGACCTTGAAGCCGCAGCGAACCTCTCTTTCCCACCCCACCTCGGTGACTAATGGCGGCCGTGGCGTCTCCCAGCCCGGACCCCGCCGGCACCCGGGTCTCCCGACCCAAGCCTCGACGAAACCCCCGCAGAGCCGCCGGGACGCAGCGCCTTTGGGCGGCGCTGGGCGTGGTGGGCCGGGAAGTATGGCGGCAGCTCGAACGCCGCGCGGCGGAGGCCATTAAGGCGTGGACGGCCCGGGAAGGCGGCCTAGGGACGCAAGCAGGCTCGGCCGCCTCTTTAGGCCACGGAGCCGCGCAGATCCGGTTCCCGGGTGACCACTCTGTCGCCATTGGGCGAGACCTACCTAGTCCTGACGACAACGGACAAAGGCCTTAAGGGGCCTGGAAGGTGAGCGAAGTCCCGAACGACGACGGGTGGAACGGTTAGCGGCCATCGGGCGGTTGGTCTTCATTCTACCAGACTTTGCTGTCGGAAGAGAGAAATGGTAGAATGACAGGCCACGTTTGGCCCGTTGGAAATGCCCACCACCCTCTGGGAAGATTTACTGGCCGTTTATGGAAGGCCTGTGTATATAATATGAAAAAGCTGCTCTCAACTCCACCCCAACCTTTTAATAGAAAACATTTGTCACATCTAGCCCTTCTAGATGGAAAGAGGTTGCCGACGTATGATAAAATAGAGTTAGAAAGTTACACATCTTGTAAATTCTCATTTGTTTAAAAGAAATCATAGAAAATACATGTCTTCTGGAGATGACTTTTGGAAATGGAGTTGTTAAGACGGCCTCTGGAAGCGATACGTCCACGTTTGTTAAGTGGGTTAGATGACATGGAGCTGGAAGACCTGAGAAGGAAGAGAAGAAGGTTCTATGCTAGACTGGTCATATTTAGAAGACATTTTCATATTCTATCCATTGTTTTGTGTGCATTTTATTCCTCACTACTGTGTATATAGTTGACAATGCTAAGCTTTTTTGAAATGTCTCTTCTTTTTAGATGTTCTGAAGTGCCTGATATGTTAAAATTAGAGGTAGCAAAATCACATTTTGTAAATACCTTTTTGTTACAATTCATAGGAAATATTTTTGGGGGGGAATGGCCAAATCACCTGTTGAGTAATACTCATTGTGTTTGTGCAGTGGTTCAGGGGAGGAGAGAGGAGGGGGAGGTGCAGAGAGCTCTATGCCATCCTGTTTACAGCGAGGCAAGATGAATCATTATGTCTGTGCATTTTGTTTTACTTATCTGTGTATATAGTGTACATAAAGGACAGACGAGTCCTAATTGACAACATCTAGTCTTTCTGGATGTTAAAGAGGTTGCCAGTGTATGACAAAAGTAGAGTTAGTAAACTAATATATTTTGTACATTTTGTTTTACAAGTCCTAGGAAAGATTGTCTTCTGAAAATTTGATGTCTTCTGGGTTGATGGAGATGGGAAGGGTTCTAGGCCAGAATGTTCACATTTGGAAGACTCTTTCAAATTATAACTGTTGTTACATGTTTGCAGTTTATTCAAGACTGCTGTATACATAGTAGACAAATTAACTCCTTACTTGAAACATCTAGTCTATCTAGATGTTTAGAAGTGCCCGATGTATGTTAAATGTATAGGTAGTAAAATACCACTTTGTAAATATCTTTTTGCTAAAATTCATAGGAAATGCTTTTGGAAATTGAATTGTGAAGCCACCTTTGTGAACAGTATAGTAATGTCTATACTTGTTCAATAGTTTAGAGGAGGTAGGAGGGAAGAAATTGCAAAAGGTAATATTACTAGTGTGTTCATACTTGGACATTTTCAGACACCATTTTTCTATATGTTTTGTGCATTTTGTTTTGCTCTGTATATAGTATATATAATGGACAAATAGTCCTAATTTTTCAACATCTAGTCTCTAGATGTTAAAGAGGTTGCCAGTGTATGACAAAGGAGTAAAATTAGCATATTTTGTACACTTTGTGTTGAAATTCGTAGGAAAACTTGTCTTCTGTAAAGACTTTTGCATAGGAATTTGTTTGACCATCTCTAAGCATTACACGTGCCTGTACTTGTCCACTGGATTGAAGGCAGAGAAGGAAGGGAGGAGGGAATGATTCAAGGCCAAAATGGCCACATTTAGAAGATACCTCAGATGATAACCATTGTTATGTGTGTGCAATTTTATTTAACAGTGCTGTGTATGTGGTGGACAAGTTATATGAAATATCTAGTCTTTCTAGATATTTGGAAGTGCTTGATGTATTTAAAAGTGGTAGTAGAATAACACTTTGTAAATAGCTTTTAAAAACTGATGGGAAATGCTGTTTGGAAGTGGAATTGTTGAACCACCTGGGAGGTGGGAGGGAAGAAATTGCAAATGGTGTTTTGCCATTGTTTATTAGAAAATTTCAGCTTAATCCATTGTGTATATGTTACATGCATTTCATTTAACTTTGCTATACTGTATATATTGTATATATAACGGACAAATTAGTCCCGATTTTATAATATCTAGTCTCTAGATATTAAAGAGGTTGCCAATGTATGACAGAAGTAGAGTTAGTAAACTAACACATTTTGTACACTTTGTTAAAATTTGTAGAAAGGCTGTCTTCTGAAAAGGACTTTTGGAAGTGAGATAACATCAGCTCTAAGTGACACGTGCCTATATCCATCAGGTTGGTGGTGGAGAGGAGTTGGAAGGAATGAAGGGTTCTAGACCAGAATGTTCGTATTTAGAAGACACTATCAGATATAACCATTGTTACATGTGTGTAGTTTATTCAACCCTACTGTGTATATAGCGGACAAACTTAAGTCCTTATTTGAAACATCTAGTCTTTCTAGATGTTTAGAAGTGCACAAAGTATGTTAAAAGTAGAGGTAGTAAATAACACATTTTGTAGCTATCCTTTTGATATGAAATATTGTCTTGGAAATTGATCAATTCTCTGAGCAGTACCCATTTTGATATTTGTGCTGGTTCAGGGGGAAGGAGGAGCACAAAGTGCAAAGGGCTTTCTACCAGTGTCCAGTGTGTTTATGAGGAGGCACATTGACCATTGTCCCTTATGTCTGCATTTTCATTTACTGTGCTGTGTATATAGTGTATATAAGCGGACATAGGAGTCCTAATTTACGTCTAGTCGATGTTAAAAAGGTTGCCAGTATATGACAAAAGTAGAATTAGTAAACTACTACATTGAGTACACTTTGTGTTAAAATTCATAGGGAAGACTTCTTAAAAACAAGTGAAATTGTTAAAACCCCCCCTAAGCATTACAGATGGCTTATAGCTGTCCACGGGGTTGGTAGAGGTGGGAAAGGGAAGGGTTCTAGGCCAGAATGTTCCTATTTAGAAGACACTCAAATTACAGTCTGTGTTATGTATGTATACCATTTATTCAATGCTACTGTGTATATAATGGAAAACTTAAGTCCAGTTTGAAACATCTAGTCTTTCTAGGTGTTTAAAAGTGTACAACGGCCTGTCGCAGTGGCGCATGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCAGGCGGATCACGAGGTCAAGAGATCAGGACCATCTTGGCCAACATGGTGAAACCCCATCTTTACTAAAAATACAAAAATTAGCTGGTCGTGGTGGTGCCCACCTGTAGCCCCAGTTACTCGAGAGGCTGAGGCAGGAGAATCGCTTGAACTTGGGAGGCGGAAGTTGCAGTGAGCCAAGATCGCACCACTGCACTCCAGCCTGGCGACAGAGCGAGGCTCCGTTTCAAAAAAAAAAGTGCACAATGTAGGTTAACAGTAGAGGGCTTAAGTAACACCCCTCTAAGCATTTGTTTTCAGTACTTCCTAGGAGTGGTTGCATTTGGGAATGGAATTGTTAAAACTTGATGCTTAGGAGCGAATGCAGACTATTCATTGGGTGTTTGGGGTGGGGGAAGGGGGGGTGGGCAGAGGAGGTATGCAGGGAGAGGGGTTCTGTGCTCCTGAGATTAGTTCAGATGGTCTAACCATTGTTCTATATGTGCATTTTAGTTAATATTGTGTATTAAAGGATAAGTCTTAATGCTCAAAGTATGTTAAAAATAGATGTAGTAAATCAGTCCCTTTGTGAATGTCCTTTTGTTAGTTTTTAGGAAGGCCTGTCCTCTGGGAGTGACCTTTATTAGTCCACCCCTTGGAGCTAGACATCCTGTACTTAGTCACGGGGATGGTGGAAGAGGGAGAAGAGGAAGGGTGAAGGGAAGGGCTCTTTGCTAGTATCTCCATATCTAGACGATGGTTTTAGATGATAACCACAGGTCTACAAGAGCGTTTTTAGTAAAGTGCCTGTGTTCATTGTGGACAAAGTTATTATTTTGCAACATCTAAGCTTTACGAATGGGGTGACAACTTATGATAAAAACTAGAGCTAGTGAATTAGCCTATTTGTAAATACCTTTGTTATAATTGATAGGATACATCTTGGACATGGAATTGTTAAGCCACCTCTGAGCAGTGTATGTCAGGACTTGTTCATTAGGTTGGCAGCAGAGGGGCAGAAGGAATTATACAGGTAGAGATGTATGCAGATGTGTCCATATATGTCCATATTTACATTTTGATAGCCATTGATGTATGCATCTCTTGGCTGTACTATAAGAACACATTAATTCAATGGAAATACACTTTGCTAATATTTTAATGGTATAGATCTGCTAATGAATTCTCTTAAAAACATACTGTATTCTGTTGCTGTGTGTTTCATTTTAAATTGAGCATTAAGGGAATGCAGCATTTAAATCAGAACTCTGCCAATGCTTTTATCTAGAGGCGTGTTGCCATTTTTGTCTTATATGAAATTTCTGTCCCAAGAAAGGCAGGATTACATCTTTTTTTTTTTTTTTAGCAGTTTGAGTTGGTGTAGTGTATTCTTGGTTATCAGAATACTCATATAGCTTTGGGATTTTGAATTGGTAAATATTCATGATGTGTGAAAAATCATGATACATACTGTACAGTCTCAGTCCCATAAAATTGGATGTTGTGCCTACACACAGGATCTAGAAGAATATGTCAAACTATAAACTGCTTGTGATTGTGAATGACTTTGTTCTTTGCTTGTGTTTTTCAATTTCCTATAATGCACATACTAACTTTTAAAAAATAAAGGTTATTTTAAAAGCCTGTA'


NORAD5='AGTTCCGGTCCGGCAGAGATCGCGGAGAGACGCAGAACGCAGCCCGCTCCTCCAGGGCCCTCCAGGCCCTCCGGCCCCGGGCCGGCGGGTGAACTGGGGGGCCCCGGGACAGGCCGAGCCCTCTGCCCTGCAGATAACGGAGGCCTCTGCTGTGGCTGCCCACTGGCTGTGCCCGCCCACTGGCTGTGCCCAGACCTTGAAGCCGCAGCGAACCTCTCTTTCCCACCCCACCTCGGTGACTAATGGCGGCCGTGGCGTCTCCCAGCCCGGACCCCGCCGGCACCCGGGTCTCCCGACCCAAGCCTCGACGAAACCCCCGCAGAGCCGCCGGGACGCAGCGCCTTTGGGCGGCGCTGGGCGTGGTGGGCCGGGAAGTATGGCGGCAGCTCGAACGCCGCGCGGCGGAGGCCATTAAGGCGTGGACGGCCCGGGAAGGCGGCCTAGGGACGCAAGCAGGCTCGGCCGCCTCTTTAGGCCACGGAGCCGCGCAGATCCGGTTCCCGGGTGACCACTCTGTCGCCATTGGGCGAGACCTACCTAGTCCTGACGACAACGGACAAAGGCCTTAAGGGGCCTGGAAGGTGAGCGAAGTCCCGAACGACGACGGGTGGAACGGTTAGCGGCCATCGGGCGGTTGGTCTTCATTCTACCAGACTTTGCTGTCGGAAGAGAGAAATGGTAGAATGACAGGCCACGTTTGGCCCGTTGGAAATGCCCACCACCCTCTGGGAAGATTTACTGGCCGTTTATGGAAGGCCTGTGTATATAATATGAAAAAGCTGCTCTCAACTCCACCCCAACCTTTTAATAGAAAACATTTGTCACATCTAGCCCTTCTAGATGGAAAGAGGTTGCCGACGTATGATAAAATAGAGTTAGAAAGTTACACATCTTGTAAATTCTCATTTGTTTAAAAGAAATCATAGAAAATACATGTCTTCTGGAGATGACTTTTGGAAATGGAGTTGTTAAGACGGCCTCTGGAAGCGATACGTCCACGTTTGTTAAGTGGGTTAGATGACATGGAGCTGGAAGACCTGAGAAGGAAGAGAAGAAGGTTCTATGCTAGACTGGTCATATTTAGAAGACATTTTCATATTCTATCCATTGTTTTGTGTGCATTTTATTCCTCACTACTGTGTATATAGTTGACAATGCTAAGCTTTTTTGAAATGTCTCTTCTTTTTAGATGTTCTGAAGTGCCTGATATGTTAAAATTAGAGGTAGCAAAATCACATTTTGTAAATACCTTTTTGTTACAATTCATAGGAAATATTTTTGGGGGGGAATGGCCAAATCACCTGTTGAGTAATACTCATTGTGTTTGTGCAGTGGTTCAGGGGAGGAGAGAGGAGGGGGAGGTGCAGAGAGCTCTATGCCATCCTGTTTACAGCGAGGCAAGATGAATCATTATGTCTGTGCATTTTGTTTTACTTATCTGTGTATATAGTGTACATAAAGGACAGACGAGTCCTAATTGACAACATCTAGTCTTTCTGGATGTTAAAGAGGTTGCCAGTGTATGACAAAAGTAGAGTTAGTAAACTAATATATTTTGTACATTTTGTTTTACAAGTCCTAGGAAAGATTGTCTTCTGAAAATTTGATGTCTTCTGGGTTGATGGAGATGGGAAGGGTTCTAGGCCAGAATGTTCACATTTGGAAGACTCTTTCAAATTATAACTGTTGTTACATGTTTGCAGTTTATTCAAGACTGCTGTATACATAGTAGACAAATTAACTCCTTACTTGAAACATCTAGTCTATCTAGATGTTTAGAAGTGCCCGATGTATGTTAAATGTATAGGTAGTAAAATACCACTTTGTAAATATCTTTTTGCTAAAATTCATAGGAAATGCTTTTGGAAATTGAATTGTGAAGCCACCTTTGTGAACAGTATAG'
NORAD6='CCACCTTTGTGAACAGTATAGTAATGTCTATACTTGTTCAATAGTTTAGAGGAGGTAGGAGGGAAGAAATTGCAAAAGGTAATATTACTAGTGTGTTCATACTTGGACATTTTCAGACACCATTTTTCTATATGTTTTGTGCATTTTGTTTTGCTCTGTATATAGTATATATAATGGACAAATAGTCCTAATTTTTCAACATCTAGTCTCTAGATGTTAAAGAGGTTGCCAGTGTATGACAAAGGAGTAAAATTAGCATATTTTGTACACTTTGTGTTGAAATTCGTAGGAAAACTTGTCTTCTGTAAAGACTTTTGCATAGGAATTTGTTTGACCATCTCTAAGCATTACACGTGCCTGTACTTGTCCACTGGATTGAAGGCAGAGAAGGAAGGGAGGAGGGAATGATTCAAGGCCAAAATGGCCACATTTAGAAGATACCTCAGATGATAACCATTGTTATGTGTGTGCAATTTTATTTAACAGTGCTGTGTATGTGGTGGACAAGTTATATGAAATATCTAGTCTTTCTAGATATTTGGAAGTGCTTGATGTATTTAAAAGTGGTAGTAGAATAACACTTTGTAAATAGCTTTTAAAAACTGATGGGAAATGCTGTTTGGAAGTGGAATTGTTGAACCACCTGGGAGGTGGGAGGGAAGAAATTGCAAATGGTGTTTTGCCATTGTTTATTAGAAAATTTCAGCTTAATCCATTGTGTATATGTTACATGCATTTCATTTAACTTTGCTATACTGTATATATTGTATATATAACGGACAAATTAGTCCCGATTTTATAATATCTAGTCTCTAGATATTAAAGAGGTTGCCAATGTATGACAGAAGTAGAGTTAGTAAACTAACACATTTTGTACACTTTGTTAAAATTTGTAGAAAGGCTGTCTTCTGAAAAGGACTTTTGGAAGTGAGATAACATCAGCTCTAAGTGACACGTGCCTATATCCATCAGGTTGGTGGTGGAGAGGAGTTGGAAGGAATGAAGGGTTCTAGACCAGAATGTTCGTATTTAGAAGACACTATCAGATATAACCATTGTTACATGTGTGTAGTTTATTCAACCCTACTGTGTATATAGCGGACAAACTTAAGTCCTTATTTGAAACATCTAGTCTTTCTAGATGTTTAGAAGTGCACAAAGTATGTTAAAAGTAGAGGTAGTAAATAACACATTTTGTAGCTATCCTTTTGATATGAAATATTGTCTTGGAAATTGATCAATTCTCTGAGCAGTACCCATTTTGATATTTGTGCTGGTTCAGGGGGAAGGAGGAGCACAAAGTGCAAAGGGCTTTCTACCAGTGTCCAGTGTGTTTATGAGGAGGCACATTGACCATTGTCCCTTATGTCTGCATTTTCATTTACTGTGCTGTGTATATAGTGTATATAAGCGGACATAGGAGTCCTAATTTACGTCTAGTCGATGTTAAAAAGGTTGCCAGTATATGACAAAAGTAGAATTAGTAAACTACTACATTGAGTACACTTTGTGTTAAAATTCATAGGGAAGACTTCTTAAAAACAAGTGAAATTGTTAAAACCCCCCCTAAGCATTACAGATGGCTTATAGCTGTCCACGGGGTTGGTAGAGGTGGGAAAGGGAAGGGTTCTAGGCCAGAATGTTCCTATTTAGAAGACACTCAAATTACAGTCTGTGTTATGTATGTATACCATTTATTCAATGCTACTGTGTATATAATGGAAAACTTAAGTCCAGTTTGAAACATCTAGTCTTTCTAGGTGTTTAAAAGTGTACAACGGCCTGTCGCAGTGGCGCATGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCAGGCGGATCACGAGGTCAAGAGATCAGGAC'
NORAD7='GAGGTCAAGAGATCAGGACCATCTTGGCCAACATGGTGAAACCCCATCTTTACTAAAAATACAAAAATTAGCTGGTCGTGGTGGTGCCCACCTGTAGCCCCAGTTACTCGAGAGGCTGAGGCAGGAGAATCGCTTGAACTTGGGAGGCGGAAGTTGCAGTGAGCCAAGATCGCACCACTGCACTCCAGCCTGGCGACAGAGCGAGGCTCCGTTTCAAAAAAAAAAGTGCACAATGTAGGTTAACAGTAGAGGGCTTAAGTAACACCCCTCTAAGCATTTGTTTTCAGTACTTCCTAGGAGTGGTTGCATTTGGGAATGGAATTGTTAAAACTTGATGCTTAGGAGCGAATGCAGACTATTCATTGGGTGTTTGGGGTGGGGGAAGGGGGGGTGGGCAGAGGAGGTATGCAGGGAGAGGGGTTCTGTGCTCCTGAGATTAGTTCAGATGGTCTAACCATTGTTCTATATGTGCATTTTAGTTAATATTGTGTATTAAAGGATAAGTCTTAATGCTCAAAGTATGTTAAAAATAGATGTAGTAAATCAGTCCCTTTGTGAATGTCCTTTTGTTAGTTTTTAGGAAGGCCTGTCCTCTGGGAGTGACCTTTATTAGTCCACCCCTTGGAGCTAGACATCCTGTACTTAGTCACGGGGATGGTGGAAGAGGGAGAAGAGGAAGGGTGAAGGGAAGGGCTCTTTGCTAGTATCTCCATATCTAGACGATGGTTTTAGATGATAACCACAGGTCTACAAGAGCGTTTTTAGTAAAGTGCCTGTGTTCATTGTGGACAAAGTTATTATTTTGCAACATCTAAGCTTTACGAATGGGGTGACAACTTATGATAAAAACTAGAGCTAGTGAATTAGCCTATTTGTAAATACCTTTGTTATAATTGATAGGATACATCTTGGACATGGAATTGTTAAGCCACCTCTGAGCAGTGTATGTCAGGACTTGTTCATTAGGTTGGCAGCAGAGGGGCAGAAGGAATTATACAGGTAGAGATGTATGCAGATGTGTCCATATATGTCCATATTTACATTTTGATAGCCATTGATGTATGCATCTCTTGGCTGTACTATAAGAACACATTAATTCAATGGAAATACACTTTGCTAATATTTTAATGGTATAGATCTGCTAATGAATTCTCTTAAAAACATACTGTATTCTGTTGCTGTGTGTTTCATTTTAAATTGAGCATTAAGGGAATGCAGCATTTAAATCAGAACTCTGCCAATGCTTTTATCTAGAGGCGTGTTGCCATTTTTGTCTTATATGAAATTTCTGTCCCAAGAAAGGCAGGATTACATCTTTTTTTTTTTTTTTAGCAGTTTGAGTTGGTGTAGTGTATTCTTGGTTATCAGAATACTCATATAGCTTTGGGATTTTGAATTGGTAAATATTCATGATGTGTGAAAAATCATGATACATACTGTACAGTCTCAGTCCCATAAAATTGGATGTTGTGCCTACACACAGGATCTAGAAGAATATGTCAAACTATAAACTGCTTGTGATTGTGAATGACTTTGTTCTTTGCTTGTGTTTTTCAATTTCCTATAATGCACATACTAACTTTTAAAAAATAAAGGTTATTTTAAAAGCCTGTA'

NORADS = NORAD5,NORAD6,NORAD7 
NORADS_NAME = ('NORAD5','NORAD6','NORAD7')
NORAD_units_fasta = directory + "/NORAD_units.fasta"



# Relative to the NRUs
# This is the position: 45 to 134	(0 to 179 by default)
# This is the position: 47 to 63	(PRE zone)
# ~ st_pos = 47
# ~ end_pos = 63
st_pos = 0
end_pos = 179

# http://genome.ucsc.edu/cgi-bin/hgBlat chrName = 'chr20' chrSt = 36045622 chrEnd = 36050960

chrName = 'chr20'
chrSt = 36045622
chrEnd = 36050960


BigWigFile_list = ('UCSC_data/hg38.phyloP100way.bw','UCSC_data/hg38.phyloP17way.bw','UCSC_data/hg38.phastCons7way.bw','UCSC_data/hg38.phastCons17way.bw','UCSC_data/hg38.phyloP4way.bw','UCSC_data/hg38.phastCons4way.bw')

# Check if it is reverse (1) or forward (0)
reverse = 1
name = 'NORAD'

# ~ pattern = '' # Default
# ~ pattern = 'TAAA' # Sam68
# ~ pattern = 'AATATCTAG' # STEM NRU5 and NRU6
# ~ pattern = 'TGT[AG]TATA' # Pum site UGURUAUA (R = A/G)
pattern = 'CTGT[GA]T[AGT][TC]' # MEME motif find by me #close to TGTATATA


stability_pattern(NORAD_units_fasta,temperatures,pattern,st_pos,end_pos,chrName, name, reverse, chrSt, chrEnd,BigWigFile_list)

# TODO: window correlation
