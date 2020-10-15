import pandas as pd 
import numpy.ma as ma
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
from string import ascii_letters
import seaborn as sns
import os.path
from pathlib import Path
import itertools    
import re
import os

from string import ascii_letters

def get_NORAD_score(target,position,temp):
	file_open = directory + "/" +  target +  "_" + str(temp) +  ".csv"

	NORAD_score  = pd.read_csv(file_open,sep=';',header=None,index_col = 0) 

	# Remove last position
	return NORAD_score.iloc[:, position:position + 179]

def get_NORAD_full_score(target,temp):
	file_open = directory + "/" +  target +  "_" + str(temp) +  ".csv"

	NORAD  = pd.read_csv(file_open,sep=';',header=None,index_col = 0) 

	# Remove last position
	return NORAD

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
	
	output_file = dir_pattern_plot + MOLECULE + '.pdf'
	figure = ax.get_figure()    
	figure.savefig(output_file, dpi=400)
	figure.clf()
	
	plt.close(f)
	plt.close('all')
	

def plot_box_plot(df_pos,df_neg,MOLECULE,BigWigFile_list,df_phastCons_all):
	
	# Figure
	sns.set_style("whitegrid") 
	
	df_pos['type'] = 'NRU_pattern'
	
	df_neg = df_neg
	df_neg['type'] = 'NRU_rest'


	if not (df_phastCons_all.empty):
		df_phastCons_all = df_phastCons_all.T
		df_phastCons_all.columns=list(BigWigFile_list)
		df_phastCons_all['type'] = 'whole_NORAD'
		# Drop first row
		df_phastCons_all = df_phastCons_all.iloc[1:]
		
		frames = [df_pos, df_neg,df_phastCons_all]
	else:
		frames = [df_pos, df_neg]
		
	df = pd.concat(frames,sort=False)


	# Default regular plot
	df.boxplot(by='type',column=list(BigWigFile_list), grid = False,figsize=(12,9))
	pdf_out = dir_pattern_plot + MOLECULE + '_boxplot.pdf'
	plt.savefig(pdf_out, dpi=800)
	plt.clf()
	
	
def plot_box_plot_final(df_pos,df_phastCons_all,name,fields):
	
	# Figure
	sns.set_style("whitegrid") 
	
	
	####
	## PhastCons 
	
	# Set values for phastCons and phyloP
	df_pos_cons = df_pos[fields]
	df_pos_cons['type'] = 'NRU_pattern'
	

	# Set values for phastCons and phyloP
	df_phastCons_all = df_phastCons_all.T
	df_phastCons_all.columns=list(fields)
	df_phastCons_all['type'] = 'whole_NORAD'

	# Drop first row
	df_phastCons_all = df_phastCons_all.iloc[1:]


	frames = [df_pos_cons, df_phastCons_all]
	df = pd.concat(frames,sort=False)
	
	
	
	df.boxplot(by='type',column=list(fields), grid = False,figsize=(12,9))
	
	pdf_out = dir_pattern_plot + name + '_boxplot.pdf'
	plt.savefig(pdf_out, dpi=800)
	plt.clf()
	
	save = dir_pattern_csv + name + '_boxplot.csv'
	df.to_csv(save,index=False)
	
def plot_box_plot_final_field(df_pos,df_neg,whole_NORAD_df_temp,name,fields):
	
	# Figure
	sns.set_style("whitegrid") 
	

	####
	## PhastCons 
	
	# Set values for phastCons and phyloP
	df_pos_cons = df_pos[fields]
	df_pos_cons['type'] = 'NRU_pattern'

	
	# ~ df_neg = df_neg.T
	df_neg = df_neg[fields]
	df_neg['type'] = 'NRU_rest'

	# Drop first row
	df_neg = df_neg.iloc[1:]
	
	# Set value for norad
	whole_NORAD_df_temp = whole_NORAD_df_temp[fields]
	whole_NORAD_df_temp['type'] = 'whole_NORAD'

	frames = [df_pos_cons, df_neg,whole_NORAD_df_temp]
	df = pd.concat(frames,sort=False)
	df.boxplot(by='type',column=list(fields), grid = False,figsize=(12,9))
	
	pdf_out = dir_pattern_plot + name + '_boxplot.pdf'
	plt.savefig(pdf_out, dpi=800)
	plt.clf()
	
	save = dir_pattern_csv + name + '_boxplot.csv'
	df.to_csv(save,index=False)

def get_full_NORAD(temp):
	new_name = 'NORAD_' + str(temp)
	NORAD_full_1 = get_NORAD_full_score(NORADS_NAME[0],temp)
	len_NORAD_full_1 = len(NORAD_full_1.columns)
	NORAD_full_1.columns = np.arange(1,len_NORAD_full_1+1)
	NORAD_full_1.rename(index={ NORAD_full_1.index[0]: new_name },inplace=True)


	NORAD_full_2 = get_NORAD_full_score(NORADS_NAME[1],temp)
	len_NORAD_full_2 = len_NORAD_full_1 + len(NORAD_full_2.columns)
	NORAD_full_2.columns = np.arange(len_NORAD_full_1,len_NORAD_full_2)
	NORAD_full_2.rename(index={ NORAD_full_2.index[0]: new_name },inplace=True)

	
	NORAD_full_3 = get_NORAD_full_score(NORADS_NAME[2],temp)
	len_NORAD_full_3 = len_NORAD_full_2+ len(NORAD_full_3.columns)
	NORAD_full_3.columns = np.arange(len_NORAD_full_2,len_NORAD_full_3)
	NORAD_full_3.rename(index={ NORAD_full_3.index[0]: new_name },inplace=True)

	frames = [NORAD_full_1.T,NORAD_full_2.T,NORAD_full_3.T]
	df_temp = pd.concat(frames)
	return df_temp
	
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
	data_all_units_positive_corrected = []
	data_all_units_negative = []
	
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
			
			####################
			# Full NORAD (5 6 and 7)
			####################
			
			
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
			csv_file = dir_pattern_csv +  MOLECULE + '_' +  str(st_pos) + "_" + str(end_pos) +  "_UNIT.csv" 
			final_df.to_csv(csv_file)
						
			csv_file = dir_pattern_csv  +  MOLECULE + '_' +  str(st_pos) + "_" + str(end_pos) +  "_UNIT_corr.csv" 
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
					df_positive_corrected = final_df[final_df.columns.intersection(positive_col)]
					df_negative = final_df[final_df.columns.intersection(negative_col)]
					
					MOLECULE_pos = MOLECULE + "_" + pattern
					MOLECULE_neg = MOLECULE + "_NO_" + pattern
					plot_correlation_matrix(df_positive,MOLECULE_pos)
					plot_correlation_matrix(df_negative,MOLECULE_neg)
					
					df_phastCons_all = pd.concat(df_phastCons_original)
					
					# Save to csv # Only pattern positions
					csv_file = dir_pattern_csv + MOLECULE_pos + ".csv" 
					df_positive.to_csv(csv_file)
					df_positive = df_positive.T

					data_all_units_positive.append(df_positive)
					data_all_units_positive_corrected.append(df_positive)
						
					# TODO: Check this
					k=1
					# Rename columns
					df_negative.columns = np.arange(k,len(df_negative.columns)+1)
					k = k+len(df_negative.columns)

					
					
					# ~ df_phastCons[BigWigFile].index = list(NORAD_score_df.index)
					# ~ df_phastCons[BigWigFile].columns=idx
					df_negative = df_negative.T
					data_all_units_negative.append(df_negative)
					
					
					MOLECULE_stdev = MOLECULE + "_stdev_" + pattern
					field = ['NRU_std']
					df_empty = pd.DataFrame({'A' : []})
					
			
					plot_box_plot(df_positive,df_negative,MOLECULE_pos,BigWig_no_path,df_phastCons_all)
					plot_box_plot(df_positive,df_negative,MOLECULE_stdev,field,df_empty)
			else:
				print("No pattern was given")
			
			# Todo: remove this
			# ~ return
			
			line = fp.readline()
		
		
		# Get full NORAD and Reindex

		df_temp1 = get_full_NORAD(temp1)
		df_temp2 = get_full_NORAD(temp2)
		df_temp3 = get_full_NORAD(temp3)
		frames = [df_temp1.T,df_temp2.T,df_temp3.T]
		
		whole_NORAD_df_temp = pd.concat(frames)
		whole_NORAD_df_temp = whole_NORAD_df_temp.T
				
		# Get dataframe and calculate mean and stdev
		# I call it NRU_std but in fact is whole_NRU but i need this name for the boxplot
		whole_NORAD_df_temp['NRU_avg'] = whole_NORAD_df_temp.mean(axis=1)
		whole_NORAD_df_temp['NRU_std'] = whole_NORAD_df_temp.std(axis=1)
		
		# Rename columns
		k = 0
		p = 0
		for dup in data_all_units_positive:
			
			p = k + len(dup.index)
			dup.index = np.arange(k,p)
			k = k + p
		
		# Rename columns
		k = 0
		p = 0
		for dup in data_all_units_negative:
			
			p = k + len(dup.index)
			dup.index = np.arange(k,p)
			k = k + p
		
		df_all_units_positive = pd.concat(data_all_units_positive,sort=False)
		df_all_units_positive_corrected = pd.concat(data_all_units_positive,sort=False)
		df_all_units_negative = pd.concat(data_all_units_negative,sort=False)
		
		name = "all_cons_" + pattern 
		plot_box_plot_final(df_all_units_positive,df_phastCons_all,name,BigWig_no_path)
		
		name = "all_stdev_" + pattern 
		field = ['NRU_std']
		
		df_all_units_positive.to_csv('positive',index=False)
		df_all_units_negative.to_csv('negative',index=False)
		
		# df_all_units_positive and df_all_units_negative are shifted but i don't care because then i'll keep only NORAD_std
		plot_box_plot_final_field(df_all_units_positive,df_all_units_negative,whole_NORAD_df_temp,name,field)
		
		# Plot correlation altogether
		for dup in data_all_units_positive_corrected:
			dup.rename(columns={ dup.columns[5]: "NRU_23" },inplace=True)
			dup.rename(columns={ dup.columns[6]: "NRU_37" },inplace=True)
			dup.rename(columns={ dup.columns[7]: "NRU_55" },inplace=True)
			dup.rename(columns={ dup.columns[0]: "NORAD_23" },inplace=True)
			dup.rename(columns={ dup.columns[1]: "NORAD_37" },inplace=True)
			dup.rename(columns={ dup.columns[2]: "NORAD_55" },inplace=True)
		

		df_all_units_positive_corrected  = pd.concat(data_all_units_positive_corrected,sort=False)
		df_all_units_positive_corrected = df_all_units_positive_corrected.drop(labels='type', axis=1) # axis 1 drops columns, 0 will drop rows that match index value in labels

		name = "all_cons_" + pattern 
		save = dir_pattern_csv + "all_cons_" + pattern + '.csv'
		df_all_units_positive_corrected.to_csv(save,index=False)
	
		plot_correlation_matrix(df_all_units_positive_corrected.T,name)

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
pattern = 'TAAA' # Sam68
# ~ pattern = 'TGT[AG]TATA' # Pum site UGURUAUA (R = A/G)
# ~ pattern = 'AATATCTAG' # STEM NRU5 and NRU6
# ~ pattern = 'CTGT[GA]T[AGT][TC]' # MEME motif find by me #close to TGTATATA
# ~ pattern = 'TAGA' # RANDOM motif


# Define directories
dir_pattern = 'out_NORAD/' + pattern + '/'
Path(dir_pattern).mkdir(parents=True, exist_ok=True)

dir_pattern_csv = dir_pattern + '/csv/'
Path(dir_pattern_csv).mkdir(parents=True, exist_ok=True)

dir_pattern_plot = dir_pattern + '/plot/'
Path(dir_pattern_plot).mkdir(parents=True, exist_ok=True)
				
stability_pattern(NORAD_units_fasta,temperatures,pattern,st_pos,end_pos,chrName, name, reverse, chrSt, chrEnd,BigWigFile_list)

# TODO: window correlation
