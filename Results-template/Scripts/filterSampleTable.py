##############################################################################################
# Description:  Creates a new sampletable.txt for differential expression analysis which
#		is necessary because a user may remove samples from their groups.tab file
#		after running QC. The new sampletable.txt file is created from the
#		groups.tab file that is provided when running DE half of the pipeline.
#		Also a corresponding RawCounts_RSEM_genes.txt file is created which only
#		contains information for samples that that are specificed in groups.tab.
#		This file is written to DEG_{group1}-{group2}_{mincpm}_{minsample} folder.
#		It is necessary that the order of the groups listed in sampletable.txt and
#		and RawCounts_RSEM_genes_filtered.txt match. If they do not the reported
#		direction (sign: postive or negative) of each gene's fold-change will be
#		opposite of what it should be. This program also addresses that problem.
# USAGE:	python filterSampleTable.py -s 'rbfoxKO_1,rbfoxKO_2,control_1,control_2'\
#					    -g 'KO,KO,WT,WT' -l 'KO_1,KO_2,WT_1,WT_2' \
#					    -r 'RawCountFile_RSEM_genes.txt' \
#					    -outsample 'outfolder/sampletable.txt'\
#					    -outraw 'outfolder/RawCountFile_RSEM_genes.txt'
##############################################################################################

from __future__ import print_function, division
import os, sys
import argparse  # module load python/3.5




def filteredRawCounts(rawcountsfilename, samples, outrawfilename):
	'''Generator that yields rawcounts information for samples in the current DEG groups.tab file.
	Helper function checks formating of the headers between the newly generated sampletable.txt and
	input RawCounts_RSEM_genes.txt file. If the order of groups formatting differs, it fixes it to match
	the ordering in sampletable.txt (fail-safe to ensure later referential integrity at DE).
	'''
	def formatted(headerlist, samples, linelist):
		headerIndexes = {}
		symbol = linelist[0]
		# Create Dictionary for mapping samplename in header of RawCounts_RSEM_genes.txt to its position 
		for i in range(len(headerlist)):
			headerIndexes[headerlist[i]] = i
		indexlist = []
		formattedlist = [symbol]
		for sample in samples:
			try:
				indexlist.append(headerIndexes[sample])
				formattedlist.append(linelist[headerIndexes[sample]])
			except KeyError:
				raise Exception('Key Error: Failed to find sample {} in RawCounts_RSEM_genes.txt.\nDid you add an additional sample to groups.tab?'.format(sample))
		return formattedlist

	rawfh = open(rawcountsfilename, 'r')
	headerlist = [ fn.split('/')[-1].split('.')[0] for fn in next(rawfh).strip().split('\t')]

	newheader = '\t'.join([headerlist[0]]+ samples) + '\n'
	print(newheader)
	outfh = open(outrawfilename, 'w')
	outfh.write(newheader)
	for line in  rawfh:
		linelist = line.strip().split('\t')
		formattedlist = formatted(headerlist, samples, linelist)
		#print(formattedlist)
		outfh.write('\t'.join(formattedlist) + '\n')
	rawfh.close()
	outfh.close()

def SampleTable(samples, groups, labels,contrasts):
	'''Generator that yields sample, group, and label info from the current DEG groups.tab file'''
	for i in range(len(samples)):
		if groups[i] in contrasts:
			yield samples[i], groups[i], labels[i]



if __name__ == '__main__':

	# USAGE: python filterSampleTable.py -c 'KO WT' -s 'rbfoxKO_1,rbfoxKO_2,rbfoxKO_3,control_1,control_2,control_3' -g 'KO,KO,KO,WT,WT,WT' -l 'KO_1,KO_2,KO_3,WT_1,WT_2,WT_3' -r 'RawCountFile_RSEM_genes.txt' -outsample 'outfolder/sampletable.txt' -outraw 'outfolder/RawCountFile_RSEM_genes.txt'
	# Parse Command-line arguments
	parser = argparse.ArgumentParser(description='Filters sampletable.txt, cannot always use the sampletable that is generated in QC, sometimes DE is run with less samples.')
	parser.add_argument('-r','--inputRawCounts', type=str, required=True, help='Input RawCounts_RSEM_genes.txt file that is generated QC: (required)')
	parser.add_argument('-outraw','--outputRawCounts', type=str, required=True, help='Output RawCounts_RSEM_genes.txt filename (required)')
	parser.add_argument('-outsample','--outputSampleTable', type=str, required=True, help='Output RawCounts_RSEM_genes.txt filename (required)')
	parser.add_argument('-s','--samples', type=str, required=True, help='samples string from the run.json file  (required)')
	parser.add_argument('-c','--contrast', type=str, required=True, help='contrast string (required)')
	parser.add_argument('-g','--groups', type=str, required=True, help='groups string from the run.json file (required)')
	parser.add_argument('-l','--labels', type=str, required=True, help='labels string from the run.json file (required)')
	args = parser.parse_args()

	sampleslist = args.samples.split(',')
	groupslist = args.groups.split(',')
	labelslist = args.labels.split(',')
	contrastlist = args.contrast.split(' ')
	print(sampleslist)
	print(groupslist)
	print(labelslist)


	#Creating the new sampletable.txt file based off of the current DEG groups.tab information
	samptablefh = open(args.outputSampleTable, 'w')
	samptablefh.write('sampleName\tfileName\tcondition\tlabel\n')
	samplesincontrast = []
	groupset = []
	for sample, group, label in SampleTable(sampleslist, groupslist, labelslist, contrastlist):
		print('{0}\t{0}\t{1}\t{2}'.format(sample, group, label))
		samplesincontrast.append(sample)
		groupset.append(group)
		samptablefh.write('{0}\t{0}\t{1}\t{2}\n'.format(sample, group, label))
	samptablefh.close()

	if len(set(groupset)) < 2:
		raise Exception('Error in groups.tab file. Contrast does not contain two groups!')

	#Creating new RawCounts_RSEM_genes.txt file based off of the current DEG groups.tab information
	rawfn = args.inputRawCounts
	outrawfn = args.outputRawCounts
	filteredRawCounts(rawfn, samplesincontrast, outrawfn)


