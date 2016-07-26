#

# /Users/mcintoshc/bin/fastqc/FastQC.app/Contents/Resources/Java/fastqc \
# --outdir


####################################################################################################
#  summarize_NGSQCToolkit_reports.py
#  Developer: Carl McIntosh
#  Email:     mcintoshc@nih.gov
#  Date:      2015-05
####################################################################################################
# This pipeline was developed to run PandaSeq on paired-end reads to generate "single-end reads" as
# specified by PandaSeq's developer Josh Neufeld and can be obtained at GitHub:
#      https://github.com/neufeld/pandaseq/wiki/PANDAseq-Assembler
#
# For this pipeline, we selected rather restrictive parameters for PandaSeq which resulted in < 5%

## Import Section
import sys
import subprocess
import os
import getopt
import ntpath
import glob
from collections import OrderedDict
#davidwheeler check_output had problems
#from subprocess import check_output
import ntpath
import fnmatch

########### Dependency Locations ###################################################################


########### Hardwired Program Variables ############################################################
# graphics_query_dict = {'summary.png':'*summary.png','R1_avgQual.png':'*R1*avgQual.png','R2_avgQual.png':'*R2*avgQual.png','R1_QualRangePerBase.png':'*R1*QualRangePerBase.png','R1_filtered_QualRangePerBase.png':'*R1*filtered_QualRangePerBase.png','R2_QualRangePerBase.png':'*R2*QualRangePerBase.png','R2_filtered_QualRangePerBase.png':'*R2*filtered_QualRangePerBase.png','R1_gcDistribution.png':'*R1*gcDistribution.png','R2_gcDistribution.png':'*R2*gcDistribution.png','R1_qualDistribution.png':'*R1*qualDistribution.png','R2_qualDistribution.png':'*R2*qualDistribution.png','R1_baseCompostion.png':'*R1*baseCompostion.png','R2_baseCompostion.png':'*R2*baseCompostion.png'}
# graphics_lable_dict = {'summary.png':'Summary of QC','R1_avgQual.png':'Per base average quality scores for file 1','R2_avgQual.png':'Per base average quality scores for file 2','R1_QualRangePerBase.png':'Read count (%) per base for different quality ranges for file 1 (Raw)','R1_filtered_QualRangePerBase.png':'Read count (%) per base for different quality ranges for file 1 (Filtered)','R2_QualRangePerBase.png':'Read count (%) per base for different quality ranges for file 2 (Raw)','R2_filtered_QualRangePerBase.png':'Read count (%) per base for different quality ranges for file 2 (Filtered)','R1_gcDistribution.png':'GC content distribution for file 1','R2_gcDistribution.png':'GC content distribution for file 2','R1_qualDistribution.png':'Quality distribution for file 1','R2_qualDistribution.png':'Quality distribution for file 2','R1_baseCompostion.png':'Base composition for file 1','R2_baseCompostion.png':'Base composition for file 2'}

graphics_query_dict = OrderedDict([('html','*html'),('summary.png','*summary.png'),('R1_avgQual.png','*R1*avgQual.png'),('R2_avgQual.png','*R2*avgQual.png'),('R1_QualRangePerBase.png','*R1*QualRangePerBase.png'),('R1_filtered_QualRangePerBase.png','*R1*filtered_QualRangePerBase.png'),('R2_QualRangePerBase.png','*R2*QualRangePerBase.png'),('R2_filtered_QualRangePerBase.png','*R2*filtered_QualRangePerBase.png'),('R1_gcDistribution.png','*R1*gcDistribution.png'),('R2_gcDistribution.png','*R2*gcDistribution.png'),('R1_qualDistribution.png','*R1*qualDistribution.png'),('R2_qualDistribution.png','*R2*qualDistribution.png'),('R1_baseCompostion.png','*R1*baseCompostion.png'),('R2_baseCompostion.png','*R2*baseCompostion.png')])
graphics_lable_dict = OrderedDict([('html','Sample HTML Report'),('summary.png','Summary of QC'),('R1_avgQual.png','Per base average quality scores for file 1'),('R2_avgQual.png','Per base average quality scores for file 2'),('R1_QualRangePerBase.png','Read count (%) per base for different quality ranges for file 1 (Raw)'),('R1_filtered_QualRangePerBase.png','Read count (%) per base for different quality ranges for file 1 (Filtered)'),('R2_QualRangePerBase.png','Read count (%) per base for different quality ranges for file 2 (Raw)'),('R2_filtered_QualRangePerBase.png','Read count (%) per base for different quality ranges for file 2 (Filtered)'),('R1_gcDistribution.png','GC content distribution for file 1'),('R2_gcDistribution.png','GC content distribution for file 2'),('R1_qualDistribution.png','Quality distribution for file 1'),('R2_qualDistribution.png','Quality distribution for file 2'),('R1_baseCompostion.png','Base composition for file 1'),('R2_baseCompostion.png','Base composition for file 2')])
QC_Stats_keys = ['File name','Total number of reads','Total number of HQ reads','Percentage of HQ reads','Total number of bases','Total number of bases in HQ reads','Total number of HQ bases in HQ reads','Percentage of HQ bases in HQ reads','Number of Primer/Adaptor contaminated HQ reads','Total number of HQ filtered reads','Percentage of HQ filtered reads']
Detailed_QC_Stats_keys = ['File name','Minimum read length','Maximum read length','Average read length','Total number of reads','Total number of reads with non-ATGC bases','Percentage of reads with non-ATGC bases','Total number of bases','Total number of HQ bases','Percentage of HQ bases','Total number of non-ATGC bases','Percentage of non-ATGC bases']


keys = ["Sample", "Basic Statistics","Per base sequence quality","Per tile sequence quality","Per sequence quality scores","Per base sequence content","Per sequence GC content","Per base N content","Sequence Length Distribution","Sequence Duplication Levels","Overrepresented sequences","Adapter Content","Kmer Content"]

basic_stats = ["Filename","File type","Encoding","Total Sequences","Sequences flagged as poor quality","Sequence length","%GC"]

########### Start function #########################################################################
def write_html_header(in_html_file):
    HTML_FILE = open(in_html_file, 'w')
    HTML_FILE.write("<!DOCTYPE html>" + "\n")
    HTML_FILE.write("<html lang=\"en\">" + "\n")



    # Style settings
    HTML_FILE.write("<!-- Style used from http://jimmybonney.com/articles/column_header_rotation_css/ -->" + "\n")
    HTML_FILE.write("<style>" + "\n")
    HTML_FILE.write(".table-header-rotated th.row-header { width: auto;}" + "\n")
    HTML_FILE.write(".table-header-rotated td{ width: 40px; border-top: 1px solid #dddddd; border-left: 1px solid #dddddd; border-right: 1px solid #dddddd; vertical-align: middle; text-align: center;}" + "\n")
    HTML_FILE.write(".table-header-rotated th.rotate-45{height: 80px; width: 40px; min-width: 40px; max-width: 40px; position: relative; vertical-align: bottom; padding: 0; font-size: 12px; line-height: 0.8;}" + "\n")
    HTML_FILE.write(".table-header-rotated th.rotate-45 > div{ position: relative; top: 0px; left: 40px; /* 80 * tan(45) / 2 = 40 where 80 is the height on the cell and 45 is the transform angle*/ height: 100%; -ms-transform:skew(-45deg,0deg); -moz-transform:skew(-45deg,0deg); -webkit-transform:skew(-45deg,0deg); -o-transform:skew(-45deg,0deg); transform:skew(-45deg,0deg); overflow: hidden; border-left: 1px solid #dddddd; border-right: 1px solid #dddddd; border-top: 1px solid #dddddd;}" + "\n")
    HTML_FILE.write(".table-header-rotated th.rotate-45 span { -ms-transform:skew(45deg,0deg) rotate(315deg); -moz-transform:skew(45deg,0deg) rotate(315deg); -webkit-transform:skew(45deg,0deg) rotate(315deg); -o-transform:skew(45deg,0deg) rotate(315deg); transform:skew(45deg,0deg) rotate(315deg); position: absolute; bottom: 30px; /* 40 cos(45) = 28 with an additional 2px margin*/ left: -25px; /*Because it looked good, but there is probably a mathematical link here as well*/ display: inline-block; // width: 100%; width: 85px; /* 80 / cos(45) - 40 cos (45) = 85 where 80 is the height of the cell, 40 the width of the cell and 45 the transform angle*/ text-align: left; // white-space: nowrap; /*whether to display in one line or not*/}" + "\n")
    HTML_FILE.write(".pass { height: auto; width: auto; max-width:7%; border=\"1\"; padding:1px; border:1px solid #021a40; }" + "\n")
    HTML_FILE.write("</style>" + "\n")

    HTML_FILE.write("<head>" + "\n")
    HTML_FILE.write("	<meta charset=\"utf-8\" />" + "\n")
    HTML_FILE.write("	<title>NGSQCToolkit Summary Webpage</title>" + "\n")
    HTML_FILE.write("</head>" + "\n")

    HTML_FILE.write("<body>" + "\n")
    HTML_FILE.write("<h1 align=\"center\" >Aggregate report of NGSQCToolkit for Project Name</h1>" + "\n")
    
    
    
    HTML_FILE.close();

def write_html_closing(in_html_file):
    HTML_FILE = open(in_html_file, 'a')
    HTML_FILE.write("</body>" + "\n")
    HTML_FILE.write("</html>" + "\n")

    HTML_FILE.close();

def pull_lines(in_html_file, QC_Statistics_Search_String, maximum_number_lines):
	end_html_table = "</table></td></tr>"

	start_html_row = "<tr><td>"
	split_column_html_table = "</td><td>"
	end_html_row = "</td></tr>"


#					html_table_lines = pull_lines(in_html_file, "<tr><td class=\"head3\">QC statistics</td></tr>", 20)


#	in_html_file = '/Users/mcintoshc/Desktop/NGSQCToolkit/example_NGSQCToolkit_reports/sample_1/output_FM403N.R1.fastq.gz_FM403N.R2.fastq.gz.html'
#	QC_Statistics_Search_String = "<tr><td class=\"head3\">QC statistics</td></tr>"

	## Get Start
#	grep_result = check_output(["grep", "--line-number", QC_Statistics_Search_String, in_html_file])
#davidwheeler non check_output version
	grep_result = os.popen("grep --line-number '{0}' {1}".format(QC_Statistics_Search_String, in_html_file)).read()
	grep_result = grep_result.split(':')
	start = int(grep_result[0]) + 2

	stop = start + maximum_number_lines

	sed_range = str(start) + "," + str(stop) + "p"
#davidwheeler non check_output version
#	grep_result = check_output(["sed", "-n", sed_range, in_html_file])
	grep_result = os.popen("sed -n '{0}' {1}".format(sed_range, in_html_file)).read()

	grep_result = grep_result.replace(start_html_row, "")
	grep_result = grep_result.replace(end_html_row, "")
	grep_result = grep_result.replace(split_column_html_table, "\t")
	grep_result_lines = grep_result.split('\n')


	indices = [index for index, line in enumerate(grep_result_lines) if '</table>' in line]

	grep_result_lines = grep_result_lines[:indices[0]]

	return grep_result_lines

def handle_files(inSEARCH_SUFFIX_FILE_NAME, inSEARCH_DIR, inAGGREGATE_REPORTS_DIR, in_html_file):
	
	end_html_table = "</table></td></tr>"

	start_html_row = "<tr><td>"
	split_column_html_table = "</td><td>"
	end_html_row = "</td></tr>"
	
	directories = []
	for root, directory_names, filenames in os.walk(inSEARCH_DIR):
		for filename in fnmatch.filter(filenames, inSEARCH_SUFFIX_FILE_NAME):
			directories.append(root)
#			directories.append(os.path.join(root, filename))


# 	directories.remove(".DS_Store")
# 	directories.remove("NGSQCToolkit_SUMMARY_REPORT.html")

	print ("Extracting " + str(len(directories)) + " files:")

	HTML_FILE = open(in_html_file, 'a')
	
	HTML_FILE.write("\n\n" + "<div class=\"outer\">" + "\n")
	print ("################################################################################")
	print ("Processing NGSQCToolkit Reports:")
	row_counter = 0
	file_name = "unassigned"

	for directory in directories:
		if os.path.isfile(directory):
			print "It's a file: " + directory			
		elif os.path.isdir(directory):
			row_counter = row_counter + 1
			
			print str(row_counter) + ". " + directory
			
	       	
	       	for key in graphics_query_dict.keys():
 				img_file = ""
 				query = graphics_query_dict[key] 				
 				title = " title=\"" + graphics_lable_dict[key] + "\" "
 				
 				## Handle parsing HTML for QC Statistics and QC Detailed Statistics
 				if graphics_query_dict[key] == "*html":
					found_files = glob.glob(directory + "/" + graphics_query_dict[key])
					last_directory_name = file_name.split ('/')
					
					last_directory_name = last_directory_name[-1]
					HTML_FILE.write("<h2><a name=\"" + last_directory_name + "\">" + str(row_counter) + ". " + last_directory_name + " </a><a href=\"" + found_files[0] +"\" target=\"_blank\"" + "\">View NGSQCToolkit Sample Report</a></h2>" + "\n")
					
					## Handle QC Statistics
					r1_qc_stats = dict.fromkeys(QC_Stats_keys)
					r2_qc_stats = dict.fromkeys(QC_Stats_keys)
					
					html_table_lines = pull_lines(found_files[0], "<tr><td class=\"head3\">QC statistics</td></tr>", 20)

					for html_table_line in html_table_lines:
							html_table_line_array = html_table_line.split('\t')
							html_line_key = html_table_line_array[0]
							if html_line_key in QC_Stats_keys:
								r1_qc_stats[html_line_key] = html_table_line_array[1]
								r2_qc_stats[html_line_key] = html_table_line_array[2]

					########### Print QC Stats table
					HTML_FILE.write("<h3>QC statistics</h3>" + "\n")
					HTML_FILE.write("<tr><td class=\"table table-striped table-header-rotated\"><table width=\"100%\" border=\"0\" class=\"cnt\">" + "\n")
										
					HTML_FILE.write("<table width=90% class=\"table table-striped table-header-rotated\" ><tr>" + "\n")
					HTML_FILE.write("<tr>")
					for qc_stats in QC_Stats_keys:
						HTML_FILE.write("<th class=\"rotate-45\"><div><span title=\"" + qc_stats + "\"><b><font size=\"-2\">" + qc_stats + "</font></b></span></div></th>")

					HTML_FILE.write("</tr>" + "\n")

					## R1 in table
					HTML_FILE.write(start_html_row)
					col = 0
					for qc_stats in QC_Stats_keys:
						col = col + 1
						if col == 1:
							HTML_FILE.write(r1_qc_stats[qc_stats])
						else:
							HTML_FILE.write(split_column_html_table + r1_qc_stats[qc_stats])
					HTML_FILE.write(end_html_row + "\n")

					## R1 in table
					HTML_FILE.write(start_html_row)
					col = 0
					for qc_stats in QC_Stats_keys:
						col = col + 1
						if col == 1:
							HTML_FILE.write(r2_qc_stats[qc_stats])
						else:
							HTML_FILE.write(split_column_html_table + r2_qc_stats[qc_stats])
					HTML_FILE.write(end_html_row + "\n")

 					# end table
 					HTML_FILE.write(end_html_table + "\n<hr>\n")

 					########## Print QC Detail table
					# HandleDetailed_QC_Stats_keys
					r1_qc_stats = dict.fromkeys(Detailed_QC_Stats_keys)
					r2_qc_stats = dict.fromkeys(Detailed_QC_Stats_keys)
					r3_qc_stats = dict.fromkeys(Detailed_QC_Stats_keys)
					r4_qc_stats = dict.fromkeys(Detailed_QC_Stats_keys)

					html_table_lines = pull_lines(found_files[0], "<tr><td class=\"head3\">Detailed QC statistics</td></tr>", 20)

					for html_table_line in html_table_lines:
							html_table_line_array = html_table_line.split('\t')
							html_line_key = html_table_line_array[0]
							if html_line_key in Detailed_QC_Stats_keys:
								r1_qc_stats[html_line_key] = html_table_line_array[1]
								r2_qc_stats[html_line_key] = html_table_line_array[2]
								r3_qc_stats[html_line_key] = html_table_line_array[3]
								r4_qc_stats[html_line_key] = html_table_line_array[4]

 					HTML_FILE.write("<h3>Detailed QC statistics</h3>" + "\n")
					HTML_FILE.write("<tr><td class=\"table table-striped table-header-rotated\"><table width=\"100%\" border=\"0\" class=\"cnt\">" + "\n")
										
					HTML_FILE.write("<table width=90% class=\"table table-striped table-header-rotated\" ><tr>" + "\n")
					HTML_FILE.write("<tr>")
 					for Detailed_QC_Stats in Detailed_QC_Stats_keys:
 						HTML_FILE.write("<th class=\"rotate-45\"><div><span title=\"" + Detailed_QC_Stats + "\"><b><font size=\"-2\">" + Detailed_QC_Stats + "</font></b></span></div></th>")

 					HTML_FILE.write("</tr>" + "\n")

					# R1 in table
					HTML_FILE.write(start_html_row)
					col = 0
					for Detailed_QC_Stats in Detailed_QC_Stats_keys:
						col = col + 1
						if col == 1:
							HTML_FILE.write(r1_qc_stats[Detailed_QC_Stats])
						else:
							HTML_FILE.write(split_column_html_table + r1_qc_stats[Detailed_QC_Stats])
					HTML_FILE.write(end_html_row + "\n")

					# R2 in table
					HTML_FILE.write(start_html_row)
					col = 0
					for Detailed_QC_Stats in Detailed_QC_Stats_keys:
						col = col + 1
						if col == 1:
							HTML_FILE.write(r2_qc_stats[Detailed_QC_Stats])
						else:
							HTML_FILE.write(split_column_html_table + r2_qc_stats[Detailed_QC_Stats])
					HTML_FILE.write(end_html_row + "\n")

					# R3 in table
					HTML_FILE.write(start_html_row)
					col = 0
					for Detailed_QC_Stats in Detailed_QC_Stats_keys:
						col = col + 1
						if col == 1:
							HTML_FILE.write(r3_qc_stats[Detailed_QC_Stats])
						else:
							HTML_FILE.write(split_column_html_table + r3_qc_stats[Detailed_QC_Stats])
					HTML_FILE.write(end_html_row + "\n")

					# R4 in table
					HTML_FILE.write(start_html_row)
					col = 0
					for Detailed_QC_Stats in Detailed_QC_Stats_keys:
						col = col + 1
						if col == 1:
							HTML_FILE.write(r4_qc_stats[Detailed_QC_Stats])
						else:
							HTML_FILE.write(split_column_html_table + r4_qc_stats[Detailed_QC_Stats])
					HTML_FILE.write(end_html_row + "\n")

 					# end table
 					HTML_FILE.write(end_html_table + "\n<hr>\n")
 					HTML_FILE.write("<h3>Graphics</h3>" + "\n")


 				else:
					found_files = glob.glob(directory + "/" + graphics_query_dict[key])
					found_files_count = len(found_files)
					if found_files_count == 1:

						img_file = str(found_files[0])

						HTML_FILE.write("<a href=\"" + img_file +"\" target=\"_blank\"" + title + " >" + "<img class=\"pass\" src=\"" + img_file + "\" /></a>" + "\n")
					
					elif found_files_count > 1:					
						if "QualRangePerBase.png" in graphics_query_dict[key]:
							for try_file in found_files:
								if not "filtered" in try_file:
									img_file = try_file
									HTML_FILE.write("<a href=\"" + img_file +"\" target=\"_blank\"" + title + " >" + "<img class=\"pass\" src=\"" + img_file + "\" /></a>" + "\n")
 						else:
 							print "ERROR: Some need to handle this: " + graphics_query_dict[key] 
	HTML_FILE.write("<span class=\"fix\"></span></div>" + "\n")

########### Main ###################################################################################
## Credit modified from https://docs.python.org/2/library/getopt.html

def main():
	# These should be parameters
	SEARCH_SUFFIX_FILE_NAME = sys.argv[1]
	SEARCH_DIR = sys.argv[2]
	AGGREGATE_REPORTS_DIR = sys.argv[3]
		
	html_report = AGGREGATE_REPORTS_DIR + "/NGSQCToolkit_SUMMARY_REPORT.html" 
	subprocess.call(["rm", html_report])

	write_html_header(html_report)
	handle_files(SEARCH_SUFFIX_FILE_NAME, SEARCH_DIR, AGGREGATE_REPORTS_DIR, html_report)
	
	write_html_closing(html_report)


if __name__ == "__main__":
    main()







