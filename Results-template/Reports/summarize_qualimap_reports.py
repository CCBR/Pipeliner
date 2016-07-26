#

# /Users/mcintoshc/bin/fastqc/FastQC.app/Contents/Resources/Java/fastqc \
# --outdir


####################################################################################################
#  summarize_qualimap_reports.py
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
#not found dwheeler
#from subprocess import check_output
import subprocess
import ntpath
import fnmatch

########### Dependency Locations ###################################################################


########### Hardwired Program Variables ############################################################
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
    HTML_FILE.write(".pass { height: auto; width: auto; max-width:15%; border=\"1\"; padding:1px; border:1px solid #021a40; }" + "\n")
    HTML_FILE.write("</style>" + "\n")

    HTML_FILE.write("<head>" + "\n")
    HTML_FILE.write("        <meta charset=\"utf-8\" />" + "\n")
    HTML_FILE.write("        <title>Qualimap Summary Webpage</title>" + "\n")
    HTML_FILE.write("</head>" + "\n")

    HTML_FILE.write("<body>" + "\n")
    HTML_FILE.write("<h1 align=\"center\" >Aggregate report of Qualimap for Project Name</h1>" + "\n")
    
    
    
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


#                                        html_table_lines = pull_lines(in_html_file, "<tr><td class=\"head3\">QC statistics</td></tr>", 20)


#        in_html_file = '/Users/mcintoshc/Desktop/NGSQCToolkit/example_NGSQCToolkit_reports/sample_1/output_FM403N.R1.fastq.gz_FM403N.R2.fastq.gz.html'
#        QC_Statistics_Search_String = "<tr><td class=\"head3\">QC statistics</td></tr>"

        ## Get Start
        grep_result = subprocess.check_output(["grep", "--line-number", QC_Statistics_Search_String, in_html_file]).decode('utf-8')
#dwheeler replace check_output
#        g= os.popen("grep --line-number {0} {1}".format(QC_Statistics_Search_String, in_html_file))
#        grep_result=g.read()
#        g.close()                       
        grep_result = grep_result.split(':')
        start = int(grep_result[0]) + 2

        stop = start + maximum_number_lines

        sed_range = str(start) + "," + str(stop) + "p"
#dwheeler
        grep_result = subprocess.check_output(["sed", "-n", sed_range, in_html_file]).decode('utf-8')
#        g= os.popen("sed -n {0} {1}".format(sed_range, in_html_file))
#        grep_result=g.read()
#        g.close()                       
        grep_result = grep_result.replace(start_html_row, "")
        grep_result = grep_result.replace(end_html_row, "")
        grep_result = grep_result.replace(split_column_html_table, "\t")
        grep_result_lines = grep_result.split('\n')


        indices = [index for index, line in enumerate(grep_result_lines) if '</table>' in line]

        grep_result_lines = grep_result_lines[:indices[0]]

        return grep_result_lines

def handle_files(inSEARCH_SUFFIX_FILE_NAME, inSEARCH_DIR, inAGGREGATE_REPORTS_DIR, in_html_file):
        sample_name = None
        end_html_table = "</table></td></tr>"

        start_html_row = "<tr><td>"
        split_column_html_table = "</td><td>"
        end_html_row = "</td></tr>"
        
        qualimap_html_reports = []
        for root, directory_names, filenames in os.walk(inSEARCH_DIR):
                for filename in fnmatch.filter(filenames, inSEARCH_SUFFIX_FILE_NAME):
                        qualimap_html_reports.append(os.path.join(root, filename))

        print ("################################################################################")
        print ("Extracting " + str(len(qualimap_html_reports)) + " files:")

        HTML_FILE = open(in_html_file, 'a')
        
        HTML_FILE.write("\n\n" + "<div class=\"outer\">" + "\n")
        print ("################################################################################")
        print ("Processing Qualimap Reports:")
        row_counter = 0
        file_name = "unassigned"

        for qualimap_html_report in qualimap_html_reports:
                row_counter = row_counter + 1
                print (str(row_counter) + ". " + qualimap_html_report)
                qualimap_img_dir = os.path.dirname(qualimap_html_report) + "/images_qualimapReport"
                s_string="<td class=column"
                grep_result = subprocess.check_output(["grep", s_string, qualimap_html_report]).decode('utf-8')
#                g = os.popen("grep '<td class=column' {0}".format(qualimap_html_report))
#                grep_result=g.read()
#                g.close()
#                print(grep_result)[1:25]
                
                grep_result = grep_result.replace("<td class=column1>", "")
                grep_result = grep_result.replace("</td>", "")
                grep_result = grep_result.replace("\n<td class=column2>", "\t")
                grep_result = grep_result.replace(": ","")

                grep_result_array = grep_result.split("\n")

                keys = []
                values = []
                cmd_line_value = None
                for grep_result in grep_result_array:
                        line = grep_result.split('\t')
                        if len(line) == 1:
                                if "qualimap bamqc -bam" in line[0]:
                                        cmd_line_value = line[0]
                                
                        if len(line) == 2:
                                (key, value) = grep_result.split('\t')
                                key = key.lstrip()
                                key = key.rstrip()
                                value = value.lstrip()
                                value = value.rstrip()
                                if "BAM file" in key:
                                        value, file_extension = os.path.splitext(os.path.basename(value))
                                        sample_name = value
#                                        print("here is the "+sample_name)
                                if "Analysis date" not in key:
                                        keys.append(key)
                                        values.append(value)
                
                if row_counter == 1:
        #                HTML_FILE.write("<a name=\"STATS_TABLE\"></a><h2>Basic Statistics and QC Test Status Table</h2>\n")
                        HTML_FILE.write("<a name=\"TABLE\"></a>" + "<table width=90% class=\"table table-striped table-header-rotated\" ><tr>" + "\n")

                        for aKey in keys:
                                HTML_FILE.write("<th class=\"rotate-45\"><div><span title=\"" + aKey + " \"><b><font size=\"-2\">" + aKey + "</font></b></span></div></th>\n")
                        HTML_FILE.write("</tr>")
                
                ## R1 in table
                HTML_FILE.write(start_html_row)
                col = 0
                
                for value in values:
                        col = col + 1
                        if col == 1:
                                HTML_FILE.write("<a href=\"#" + value + "\">"+ value +"</a>")
                                
                                #<font size=\"-1\">     "</font>"
                                
                        else:
                                HTML_FILE.write(split_column_html_table + "<font size=\"-1\">" + value + "</font>")
                                
                                
                                
                                
                        
        HTML_FILE.write(end_html_row + "</table>\n<br>\n")
                # Grahpics
        for qualimap_html_report in qualimap_html_reports:
#                print("here is report "+qualimap_html_report)
#                HTML_FILE.write("<hr><br><h3><a name=\"" + sample_name + "\"><h3>" + sample_name + " [<a href=\"" + qualimap_html_report +" \" target=\"_blank\" " + "\">Full Qualimap Report</a>, <a href=\"#TABLE\">Table</a>]</h3>")                
#dwheeler
                sample_name=qualimap_html_report.split("/")[1]
                print("Sample name is: {0}".format(sample_name))
                sample_name=sample_name.split(".")[0]+"."+sample_name.split(".")[1]
                sample_name=str(sample_name)
#                print("here is the second "+sample_name)
                HTML_FILE.write("<hr><br><h3><a name=\"" + sample_name + "\"><h3>" + sample_name + " [<a href=\"" + qualimap_html_report +" \" target=\"_blank\" " + "\">Full Qualimap Report</a>, <a href=\"#TABLE\">Table</a>]</h3>")                
                
                ## Graphics
                qualimap_img_dir = os.path.dirname(qualimap_html_report) + "/images_qualimapReport"
                if os.path.isdir(qualimap_img_dir):
#                        print ("here is " + qualimap_img_dir)
                        qualimap_graphics = []
                        for root, directory_names, filenames in os.walk(qualimap_img_dir):
                                for filename in fnmatch.filter(filenames, "*.png"):
                                        qualimap_graphics.append(os.path.join(root, filename))
                                        img_file = os.path.join(root, filename)
                                        title = filename
#                                        print("here is "+img_file+" and "+title)                         
                                
                                
                                        HTML_FILE.write("<a href=\"" + img_file +"\" target=\"_blank\"" + title + " >" + "<img class=\"pass\" src=\"" + img_file + "\" /></a>" + "\n")
                else:
                        print ("ERROR: images_qualimapReport directory not found beside qualimap html report.")
                        print ("       Program will exit.")
                        sys.exit()
                        
                        
#        sys.exit()

                        
                        
                        
                        
                        
                        
                        
                        






## Usage:
# clear
# python /Users/mcintoshc/Desktop/QC_REPORTS_CEM/SCRIPTS/summarize_qualimap_reports.py \
# qualimapReport.html \
# /Users/mcintoshc/Desktop/QC_REPORTS_CEM/qualimap_reports \
# /Users/mcintoshc/Desktop/QC_REPORTS_CEM/aggregate_reports






########### Main ###################################################################################

def main():
        # These should be parameters
        SEARCH_SUFFIX_FILE_NAME = sys.argv[1]
        SEARCH_DIR = sys.argv[2]
        AGGREGATE_REPORTS_DIR = sys.argv[3]
                
        html_report = AGGREGATE_REPORTS_DIR + "/QUALIMAP_SUMMARY_REPORT.html" 
        subprocess.call(["rm", html_report])

        write_html_header(html_report)
        handle_files(SEARCH_SUFFIX_FILE_NAME, SEARCH_DIR, AGGREGATE_REPORTS_DIR, html_report)
        
        write_html_closing(html_report)


if __name__ == "__main__":
    main()

