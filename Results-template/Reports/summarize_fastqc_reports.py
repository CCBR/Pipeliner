#

# /Users/mcintoshc/bin/fastqc/FastQC.app/Contents/Resources/Java/fastqc \
# --outdir


####################################################################################################
#  summarize_fastqc_reports_gz.py
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
import fnmatch

########### Dependency Locations ###################################################################


########### Hardwired Program Variables ############################################################
keys = ["Sample", "Basic Statistics","Per base sequence quality","Per tile sequence quality","Per sequence quality scores","Per base sequence content","Per sequence GC content","Per base N content","Sequence Length Distribution","Sequence Duplication Levels","Overrepresented sequences","Adapter Content","Kmer Content"]
extract_files = ["summary.txt", "fastqc_data.txt", "per_base_quality.png","per_tile_quality.png","per_sequence_quality.png", "per_base_sequence_content.png","per_sequence_gc_content.png","per_base_n_content.png","sequence_length_distribution.png", "duplication_levels.png" ,"adapter_content.png","kmer_profiles.png" ]
file_key_dict = {'per_base_quality.png':'Per base sequence quality', 'per_tile_quality.png':'Per tile sequence quality', 'per_sequence_quality.png':'Per sequence quality scores', 'per_base_sequence_content.png':'Per base sequence content', 'per_sequence_gc_content.png':'Per sequence GC content', 'per_base_n_content.png':'Per base N content', 'sequence_length_distribution.png':'Sequence Length Distribution', 'duplication_levels.png':'Sequence Duplication Levels', 'adapter_content.png':'Adapter Content', 'kmer_profiles.png':'Kmer Content' } 
basic_stats = ["Filename","File type","Encoding","Total Sequences","Sequences flagged as poor quality","Sequence length","%GC"]

SUPPLIMENTAL_QC_DIR_NAME = "suppliment_fastqc_dir"


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
    HTML_FILE.write(".pass { width:120px; height=90px; border=\"3\"; padding:1px; border:5px solid #021a40; border-color:#00ff00; }" + "\n")
    HTML_FILE.write(".fail { width:120px; height=90px; border=\"3\"; padding:1px; border:5px solid #021a40; border-color:#FF0000; }" + "\n")
    HTML_FILE.write(".warn { width:120px; height=90px; border=\"3\"; padding:1px; border:5px solid #021a40; border-color:#FFFF00; }" + "\n")
    HTML_FILE.write("</style>" + "\n")

    HTML_FILE.write("<head>" + "\n")
    HTML_FILE.write("	<meta charset=\"utf-8\" />" + "\n")
    HTML_FILE.write("	<title>FastQC Summary Webpage</title>" + "\n")
    HTML_FILE.write("</head>" + "\n")

    HTML_FILE.write("<body>" + "\n")
    HTML_FILE.write("<h1 align=\"center\" >Aggregate report of fastq files for Project Name</h1>" + "\n")
    
    
    
    HTML_FILE.close();

def write_html_closing(in_html_file):
    HTML_FILE = open(in_html_file, 'a')
    HTML_FILE.write("</body>" + "\n")
    HTML_FILE.write("</html>" + "\n")

    HTML_FILE.close();

########### dddd()
        

def extract_fastqc_zip_files(inSEARCH_SUFFIX_FILE_NAME, inSEARCH_DIR, inSUPPLIMENTAL_QC_DIR):
    ## http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    fastqc_ziped_files = []
    for root, directory_names, filenames in os.walk(inSEARCH_DIR):
        for filename in fnmatch.filter(filenames, '*' + inSEARCH_SUFFIX_FILE_NAME):
            fastqc_ziped_files.append(os.path.join(root, filename))

    print ("Extracting " + str(len(fastqc_ziped_files)) + " files:")

    for inZipFastqcFile in fastqc_ziped_files:
    	print "\t" + inZipFastqcFile
        file_name = ntpath.basename(inZipFastqcFile)
        file_name = file_name.split(".")[0]

        ## Create Acrhive manifest to locate two files to parse
        archive_manifest_tmp= inSUPPLIMENTAL_QC_DIR + "/" + "archive_manifest_tmp.txt"
		
        TMP_MANIFEST_FILE = open(archive_manifest_tmp, 'w')
        subprocess.call(["unzip", "-l", inZipFastqcFile], stdout=TMP_MANIFEST_FILE)
        TMP_MANIFEST_FILE.close()
        for extract_file in extract_files:
        	search_term = "[^ ]*\/" + extract_file + "[^ ]*"
        	grep_result = subprocess.Popen(["grep", "--only-matching", search_term, archive_manifest_tmp], stdout=subprocess.PIPE)
        	extract_file_loc = grep_result.stdout.read()
        	extract_file_loc = extract_file_loc.rstrip()
        	extract_file = inSUPPLIMENTAL_QC_DIR + "/" + file_name + "." + extract_file
		
        	TMP_FILE = open(extract_file, 'w')
        	subprocess.call(["unzip", "-p", inZipFastqcFile, extract_file_loc], stdout=TMP_FILE)
        	TMP_FILE.close();

        ## Clean extra files
        subprocess.call(["rm", archive_manifest_tmp])

    return fastqc_ziped_files

def process_summary_file(in_html_file, in_R_script_file, in_fastqc_ziped_files, in_SUPPLIMENTAL_QC_DIR):
    print ("Processing " + str(len(in_fastqc_ziped_files)) + " summary files:")
    HTML_FILE = open(in_html_file, 'a')

    HTML_FILE.write("<div style=\"color=black\">" + "\n")

    R_SCRIPT_FILE = open(in_R_script_file, 'a')
    R_SCRIPT_FILE.write("## Create Matrix\n" + "\n")


    ## Table
    HTML_FILE.write("<a name=\"STATS_TABLE\"></a><h2>Basic Statistics and QC Test Status Table</h2>\n")
    HTML_FILE.write("<table width=90% class=\"table table-striped table-header-rotated\" ><tr>" + "\n")

    for aKey in ["", "Sample"] + basic_stats[3:] + keys[1:]:
        HTML_FILE.write("<th class=\"rotate-45\"><div><span title=\"" + aKey + " \"><b><font size=\"-2\">" + aKey + "</font></b></span></div></th>")

        if aKey == "Sample":
            R_SCRIPT_FILE.write("keys<-c(\"" + aKey + "\"")
        else:
            R_SCRIPT_FILE.write(",\"" + aKey + "\"")

    HTML_FILE.write("</tr>")
    R_SCRIPT_FILE.write(")\r")
    first_row_name = ""

#    R_SCRIPT_FILE.write("DDDDD" + "\n")

    row_count = 0
    row_table_count = 0
    for inZipFastqcFile in in_fastqc_ziped_files:
    	print "\t" + inZipFastqcFile
        row_table_count = row_table_count + 1
        file_name = ntpath.basename(inZipFastqcFile)
        file_name = file_name.split(".")[0]
        inSummary_file = in_SUPPLIMENTAL_QC_DIR + "/" + file_name + ".summary.txt"
        inFastqc_data_file = in_SUPPLIMENTAL_QC_DIR + "/" + file_name + ".fastqc_data.txt"

        summary_dict = {}
        with open(inSummary_file) as sum_file:
            for sum_file_line in sum_file:
                (val,key,file) = sum_file_line.split("\t")
                summary_dict[key]=val
                summary_dict["Sample"]=file

#            print ("\t" + keys[0] + "\t" + summary_dict[keys[0]].rstrip())
        for aKey in keys:
            title = aKey + " (" + summary_dict[aKey].rstrip() + ")"
            if summary_dict[aKey] == "PASS":
                HTML_FILE.write("<td bgcolor=\"#00FF00\" title=\"" + title + " \"></td>")
            elif summary_dict[aKey] == "WARN":
                HTML_FILE.write("<td bgcolor=\"#FFD700\" title=\"" + title + " \"></td>")
            elif summary_dict[aKey] == "FAIL":
                HTML_FILE.write("<td bgcolor=\"#FF0000\" title=\"" + title + " \"></td>")
            else:
                HTML_FILE.write("\n<td>" + str(row_table_count) + "</td>" + "<td bgcolor=\"#FFFFFF\" title=\"" + file_name + "\">" + "<a href=#" + file_name + ">" + file_name + "</a>"+ "</td>")
                
                for grep_search_term in basic_stats[3:]:
					grep_result = subprocess.Popen(["grep", "-m", "1" ,grep_search_term, inFastqc_data_file], stdout=subprocess.PIPE)
					result = grep_result.stdout.read()
					result = result.rstrip()
					result = result.lstrip(grep_search_term)
					result = result.lstrip()
#					print (grep_search_term + ": " + result)
					HTML_FILE.write("<td>" + result + "</td>")






            if aKey == "Sample":
                row_count = row_count + 1
                if row_count == 1:
                    first_row_name = summary_dict[aKey].rstrip()
                    R_SCRIPT_FILE.write( first_row_name + "<-c(\"" + summary_dict[aKey].rstrip() + "\"")
                else:
                    R_SCRIPT_FILE.write("df$" + summary_dict[aKey].rstrip() + "<-c(\"" + summary_dict[aKey].rstrip() + "\"")
            else:
                # R_SCRIPT_FILE.write(",\"" +  summary_dict[aKey] + "\"")
                if summary_dict[aKey] == "PASS":
                    R_SCRIPT_FILE.write(",1L")
                elif summary_dict[aKey] == "WARN":
                    R_SCRIPT_FILE.write(",0L")
                elif summary_dict[aKey] == "FAIL":
                    R_SCRIPT_FILE.write(",-1L")

			
        R_SCRIPT_FILE.write(")\n")

        if row_count == 1:
            R_SCRIPT_FILE.write("\ndf = data.frame(" + first_row_name + ")\n")


        HTML_FILE.write("</th></tr>")
        

    HTML_FILE.write("</div>" + "\n")
    HTML_FILE.write("</table>")
    HTML_FILE.write("<br><br>Green  - Test Metric pass" + "\n")
    HTML_FILE.write("<br>Yellow - Test Metric warning" + "\n")
    HTML_FILE.write("<br>Red    - Test Metric fail" + "\n")
    HTML_FILE.write("<hr>" + "\n")
	
    row_counter = 0
    for inZipFastqcFile in in_fastqc_ziped_files:
        row_counter = row_counter + 1
        file_name = ntpath.basename(inZipFastqcFile)
        file_name = file_name.split(".")[0]
        
        inSummary_file = in_SUPPLIMENTAL_QC_DIR + "/" + file_name + ".summary.txt"

        summary_dict = {}

        with open(inSummary_file) as sum_file:
            for sum_file_line in sum_file:
                (val,key,file) = sum_file_line.split("\t")
                summary_dict[key]=val
                summary_dict["Sample"]=file
        
        HTML_FILE.write("\n\n" + "<div class=\"outer\">" + "\n")
        HTML_FILE.write("<h3><a name=\"" + file_name + "\">" + str(row_counter) + ": " + file_name + " </a><a href=\"#STATS_TABLE\">(Return to Basic Statistics and QC Test Status Table)</a></h3>" + "\n")
        
        
        for img in extract_files[2:]:
        	img = img.rstrip() 
        	img_file = SUPPLIMENTAL_QC_DIR_NAME + "/" + file_name + "." + img
        	if summary_dict[file_key_dict[img]] == "PASS":
        		HTML_FILE.write("<a href=\"" + img_file +"\" target=\"_blank\" >" + "<img class=\"pass\" src=\"" + img_file + "\" /></a>" + "\n")
        	elif summary_dict[file_key_dict[img]] == "WARN":
        		HTML_FILE.write("<a href=\"" + img_file +"\" target=\"_blank\" >" + "<img class=\"warn\" src=\"" + img_file + "\" /></a>" + "\n")
        	elif summary_dict[file_key_dict[img]] == "FAIL":
        		HTML_FILE.write("<a href=\"" + img_file +"\" target=\"_blank\" >" + "<img class=\"fail\" src=\"" + img_file + "\" /></a>" + "\n")
        	else:
        		print "error: " + img
    HTML_FILE.write("<span class=\"fix\"></span></div>" + "\n")

    # Work on table in R
    R_SCRIPT_FILE.write("df <- t(df)\n")
    R_SCRIPT_FILE.write("colnames(df) <- keys\n")
    R_SCRIPT_FILE.write("df <- df[,2:13]\n")
    R_SCRIPT_FILE.write("colors <- c(\"red\", \"yellow\", \"green\")\n")
    R_SCRIPT_FILE.write("sample <- rownames(df)\n")
    R_SCRIPT_FILE.write("df_melt <- melt(cbind(df))\n")
    R_SCRIPT_FILE.write("ggplot(df_melt, aes(x = df_melt$Var2, y = df_melt$Var1, fill = factor(value))) +\n")
    R_SCRIPT_FILE.write("    theme(panel.grid.minor = element_line(colour=\"white\", size=0.5)) +\n")
    R_SCRIPT_FILE.write("    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +\n")
    R_SCRIPT_FILE.write("    theme(axis.text.y = element_text(angle = 45, hjust = 1)) +\n")
    R_SCRIPT_FILE.write("    labs(title = \"Summary of NGS Quality\", x = \"Sequence Metrics\", y = \"Fastq Files\") +\n")
    R_SCRIPT_FILE.write("    geom_tile(colour = \"black\") +\n")

    R_SCRIPT_FILE.write("    scale_fill_manual(values=colors,name=\"Tests\", labels=c(\"Failing\", \"Warning\", \"Passing\"))\n")
    
    subprocess.call(["rm", inSummary_file])
    HTML_FILE.close()
    R_SCRIPT_FILE.close()

########### Main ###################################################################################
## Credit modified from https://docs.python.org/2/library/getopt.html

def main():
    try:
        # These should be parameters

        SEARCH_SUFFIX_FILE_NAME = sys.argv[1]
        SEARCH_DIR = sys.argv[2]
        AGGREGATE_REPORTS_DIR = sys.argv[3]

        SUPPLIMENTAL_QC_DIR = AGGREGATE_REPORTS_DIR + "/" + SUPPLIMENTAL_QC_DIR_NAME
        
        subprocess.call(["mkdir", "-p", SUPPLIMENTAL_QC_DIR])
                
        html_file = AGGREGATE_REPORTS_DIR + "/aggregate_fastqc_report.html"
        R_script_file = AGGREGATE_REPORTS_DIR + "/aggregate_fastqc_report.R"

        R_SCRIPT_FILE = open(R_script_file, 'w')
        R_SCRIPT_FILE.write("library(ggplot2)\n")
        R_SCRIPT_FILE.write("library(reshape2)\n")

        R_SCRIPT_FILE.close()

        write_html_header(html_file)
        out_fastqc_ziped_files = extract_fastqc_zip_files(SEARCH_SUFFIX_FILE_NAME, SEARCH_DIR, SUPPLIMENTAL_QC_DIR)
        process_summary_file(html_file, R_script_file, out_fastqc_ziped_files, SUPPLIMENTAL_QC_DIR)

        write_html_closing(html_file)
        R_SCRIPT_FILE.close()

    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)


if __name__ == "__main__":
    main()
