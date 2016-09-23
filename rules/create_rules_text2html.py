#!/usr/bin/python
# Developer: Carl McIntosh

## Usage example
# python create_rules_text2html.py \
# --rule_file test/annovar.rl \
# --rule_name annovar \
# --out_file test/rule_annovar.html

# python create_rules_text2html.py \
# -i test/annovar.rl \
# -l annovar \
# -o test/rule_annovar.html



import sys, getopt

def main(argv):
   rulefile = ''
   rulename = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:l:o:",["rule_file=", "rule_name=","out_file="])
   except getopt.GetoptError:
      print 'Error: python create_rules_text2html.py -i <rule_file> -l <rule_name> -o <out_file>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'python create_rules_text2html.py -i <rule_file> -l <rule_name> -o <out_file>'
         sys.exit()
      elif opt in ("-i", "--rule_file"):
         rulefile = arg
      elif opt in ("-l", "--rule_name"):
         rulename = arg
      elif opt in ("-o", "--out_file"):
         outputfile = arg

    ## Write subject/database fasta
   OUT_HTML_FILE = open(outputfile, 'w')


   OUT_HTML_FILE.write("<!DOCTYPE html>\n")
   OUT_HTML_FILE.write("<html lang=\"en\">\n")
   OUT_HTML_FILE.write("<head>\n")
   OUT_HTML_FILE.write("     <meta charset=\"utf-8\" />\n")
   OUT_HTML_FILE.write("     <title>Rule Documentation: " + rulename + "</title>\n")
   OUT_HTML_FILE.write("     <link href=\"../assets/bootstrap.css\" rel=\"stylesheet\">\n")
   OUT_HTML_FILE.write("<style type=\"text/css\">\n")
   OUT_HTML_FILE.write("     pre {\n")
   OUT_HTML_FILE.write("          display: block;\n")
   OUT_HTML_FILE.write("          font-family: monospace;\n")
   OUT_HTML_FILE.write("          white-space: pre;\n")
   OUT_HTML_FILE.write("          background-color:white;\n")
   OUT_HTML_FILE.write("          border: none;\n")
   OUT_HTML_FILE.write("          margin: 1em 0; }\n")


   OUT_HTML_FILE.write("     body {\n")
   OUT_HTML_FILE.write("          padding-left: 32px;\n")
   OUT_HTML_FILE.write("          padding-right: 32px;\n")
   OUT_HTML_FILE.write("          min-width: 800px;\n")

   OUT_HTML_FILE.write("</style>\n")
   OUT_HTML_FILE.write("</head>\n")
   OUT_HTML_FILE.write("<body>\n")
   OUT_HTML_FILE.write("     <!-- Code Parameter -->\n")
   OUT_HTML_FILE.write("     <div style=\"border-radius: 10px; background-color:Gainsboro; padding: 10px;\">\n")
   OUT_HTML_FILE.write("     <p><b>" + rulename + ".rl</b> - Description\n")
   OUT_HTML_FILE.write("     <br>Version: x.x</p>\n")
   OUT_HTML_FILE.write("     <div style=\"border-radius: 10px; background-color:white; ; padding: 10px;\">\n")
   OUT_HTML_FILE.write("     <font color=\"red\">\n")
   OUT_HTML_FILE.write("     <font color=\"gray\"># Comment</font>\n")
   OUT_HTML_FILE.write("     <pre><font color=\"red\">\n")

   with open (rulefile, "r") as my_rule_file:
      rule_text=my_rule_file.read()



   OUT_HTML_FILE.write(rule_text + "\n")
   OUT_HTML_FILE.write("</font></pre>\n")
   OUT_HTML_FILE.write("</div>\n")
   OUT_HTML_FILE.write("</div>\n")
   OUT_HTML_FILE.write("</body>\n")
   OUT_HTML_FILE.write("</html>\n")


if __name__ == "__main__":
   main(sys.argv[1:])
