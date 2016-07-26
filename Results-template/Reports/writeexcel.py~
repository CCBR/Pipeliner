import os,sys


samples=["FM403N_10000"]
H={}

sys.path.insert(0, "/home/dwheeler/bin")

import xlsxwriter

cols={"Sample ID":1, "Yield (Mbases)":2, "Total Reads (PF)":3, "% of >= Q30 bases":4, "Total trimmed reads":5, "Mapped Reads":6, "% Aligned":7, "%Mapped OnTarget":8, "%Unique Reads":9,"Mean Coverage on Target":10, "%Coverage >=30x":11}

workbook = xlsxwriter.Workbook('bamstats.xlsx')
worksheet = workbook.add_worksheet()


for i in cols.keys():
    worksheet.write(0,cols[i],i)


row=0
for s in samples:
    row+=1
    H["%Mapped OnTarget"]=os.popen("./cal.on.target.pl {0}.recal.bam --target_bed /data/CCBR/local/lib/SS_exome.bed".format(s)).read().split()[-1]
    H["Sample ID"]=s
    F=open("{0}{1}".format(s,".flagstats"),"r").read()

    for line in F.split("\n"):
        if line.find("total")>=0:
            H["Total Reads (PF)"]=line.split()[0]

        if line.find("mapped (")>=0:
            H["Mapped Reads"]=line.split()[0]

        print(line)

    for k in H.keys():
        worksheet.write(row,cols[k],H[k])




workbook.close()




