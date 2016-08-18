#!/usr/bin/env python3
import os
import re

### library of html sections and elements

header="""
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html>
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <title>CCBR Report</title>
    <link href="tabs.css" rel="stylesheet" type="text/css" />
</head>
<body>
"""
footer="""
<script src="activatables.js" type="text/javascript"></script>
<script type="text/javascript">
activatables('page', ['page-1', 'page-2', 'page-3']);
</script>
</body>
</html>
"""
div="""
<div class="content" id="page-1">
    <h2>Page 1</h2>
    <p>Text...</p>
</div>
"""



divs=""

li="""
    <li><a href="#page-1"><span>Page 1</span></a></li>
"""
### code

tabs=[]
page=""



def maketabset(D):
    global page
    global tabs
    divs=""
    toc="""<ol id="toc">"""
    for k in sorted(D.keys()):
        text="<table>"
        thediv=re.sub("page-1",k,div)
        thediv=re.sub("Page 1",k,thediv)
        processedDict=0
        P=""
        for k2 in sorted(D[k].keys()):
            try:
                k3=D[k][k2].keys()

                if processedDict==1:
                    pass
                else:
                    text=text+maketabset(D[k])
                    thetab={k:sorted(D[k].keys())}
                processedDict=1                        
            except:
                if re.search("include",k2):
                    for path in D[k][k2]:
                         if re.search("html",path):
                             text=text+"<iframe src='{}' width='100%' style='height: 100vh;'></iframe>".format(path)
                             thetab=k
                         if re.search("png",path):
                             text=text+"<image src='{}'></iframe>".format(path)
                             thetab=k
                         if re.search("pdf",path):
                             text=text+"<iframe src='{}'  width='100%' style='height: 100vh;'></iframe>".format(path)
                             thetab=k
                         if re.search("json",path):
                             text=text+"<iframe src='{}'  width='100%' style='height: 100vh;'></iframe>".format(path)
                             thetab=k
 
                else:
                    text=text+"<tr><td>{}</td><td>{}</td></tr>".format(k2,path)
                    thetab=k
        text=text+"</table>"
        thediv=re.sub("Text...",text,thediv)        
        theli=re.sub("page-1",k,li)
        theli=re.sub("Page 1",k,theli)
        toc=toc+theli
        divs=divs+thediv
        tabs.append(thetab)

    k="Contents"
    thetab=k
    toc2=re.sub("toc","toc2",toc)
    thediv="<div class='content' id='Contents'><h2>Contents</h2><p><table>"+toc2+"</ol></table></p>"+"</div>"
    theli=re.sub("page-1",k,li)
    theli=re.sub("Page 1",k,theli)
    toc=toc+theli
    divs=divs+thediv
    tabs.append(thetab)
    
    return toc+divs




J=eval(open("report.json","r").read())
P=maketabset(J)
page=page+P
footer=re.sub("\['page-1', 'page-2', 'page-3'\]",str(tabs),footer)
page=header+page+footer
print(page)

