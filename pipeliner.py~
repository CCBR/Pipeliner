#!/usr/bin/env python3
import sys,os,math,time
import json,re
import contextlib
import webbrowser
import time,threading

# from subprocess import Popen, PIPE, STDOUT,check_output
from subprocess import *

from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.filedialog import askdirectory
from tkinter.messagebox import showerror
from tkinter.ttk import *

try:
    # Python2
    from Tkinter import *
    import Tkinter
    import tkMessageBox
    import tkFileDialog    

except ImportError:
    # Python3
    from tkinter import *
    import tkinter
    import tkinter.messagebox
    import tkinter.filedialog

#################
#Globals
#################

UnitsBak=[]
RGbak={}
cfile="pipeliner.json"
def load_config(cfile=cfile):
    global Cfg
    Q=open(cfile,"r")
    Cfg=eval(Q.read())
    Q.close()

load_config()    

def set_colors(*args):
    global baseColor
    global entryFgColor
    global entryBgColor
    global textDarkColor
    global textLightColor
    global widgetFgColor
    global widgetBgColor
    global projectFgColor
    global projectBgColor
    global statusFgColor
    global statusBgColor
    global dryrunFgColor
    global dryrunBgColor
    global logFgColor
    global logBgColor
    global errorsFgColor
    global errorsBgColor
    global menuFgColor
    global menuBgColor
    global commentFgColor
    global commentBgColor
    
    baseColor=Cfg['colors']['base'][0]
    entryFgColor=Cfg['colors']['entry'][0]
    entryBgColor=Cfg['colors']['entry'][1]
    textDarkColor=Cfg['colors']['text'][0]
    textLightColor=Cfg['colors']['text'][1]
    widgetFgColor=Cfg['colors']['widget'][0]
    widgetBgColor=Cfg['colors']['widget'][1]
    projectFgColor=Cfg['colors']['project'][0]
    projectBgColor=Cfg['colors']['project'][1]
    statusFgColor=Cfg['colors']['status'][0]
    statusBgColor=Cfg['colors']['status'][1]
    
    dryrunFgColor=Cfg['colors']['dryrun'][0]
    dryrunBgColor=Cfg['colors']['dryrun'][1]
    logFgColor=Cfg['colors']['log'][0]
    logBgColor=Cfg['colors']['log'][1]
    
    errorsFgColor=Cfg['colors']['errors'][0]
    errorsBgColor=Cfg['colors']['errors'][1]
    
    menuFgColor=Cfg['colors']['menu'][0]
    menuBgColor=Cfg['colors']['menu'][1]

    commentFgColor=Cfg['colors']['comments'][0]
    commentBgColor=Cfg['colors']['comments'][1]

    
    #baseColor="red4"
    #widgetBgColor="goldenrod"


    
set_colors()
customRules=[]
pipeline=["initialqc","bam2recal","wgslow","exomeseq-somatic","exomeseq-germline","exomeseq-germline-recal","exomeseq-germline-partial"]
auto_check=0

whereiam=os.popen("cd ../ && pwd").read().strip()
runmode=1
defaultwork=whereiam+"/Results"

#tkinter.messagebox.showinfo("path",whereiam)

rules=os.popen("cd Rules && ls *.rl|grep -v all-").read().replace(".rl","").split("\n")

b=rules.pop()

#rules=['samtools.sort','picard.headers','picard-markdups']
PL="initialqc"
smver=os.popen("snakemake --version").read().strip()
PJ=open("project.json","r")
projectjson=PJ.read()
PJ.close()
errorreport=""
batchsize="20"
parameters=[]
smparams=[]

newlines = ['\n', '\r\n', '\r']

####################
#Functions
###################

def writeheader(*args):
    if ftype.get()=="rg.tab":
        comments.delete(1.0,END)
        comments.insert(INSERT,"rgid\trgsm\trglb\trgpl\trgpu\trgcn\n")
    if ftype.get()=="pairs":
        comments.delete(1.0,END)
        comments.insert(INSERT,"Sample1\tSample2\n")
    return

def writepaste():
    try:
        fname=workpath.get()+"/"+ftype.get()    
        F=open(fname,"w")
        F.write(comments.get('1.0',END))
        F.close()
        tkinter.messagebox.showinfo("Success","Wrote file "+fname)
    except:
        tkinter.messagebox.showinfo("Error","Did not write file "+fname+"\nIs working directory set?")
    makejson()
    return

def load_configuration():
    fname = askopenfilename(filetypes=(("json files", "*.json"),
                                       ("All files", "*.*") ))
    if fname:
        try:
            load_config(cfile=fname)
        except:                     # <- naked except is a bad idea
            showerror("Open Source File", "Failed to read config file\n'%s'" % fname)
            return
        
    set_colors()
    top.configure(bg=baseColor)
#    i=top.winfo_children()
#    showerror("", i)

    _list = top.winfo_children()

    for item in _list :
        if item.winfo_children() :
            _list.extend(item.winfo_children())
            
    for i in filter(lambda w:isinstance(w,LabelFrame), _list):
        i.config(fg=textLightColor,bg=baseColor)
    for i in filter(lambda w:isinstance(w,Frame), _list):
        i.config(fg=textLightColor,bg=baseColor)
    for i in filter(lambda w:isinstance(w,Menu), _list):
        i.config(fg=menuFgColor,bg=menuBgColor)
    for i in filter(lambda w:isinstance(w,Radiobutton), _list):
        i.config(fg=widgetFgColor,bg=widgetBgColor)

    for i in filter(lambda w:isinstance(w,OptionMenu), _list):
        i.config(fg=widgetFgColor,bg=widgetBgColor)
#        i.config["menu"](fg=widgetFgColor,bg=widgetBgColor)        
    for i in filter(lambda w:isinstance(w,Label), _list):
        i.config(fg=textLightColor,bg=textDarkColor)
    for i in filter(lambda w:isinstance(w,Checkbutton), _list):
        i.config(fg=textLightColor,bg=textDarkColor)

    return


def load_project():
    global customRules
    fname = askopenfilename(filetypes=(("json files", "*.json"),
                                       ("All files", "*.*") ))
    if fname:
        try:
#            print("here is %s" % fname)
            F=open(fname,"r")
            PD=F.read()
            PD=eval(PD)
            J=json.dumps(PD, sort_keys = True, indent = 4, ensure_ascii=TRUE)
            jsonconf.delete("1.0", END)    

            jsonconf.insert(INSERT, J)
            customRules=PD['project']['custom']
            workpath.set(PD['project']['workpath'])
            datapath.set(PD['project']['datapath'])
            annotation.set(PD['project']['annotation'])
            binset.set(PD['project']['binset'])
            efiletype.set(PD['project']['efiletype'])
            filetype.set(PD['project']['filetype'])
            Pipeline.set(PD['project']['pipeline'])
            epoc.set(PD['project']['poc'])
            eanalyst.set(PD['project']['analyst'])
            euser.set(PD['project']['username'])
            eplatform.set(PD['project']['platform'])
            epi.set(PD['project']['pi'])
            eprojectid.set(PD['project']['id'])
            eorganism.set(PD['project']['organism'])
            technique.set(PD['project']['technique'])
            description.insert(INSERT, PD['project']['description'])
            
            
#            tkinter.messagebox.showinfo(fname,customRules)            
            F.close()
        except:                     # <- naked except is a bad idea
            showerror("Open Source File", "Failed to read file\n'%s'" % fname)
    makejson()
    #MkaS=os.popen("./makeasnake.py 2>&1 | tee -a "+workpath.get()+"/Reports/makeasnake.log").read()
    return

def set_working_directory():
    fname = tkinter.filedialog.askdirectory(initialdir=whereiam,title="Select Working Directory")    
    workpath.set(fname)                                    

def set_data_directory():
    fname = tkinter.filedialog.askdirectory(initialdir=whereiam,title="Select Data Directory")
    datapath.set(fname)                                    

def save_project():
    try:
        fname = asksaveasfilename(filetypes=(("json files", "*.json"),("All files","*.*")),initialdir=workpath.get())
    except:
        tkinter.messagebox.showinfo("Error","Project not Saved\n Is the working directory set?")
        return
    
    if fname:
        P=eval(jsonconf.get("1.0",END))
 #       tkinter.messagebox.showinfo("Json",str(P))
        try:
            with open(fname, 'w') as F:
                json.dump(P, F, sort_keys = True, indent = 4, ensure_ascii=False)
            F.close()
            tkinter.messagebox.showinfo("Project Write","Project written to %s."% fname)
        except Exception as e:
            tkinter.messagebox.showerror("Error","Name=%s  Error= %s"%(fname,e))

    return

def progress():
    o=os.popen("cd {0} && snakemake --dryrun --rerun-incomplete > {0}/Reports/checkpoint".format(workpath.get()))
    F=open("{0}/Reports/checkpoint".format(workpath.get()),"r").read()
    rules2={}
    rules=re.findall(r'rule .+:',F)
    for i in rules:
        i=re.sub("rule ","",i)
        i=re.sub(":","",i)
        rules2[i]=0
    
        
    F=open("{0}/Reports/{1}.dot".format(workpath.get(),Pipeline.get()),"r").read()
    for i in rules2.keys():

        F=re.sub(r'('+i+')(\".+?)\".+?\"',r'\1_pending\2"0.0 0.0 0.0"',F)
#        F=re.sub(i,"",F)        
        

    G=open("{0}/Reports/{1}-{2}.dot".format(workpath.get(),Pipeline.get(),"progress"),"w")
    G.write(F)
    G.close()

    o=os.popen("cd {0}/Reports && dot -Tpng -o {0}/Reports/{1}-progress.png {0}/Reports/{1}-progress.dot;convert {0}/Reports/{1}-progress.png {0}/Reports/{1}-progress.gif".format(workpath.get(),Pipeline.get()))

#    tkinter.messagebox.showerror("o",o)

    PL=Pipeline.get()
    gf=Toplevel()

    gf.title("CCBR Pipeliner: {0} Progress Graph".format(PL))
    cgf = Canvas(gf,bg="white")
#    gff=Frame(cgf,width=300,height=300)
    xscrollbar = Scrollbar(gf, orient=HORIZONTAL)
    xscrollbar.pack(side = BOTTOM, fill=X )
    xscrollbar.config(command=cgf.xview)
    
    yscrollbar = Scrollbar(gf,orient=VERTICAL)
    yscrollbar.pack(side = RIGHT, fill=Y )
    yscrollbar.config(command=cgf.yview)
    
    cgf.config(xscrollcommand=xscrollbar.set, yscrollcommand=yscrollbar.set)
    cgf.config(width=600,height=600)
    cgf.pack(expand=1,fill=BOTH,side=RIGHT)
    cgf.config(scrollregion=(0,0,1000,5000))
    try:
        time.sleep(5)
        img = PhotoImage(file="{0}/Reports/{1}-progress.gif".format(workpath.get(),PL))
    except:
        time.sleep(5)
        img = PhotoImage(file="{0}/Reports/{1}-progress.gif".format(workpath.get(),PL))
    cgf.create_image(0,0,image=img, anchor="nw")
    cgf.image=img
    

def unbuffered(proc, stream='stdout'):
    stream = getattr(proc, stream)
    with contextlib.closing(stream):
        while True:
            out = []
            last = stream.read(1)
            # Don't loop forever
            if last == '' and proc.poll() is not None:
                break
            while last not in newlines:
                # Don't loop forever
                if last == '' and proc.poll() is not None:
                    break
                out.append(last)
                last = stream.read(1)
            out = ''.join(out)
            yield out


def makejson(*args):
    global PD
    global UnitsBak
    global RGbak
    D=dict()
    try:
        F=open(workpath.get()+"/samples","r")
        f=F.read().split("\n")
        F.close()
        for line in f:
            L=line.split()
            a=L.pop(0)
            D[a]=L
            samples=D
    except:
        samples={"na":"na"}
    
    
    D=dict()
    try:
        F=open(workpath.get()+"/pairs","r")
        f=F.read().split()
        F.close()
        for i in range(0,len(f),2):
#            a=f[i].split(".")[0]
#            b=f[i+1].split(".")[0]
            a=f[i]
            b=f[i+1]

            D[a+"+"+b]=[a,b]

        pairs=D
    except:
        pairs={"na":"na"}   

   
    D=dict() 
    FT=filetype.get()
#    p = Popen("ls "+workpath.get()+"/*."+FT, shell=True, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, close_fds=True)
    p = Popen("find "+workpath.get()+" -maxdepth 1 -type l -printf '%f\n' ", shell=True, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, close_fds=True)
    a = p.stdout.read().decode(encoding='UTF-8').split("\n")

    RG=dict()   
    b=a.pop()
#    tkinter.messagebox.showerror("",a)
#    if freezeunits.get()=="no":
    for i in a:

        key=re.sub(".realign","",i.split("/")[-1])
        key=re.sub(".bai","",key)
        key=re.sub(".bam","",key)
        key=re.sub(".sam","",key)        
        key=re.sub(".recal","",key)
        key=re.sub(".dedup","",key)
        key=re.sub(".sorted","",key)
        key=re.sub(".fin","",key)
        key=re.sub(".R[12]","",key)
        key=re.sub("_R[12]","",key)
        key=re.sub(".fastq","",key)
        key=re.sub(".gz","",key)                                
#        key=re.sub("[\._](R[12]\.)*"+FT+"$","",i.split("/")[-1])        
#        key=re.sub(".R[12]."+FT+"$","",i.split("/")[-1])
#        key=re.sub("([._]R[12][._])*([fin|sorted|dedup|recal|realign])*\.{0}$".format(FT),"",i.split("/")[-1])
        D[key]=key
        RG[key]={'rgsm':key,'rglb':'na','rgpu':'na','rgpl':'ILLUMINA','rgcn':'na'}
    units=D
    UnitsBak=D

    try:
        F=open(workpath.get()+"/rg.tab","r")
        f=F.read().splitlines()
        F.close()
        for theLine in f:
            if not re.match("^ID",theLine):
                (rgid,rgsm,rglb,rgpl,rgpu,rgcn)=theLine.split("\t")
                RG[rgid]={'rgsm':rgsm,'rglb':rglb,'rgpu':rgpu,'rgpl':rgpl,'rgcn':rgcn}
    except:
        pass
    RGbak=RG
#     else:
#         units=UnitsBak
#         RG=RGbak
# 
    PD=dict()

    smparams=[]

    for i in range(len(parameters)):

        if cp[i].var.get()=="1":
            smparams.append(parameters[i])

    
    PD={'project':{'pfamily':pfamily.get(),'units':units,'samples':samples,'pairs':pairs,
                   'id':eprojectid.get(),'pi':epi.get(),'organism':eorganism.get(),
                   'analyst':eanalyst.get(),'poc':epoc.get(),'pipeline':Pipeline.get(),'version':"1.0",
                   'annotation':annotation.get(),'datapath':datapath.get(),'filetype':filetype.get(), 'binset':binset.get(),'username':euser.get(),'platform':eplatform.get(),'custom':customRules,'efiletype':efiletype.get(),'workpath':workpath.get(),'batchsize':batchsize,"smparams":smparams,"rgid":RG,"cluster":"cluster_medium.json","description":description.get('1.0',END),"technique":technique.get()}}

    
    J=json.dumps(PD, sort_keys = True, indent = 4, ensure_ascii=TRUE)
    jsonconf.delete("1.0", END)    

    jsonconf.insert(INSERT, J)
    saveproject(jsonconf.get("1.0",END))
    #MkaS=os.popen("./makeasnake.py 2>&1 | tee -a "+workpath.get()+"/Reports/makeasnake.log").read()
#    tkinter.messagebox.showinfo("Project Json Build","Project Json Built.")    

def initialize():  
    global initLock  

    if initLock.get()=="unlocked":
        pass
    else:

        p=os.popen("if [ ! -d {0} ]; then mkdir {0};fi".format(workpath.get()))

        p=os.popen("rm -rf {0}/*stats rm -rf {0}/QC;rm -rf {0}/Reports;rm -rf {0}/*dedup_stats; rm {0}/*; rm -rf {0}/.*; mkdir {0}/QC;touch {0}/pairs;touch {0}/samples;cp -rf ".format(workpath.get())+whereiam+"/Pipeliner/Results-template/* {0}".format(workpath.get()))

    
def initialize_results():
    global initLock
#    tkinter.messagebox.showinfo("initLock","initLock={0}".format(initLock))

    if initLock.get()=="unlocked":
        tkinter.messagebox.showinfo("Locked","Initialize button is locked. Uncheck to unlock.")
        pass
    else:
        result=tkinter.messagebox.askquestion("Initialize Directory Warning", "Initialize working directory %s?  ALL FILES IN THIS DIRECTORY WILL BE DELETED!"%workpath.get(), icon='warning')
        if result=='yes':
            p=os.popen("ls {0}".format(workpath.get())).read()
            result=tkinter.messagebox.askquestion("Initialize Directory Warning", "THESE FILES WILL BE DELETED! Continue? \n {0}".format(p), icon='warning')
            if result=='yes':
                initialize()
                tkinter.messagebox.showinfo("Initializing Directory","Directory Initialized")
            else:
                tkinter.messagebox.showinfo("Aborted Initialing Directory","Directory Not Initialized")
        
        else:
            tkinter.messagebox.showinfo("Aborted Initializing Directory","Directory Not Initialized")

    
def symlink(data):
    data=data+"/"
    data=re.sub("/+$","/",data)
    FT=filetype.get()
    try:
        #cmd="for f in `ls "+data+"*.fastq`;do ln -s $f;done"
#        cmd="for f in `ls "+data+"*."+FT+"`;do ln -s $f ../;done"
        cmd="for f in `ls {0}*[._]{1}`;do ln -s $f {2};done".format(data,FT,workpath.get())        
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        Out = p.stdout.read()
    except Exception as e:
        tkinter.messagebox.showinfo("Error",str(e))
    if str(Out).find("No such file")==-1:    
        tkinter.messagebox.showinfo(FT,"Symlinks Created")
    else:
        tkinter.messagebox.showinfo("Error","Symlinks Not Created")
    if re.sub("bam","",FT):
        FT2=re.sub("bam","",FT)
        try:
            cmd="for f in `ls "+data+"*."+FT2+"bai`;do ln -s $f %s;done"%workpath.get()
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            Out = p.stdout.read()
        except Exception as e:
            tkinter.messagebox.showinfo("Error",str(e))
        if str(Out).find("No such file")==-1:    
            tkinter.messagebox.showinfo(FT2+"bai","Symlinks Created")
        #else:
            #tkinter.messagebox.showinfo(FT2+"bai","Index Symlinks Not Created")

    p = os.popen("for f in `ls {0}/*fastq*`;do mv $f `echo $f|sed s/_fastq/.fastq/g`;done".format(workpath.get()))
    p = os.popen("for f in `ls {0}/*fastq*`;do mv $f `echo $f|sed s/_R1.fastq/.R1.fastq/g`;done".format(workpath.get()))
    p = os.popen("for f in `ls {0}/*fastq*`;do mv $f `echo $f|sed s/_R2.fastq/.R2.fastq/g`;done".format(workpath.get()))




    makejson()
    
def show(x,y,fg,bg,height,width):
    sm=Toplevel()
    sm.title(y)
    smyscrollbar = Scrollbar(sm, orient=VERTICAL)
    smyscrollbar.pack( side = RIGHT, fill=Y )
    smxscrollbar = Scrollbar(sm, orient=HORIZONTAL)
    smxscrollbar.pack( side = BOTTOM, fill=X )
    smout = Text(sm,width=width,height=height,bg=bg,fg=fg,yscrollcommand = smyscrollbar.set, xscrollcommand = smxscrollbar.set)
#    smout = Text(sm,width=100,height=100,bg="gray30",fg="yellow",yscrollcommand = smyscrollbar.set, xscrollcommand = smxscrollbar.set)
    smout.insert(INSERT, x)
    smout.pack(expand=1,fill=BOTH,side=RIGHT)
    smyscrollbar.config(command=smout.yview)
    smxscrollbar.config(command=smout.xview)
    
def makeunitswarm():
    S=""
    F=open("project.json","r")
    config=eval(F.read())
    for u in config['project']['units']:
        S=S+"module load python/3.4.0; samples=\""+ u + "\"; snakemake -s Snakefile\n"
    tkinter.messagebox.showinfo("Swarm",S)
    F=open("../"+PL+".swarm","w")
    F.write(S)
    F.close()    

def makepairswarm():
    S=""
    F=open("project.json","r")
    config=eval(F.read())
    for u in config['project']['pairs'].values():
        S=S+"module load python/3.4.0; samples=\""+ ' '.join(u) + "\"; snakemake -s Snakefile\n"
    tkinter.messagebox.showinfo("Swarm",S)
    F=open("../"+PL+".swarm","w")
    F.write(S)
    F.close()    

def getrootnames():
    D=dict()
    a=os.popen("ls ../*.fastq").read().split("\n")
    b=a.pop()
    for i in a:
        D[i.split(".")[0]]=i.split(".")[0]

    show(D,"Root Fastq Names to Units","black","white",50,50)
    
def twocol2pairs():
    D=dict()
    F=open("../pairs","r")
    f=F.read().split()
    F.close()
    for i in range(0,len(f),2):
        a=f[i].split(".")[0]
        b=f[i+1].split(".")[0]
        D[a+"+"+b]=[a,b]
    show(D,"Pairs","black","white",50,50)


def tab2samples():
    D=dict()
    F=open("../samples","r")
    f=F.read().split("\n")
    F.close()
    try:
        for line in f:
            L=line.split()
            a=L.pop(0)
            D[a]=L
    except:
        pass
    
    show(D,"Samples","black","white",50,50)
    
    

def about():
    info="""
    CCBR Pipeliner
    Version 1.0
    August 28, 2015
    """ 
    tkinter.messagebox.showinfo("CCBR Pipeliner\nVersion 1.0",info)

def hidew(self):
        """"""
        self.top.withdraw()

def showw(self):
        """"""
        self.top.update()
        self.top.deiconify()
 
def jsonform():
    return

def saveproject(proj):
    P=eval(proj)

    try:
        with open('project.json', 'w') as F:
            json.dump(P, F, sort_keys = True, indent = 4, ensure_ascii=False)
        F.close()
#        tkinter.messagebox.showinfo("Project Json Write","Project Json file written.")
    except:
        tkinter.messagebox.showerror("Error","Project Json file not written.")


def viewreport():
    PL=Pipeline.get()
    try:
        webbrowser.open(workpath.get()+'/Reports/'+PL+'.html')
    except Exception as e:
        tkinter.messagebox.showinfo("View Report",str(e))
        
def getworkflow():
    PL=Pipeline.get()
    gf=Toplevel()
    #MkaS=os.popen("./makeasnake.py 2>&1 | tee -a "+workpath.get()+"/Reports/makeasnake.log").read()

    gf.title("CCBR Pipeliner: "+ PL + " Workflow Graph")
    cgf = Canvas(gf,bg="white")
#    gff=Frame(cgf,width=300,height=300)
    xscrollbar = Scrollbar(gf, orient=HORIZONTAL)
    xscrollbar.pack(side = BOTTOM, fill=X )
    xscrollbar.config(command=cgf.xview)
    
    yscrollbar = Scrollbar(gf,orient=VERTICAL)
    yscrollbar.pack(side = RIGHT, fill=Y )
    yscrollbar.config(command=cgf.yview)
    
    cgf.config(xscrollcommand=xscrollbar.set, yscrollcommand=yscrollbar.set)
    cgf.config(width=600,height=600)
    cgf.pack(expand=1,fill=BOTH,side=RIGHT)
    cgf.config(scrollregion=(0,0,1000,5000))
    img = PhotoImage(file=workpath.get()+"/Reports/"+PL+".gif")
    cgf.create_image(0,0,image=img, anchor="nw")
    cgf.image=img

def cmd(smcommand):
    global img
    PL=Pipeline.get()
#    makejson()
    saveproject(jsonconf.get("1.0",END))
    MkaS=os.popen("./makeasnake.py 2>&1 | tee -a "+workpath.get()+"/Reports/makeasnake.log").read()
    Out=os.popen("cd "+workpath.get()+" && snakemake " + smcommand +" 2>&1").read()
    show(Out,"CCBR Pipeliner: "+ PL + " "+smcommand,dryrunFgColor,dryrunBgColor,200,100)

def swarm():
    return

def qsub():
    return
def custompipe():
    global cb
    custom=Toplevel(bg=baseColor)
    custom.title("Construct Custom Pipeline")    
    cb=[0 for x in range(len(rules))]
    y=0
    for i in range(len(rules)):
        if i%8==0:
            y=y+1
            x=0
        v=StringVar(value="0")

        cb[i] = Checkbutton(custom, text=rules[i],variable=v,bg=baseColor,fg=textLightColor,offvalue="0",onvalue="1")        
        cb[i].var=v
        cb[i].var.trace("w",setrules)
        cb[i].grid(row=x,column=y,sticky=W)        
        x=x+1
    custom.update()
    custom.deiconify()

def setrules(*args):
    global customRules
    customRules=[]
    for i in range(len(rules)):
        if cb[i].var.get()=="1":
#                tkinter.messagebox.showinfo(str(i),cb[i].var.get())
            customRules.append(rules[i])
    makejson()
    
def updateParams(*args):
    global batchsize
    batchsize=BS.get()    
    makejson()
    
def setparameters():
    global batchsize
    global parameters
    global cp
    #batch size for GVCFs
    custom=Toplevel(bg=baseColor)
    custom.title("Additional Parameters")    
    ppar = LabelFrame(custom,text="Additional Pipeline Parameters")
    ppar.grid(row=0,column=0,sticky=W,padx=10,pady=10)
    ppar.configure(bg=baseColor,fg=textLightColor)

    L=Label(ppar,text="Batch size for Combine GVCFs",anchor="ne",bg=widgetBgColor,fg=widgetFgColor)
    E = Entry(ppar, bd =2, width=15,bg=entryBgColor,fg=entryFgColor,textvariable=BS)
    L.grid(row=1,column=0,sticky=W)
    E.grid(row=2,column=0,sticky=W)
    BS.trace('w', updateParams)
    
    par = LabelFrame(custom,text="Additional Snakemake Parameters")
    par.grid(row=3,column=0,sticky=W,padx=10,pady=10)
    par.configure(bg=baseColor,fg=textLightColor)

    parameters=["--debug","--notemp","--reason","--forceall","--keep-going"]
    plabels=["Run in debug mode","Ignore temp directives in rules","Show reasons for rules","Force execution of all rules","Keep going if non-essential rule fails"]
    
    cp=[0 for x in range(len(parameters))]
    y=3

    for i in range(len(parameters)):
        if i%1==0:
            y=y+1
            x=0
        v=StringVar(value="0")

        cp[i] = Checkbutton(par, text=plabels[i],variable=v,bg=baseColor,fg=textLightColor,offvalue="0",onvalue="1")        
        cp[i].var=v
        cp[i].var.trace("w",makejson)
        cp[i].grid(row=y,column=x,sticky=W)        
        x=x+1
    



    custom.update()
    custom.deiconify()


    
def seelog():
    sm=Toplevel()
    sm.title("Snakemake Log")

    smyscrollbar = Scrollbar(sm, orient=VERTICAL)
    smyscrollbar.pack( side = RIGHT, fill=Y )
    smxscrollbar = Scrollbar(sm, orient=HORIZONTAL)
    smxscrollbar.pack( side = BOTTOM, fill=X )
    smout = Text(sm,width=80,font=("nimbus mono","10"),height=30,bg=logBgColor, fg=logFgColor,yscrollcommand = smyscrollbar.set, xscrollcommand = smxscrollbar.set)
    smout.pack(expand=1,fill=BOTH,side=RIGHT)
    smyscrollbar.config(command=smout.yview)
    smxscrollbar.config(command=smout.xview)

    smout.tag_config('rule',background="orange",foreground="black")
    smout.tag_config('input:',background="cyan",foreground="black")
    smout.tag_config('output:',background="yellow",foreground="black")
    smout.tag_config('Error',background="red",foreground="black")
    smout.tag_config('java',background="purple",foreground="black")
    smout.tag_config('threads:',background="white",foreground="black")
    smout.tag_config('%) done',background="blue",foreground="white")        
    
    #smout.insert(INSERT, f)
        
#    sm.pack()
# force drawing of the window
    #smout.update_idletasks()
    #smout.insert(INSERT,"Starting Snakemake Run "+PL+"\n")

    def refreshlog():
        F=open(workpath.get()+"/Reports/snakemake.log","r")
        f=F.read()
        F.close()
        smout.delete("1.0", END)
        smout.insert(INSERT, f)

        for R in ["rule","input:","output:","Error","threads:","%) done"]:
            start="1.0"            
            while 1:
                pos=smout.search(R,start,stopindex=END)
                if not pos:
                    break
                pos2=re.sub("\.\d+","",pos)+".end"
                pos3=re.sub("\.\d+","",pos)+".0"                
                smout.tag_add(R,pos3,pos2)                
                start=pos+"+1c"
            
        smout.update_idletasks()
        
        smout.see(END)
        
    refreshlog()
    
    button = Button(sm, text="Update", fg=widgetFgColor,bg=widgetBgColor,command=refreshlog)
    button.pack( side = BOTTOM,padx=5)
    button.pack( side = TOP,padx=5)

    
def run():
    global snakemakeRun
    PL=Pipeline.get()
    MkaS=os.popen("./makeasnake.py "+PL+" 2>&1 | tee -a makeasnake.log").read()
    tkinter.messagebox.showinfo("Run by Simple Script","Starting Snakemake Run "+PL+"\n")
    snakemakeRun=Popen("../run.sh")
    seelog()

def runqsub():
    global snakemakeRun
    global runmode
    runmode=1
    PL=Pipeline.get()
    MkaS=os.popen("./makeasnake.py "+PL+" 2>&1 | tee -a "+workpath.get()+"/Reports/makeasnake.log").read()
    tkinter.messagebox.showinfo("Qsub","Starting Qsub Snakemake Run "+PL+"\n")
    #snakemakeRun=Popen("../qsub.sh")
    snakemakeRun=Popen(workpath.get()+"/submit.sh")    
    #seelog()    

def runslurm():
    global snakemakeRun
    global runmode
    runmode=2
    PL=Pipeline.get()
    MkaS=os.popen("./makeasnake.py "+PL+" 2>&1 | tee -a "+workpath.get()+"/Reports/makeasnake.log").read()
    tkinter.messagebox.showinfo("Slurm","Starting Slurm Snakemake Run "+PL+"\n")
    #snakemakeRun=Popen("../qsub.sh")
    snakemakeRun=Popen(workpath.get()+"/submit_slurm.sh")    
    #seelog()    


def snakemakeoptions():
    sm=Toplevel()
    sm.title("Snakemake Options")
    flags=['--forceall','--touch','--rerun-incomplete']

    c = Checkbutton(sm, text="Forceall", variable=FcAl,fg=textLightColor)
    c.pack()

    c = Checkbutton(sm, text="Rerun Incomplete", fg=textLightColor, variable=ReIn)
    c.pack()
    
def donothing():
    return


def stats2html():



    date = os.popen('date').read()

    V = os.popen('snakemake -v')
    v = V.read()
    
        
    F = open('stats.css', 'r')
    CSS = F.read()
    F.close()
    
    F = open('toc.js', 'r')
    TOC = F.read()
    F.close()
    
    F = open('sorttable.js', 'r')
    SRT = F.read()
    F.close()
    
    
    F = open(workpath.get()+'/run.json', 'r')
    P = eval(F.read())
    F.close()
    pid = P['project']['id']
    org = P['project']['organism']
    user = P['project']['analyst']
    pipeline = P['project']['pipeline']
    pipeline_ver = P['project']['version']
    workpath=P['project']['workpath']
    
    os.popen("cd ../ && snakemake --detailed-summary > "+workpath+"/Reports/"+pipeline+".summary")
    
    
    try:
        F = open(workpath.get()+"/Reports/"+pipeline+".stats", 'r')
        f=F.read()
        D=eval(f)
        F.close()
    except:
        print("No "+pipeline+".stats file found. Exiting.\n")
        quit()
    F.close()
    
    title = 'Pipeline '+pipeline+' Stats'
    page="<html><head><style>"+CSS+"</style>"
    page=page+"<title>"+title+"</title>"
    
    page=page+"<script type='text/javascript' src='sorttable.js'></script>"
    page=page+"<script type='text/javascript' src='toc.js'></script>"
    page=page+"<script type='text/javascript' src='pipeliner.js'></script>"
    
    
    page=page+"</head>"
    
    page=page+"<body onload=generateTOC(document.getElementById('toc'));>"
    
    
    nsamp = len(P['project']['units'])
    
    page=page+'<h1>Project: <i>{}</i> with {} samples</h1>'.format(pid, nsamp)
    page=page+'<h5>Date: <i>{}</i></h5>'.format(date)
    page=page+'<h5>Organism: <i>{}</i></h5>'.format(org)
    page=page+'<h5>Analyst: <i>{}</i></h5>'.format(user)
    page=page+'<h5>Pipeline: <i>{}</i></h5>'.format(pipeline)
    page=page+'<h5>Pipeline Version: <i>{}</i></h5>'.format(pipeline_ver)
    page=page+'<h5>Snakemake Version: <i>{}</i></h5>'.format(v)
    
    page=page+'<h2>Contents</h2>'
    page=page+"<div id='toc'></div>"
    
    
    page=page+"<div id='main'>"
    page=page+'<h2>Pipeline Flow Diagram</h2>'
    #page=page+"<img src='"+pipeline+".gif'>"
    
    page=page+"<object type='image/svg+xml' data="+workpath+"/Reports/"+pipeline+".svg>Your browser does not support SVG</object>"
    
    #page=page+"<iframe src="+pipeline+".svg>Your browser does not support SVG</iframe>"
    
    page=page+'<h3>Total Computing Time for Pipeline: {} hours'.format(D['total_runtime']/3600)
    page=page+'<h2>Execution Times for Rules</h2>'
    page=page+"<table cols=2 class='sortable'>"
    page=page+'<th>Rule</th><th>Average Duration</th><th>Maximum Duration</th><th>Minimum Duration</th>'
    for k in D['rules'].keys():
        page=page+'<ul><tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(k,
                D['rules'][k]['mean-runtime'], D['rules'][k]['max-runtime'
                ], D['rules'][k]['min-runtime'])
    
    page=page+"</table>"
    
    page=page+'<h2>Starts, Stops, and Processing Times for Files'
    
    page=page+"<table cols=4 class='sortable'>"
    
    page=page+'<tr>'
    
    page=page+'<th>File</th><th>Start</th><th>Stop</th><th>Duration</th>'
    
    page=page+'</tr>'
    
    for k in D['files'].keys():
        page=page+'<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>'.format(k,
                 D['files'][k]['start-time'], D['files'][k]['stop-time'],
                 D['files'][k]['duration'])
    
    page=page+"</table>"
    
    
    
    page=page+'Contents of <i>project.json</i>'
    page=page+"<table cols=2 class='sortable'>"
    
    page=page+'<tr>'
    
    page=page+'<th>Parameter</th><th>Value</th>'
    
    page=page+'</tr>'
    
    for k in P.keys():
        T = '<table>'
        try:
            for k2 in P[k].keys():
                T2 = '<table>'
                try:
                    for k3 in P[k][k2].keys():
                        T2 = T2 \
                            + '<tr><td>{}</td><td>{}</td></tr>'.format(k3,
                                P[k][k2][k3])
                    T2 = T2 + '</table>'
                except:
                    T2 = P[k][k2]
    
                T = T + '<tr><td>{}</td><td>{}</td></tr>'.format(k2, T2)
            T = T + '</table>'
        except:
            T = P[k]
    
        page=page+'<tr><td>{}</td><td>{}</td></tr>'.format(k, T)
    
    page=page+"</table>"
    
    
    page=page+"<form name='fastqcform'>"
    
    page=page+'<h2>FastQC Reports</h2>'
    
    
    page=page+"<h5> Select Sample ID</h5><select name='fastqcselect' onchange='fastqc()'>"
    
    first=sorted(P['project']['units'].keys())[0]
                 
    for unit in sorted(P['project']['units'].keys()):
        page=page+"<option value='%s'>%s</option>" % (unit,unit)
    page=page+"</select>"
    
    
    page=page+"<object type='text/html' name='fastqcr1html' data='../QC/%s.R1_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first
    
    page=page+"<object type='text/html' name='fastqcr2html' data='../QC/%s.R2_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first
    
    
    page=page+'<h2>FastQC Reports for Trimmed Reads</h2>'
    
    
    page=page+"<h5> Select Sample ID</h5><select name='fastqcselecttrimmed' onchange='fastqctrimmed()'>"
    
    
    for unit in sorted(P['project']['units'].keys()):
        page=page+"<option value='%s'>%s</option>" % (unit,unit)
    page=page+"</select>"
    
    
    page=page+"<object type='text/html' name='fastqcr1trimmedhtml' data='../QC/%s.R1.trimmed_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first
    
    page=page+"<object type='text/html' name='fastqcr2trimmedhtml' data='../QC/%s.R2.trimmed_fastqc.html' width='800' height='800' scrolling='yes'></object>" % first
    
    
    
    
    
    page=page+'<h2>QualiMap Reports</h2>'
    
    
    page=page+"<h5> Select Sample ID</h5><select name='qualimapselect' onchange='qualimap()'>"
    for unit in sorted(P['project']['units'].keys()):
        page=page+"<option value='%s'>%s</option>" % (unit,unit)
    page=page+"</select>"
    
    
    page=page+"<object type='text/html' name='qualimaphtml' data='../QC/%s.qualimapReport.pdf' width='800' height='400' scrolling='yes'></object>" % first
    
    
    page=page+"<h2>Snakemake Files Summary</h2>"
    
    
    page=page+"<object type='text/html' name='summaryhtml' data='../Reports/%s.summary' width='800' height='400' scrolling='yes'></object>" % pipeline
    
    page=page+"</form>"
    
    page=page+"</body></html>"
    
    G = open(workpath.get()+"/Results/"+pipeline+".html", 'w')
    
    G.write(str(page))
    G.close()
    


#######################
## Gui components
#######################


top = Tk()
top.title("CCBR Pipeliner   (using Snakemake version "+smver+")")
top.resizable(1,1)
top.configure(bg=baseColor)

####################
#The Menu
####################

menubar = Menu(top)
menubar.configure(bg=menuBgColor,fg=menuFgColor)

############ File ############
mainmenu = Menu(menubar, tearoff=0)
mainmenu.add_command(label="Load Project",command=load_project)
mainmenu.add_command(label="Save Project",command=save_project)
mainmenu.add_command(label="Choose Data Source",command=set_data_directory)
mainmenu.add_command(label="Choose Working Directory",command=set_working_directory)
mainmenu.add_command(label="Load Configuration",command=load_configuration)
mainmenu.add_command(label="Exit", command=top.quit)
menubar.add_cascade(label="File", menu=mainmenu)

####### View ###############
viewmenu = Menu(menubar, tearoff=0)
viewmenu.add_command(label="Final Report", command=viewreport)
viewmenu.add_command(label="Snakemake Log", command=seelog)
viewmenu.add_command(label="Workflow Graphic", command=getworkflow)
viewmenu.add_command(label="Progress Graphic", command=progress)
viewmenu.add_command(label="Summary", command=lambda: cmd("--summary"))
viewmenu.add_command(label="Detailed Summary", command=lambda: cmd("--detailed-summary"))
viewmenu.add_command(label="Rules", command=lambda: cmd("--list"))
viewmenu.add_command(label="Target Rules", command=lambda: cmd("--list-target-rules"))
menubar.add_cascade(label="View", menu=viewmenu)


########## Run ##############
runmenu = Menu(menubar, tearoff=0)
runmenu.add_command(label="Initialize Working Directory", command=initialize_results)
runmenu.add_command(label="Symlink Data", command=lambda:symlink(datapath.get()))
runmenu.add_command(label="Dry Run", command=lambda: cmd("--dryrun -s %s/Snakefile --rerun-incomplete -d %s" % (workpath.get(),workpath.get())))               
#runmenu.add_command(label="Run Qsub", command=runqsub)
runmenu.add_command(label="Run Slurm", command=runslurm)
menubar.add_cascade(label="Run", menu=runmenu)

########## Tools ##############
toolmenu = Menu(menubar, tearoff=0)
toolmenu.add_command(label="Make Custom Pipeline", command=lambda: custompipe())
toolmenu.add_command(label="Set Additional Parameters", command=lambda: setparameters())
toolmenu.add_command(label="Generate Final Report", command=lambda: stats2html())
#toolmenu.add_command(label="Make Unit Swarm File", command=makeunitswarm)
#toolmenu.add_command(label="Make Pairs Swarm File", command=makepairswarm)

menubar.add_cascade(label="Tools", menu=toolmenu)

########## Help #############
aboutmenu = Menu(menubar, tearoff=0)
#aboutmenu.add_command(label="Procedure Summary", command=about)
menubar.add_command(label="About",command=about)
top.config(menu=menubar)

########## Dev ##############
devmenu = Menu(menubar, tearoff=0)
devmenu.add_command(label="Build Project Json", command=lambda:makejson())
devmenu.add_command(label="Save Project Json", command=lambda:saveproject(jsonconf.get("1.0",END)))
devmenu.add_command(label="Two Column Names to Pairs", command=twocol2pairs)
devmenu.add_command(label="Table to Samples", command=tab2samples)
toolmenu.add_cascade(label="Advanced", menu=devmenu)

##################
#The Frames
#################


opts = LabelFrame(top,text="")
opts.pack(side=TOP,fill=X)
#opts.configure(bg=baseColor,fg="white")

#opts3 = LabelFrame(opts,text="")
#opts3.pack(side=TOP,fill=X )
#opts3.configure(bg=baseColor,fg=textLightColor)

#opts2 = LabelFrame(opts,text="")
#opts2.pack(side=TOP,fill=X )
#opts2.configure(bg=baseColor,fg=textLightColor)


#editframe = LabelFrame(top,text="Project Information",bg=baseColor,fg=textLightColor)
#editframe.pack( side = LEFT,fill=BOTH )

nbook = ttk.Notebook(top)
projframe = ttk.Frame(nbook)
#filesframe = ttk.Frame(nbook)
editframe = ttk.Frame(nbook)
opts2 = ttk.Frame(nbook)
pastewriteframe = ttk.Frame(nbook)
runframe=ttk.Frame(nbook)
#manualframe=ttk.Frame(nbook)
rnaseqframe=ttk.Frame(nbook)
exomeseqframe=ttk.Frame(nbook)
mirseqframe=ttk.Frame(nbook)
chipseqframe=ttk.Frame(nbook)
customframe=ttk.Frame(nbook)

nbook.add(editframe, text='Project Information')
nbook.add(opts2, text='Global Options')
nbook.add(projframe, text='Project Json')
nbook.add(runframe,text="Run Sequence")
#nbook.add(filesframe, text='Job Status')
nbook.add(pastewriteframe, text='Paste/Write Files')
#nbook.add(manualframe,text="Manual")
nbook.add(rnaseqframe,text="RNASeq Options")
nbook.add(exomeseqframe,text="ExomeSeq Options")
nbook.add(mirseqframe,text="MirSeq Options")
nbook.add(chipseqframe,text="ChIPSeq Options")
nbook.add(customframe,text="Custom Pipeline")

nbook.pack( side = LEFT,fill=BOTH,padx=10,pady=10,expand=YES )

# projframe = LabelFrame(top,text="Project Json",fg=textLightColor,bg=baseColor)
# projframe.pack( side = LEFT,fill=BOTH,padx=10,pady=10,expand=YES )
# 
# filesframe = LabelFrame(top,text="Job Status",fg=textLightColor, bg=baseColor)
# filesframe.pack( side = LEFT,fill=BOTH,padx=10,pady=10,expand=YES )
# 
# 

#####################
#The Manual
#####################

# manquery=StringVar()

# E = Entry(manualframe, textvariable=manquery, fg=widgetFgColor,bg=widgetBgColor)
# E.pack(side=TOP,padx=10,pady=10)

# manscrollbar = Scrollbar(manualframe)
# manscrollbar.pack( side = RIGHT,fill=BOTH,padx=10,pady=10,expand=YES )
 
# mandisplay = Text(manualframe,width=80,height=38,bg=statusBgColor,fg=statusFgColor,font=("nimbus mono","9"),yscrollcommand = manscrollbar.set)
# manpage=open("pipeliner.manual","r").read()
# photo=PhotoImage(file='./pipeliner-logo.gif')
# mandisplay.insert(END,'\n')
# mandisplay.image_create(END, image=photo)
 
# mandisplay.insert(INSERT, manpage)

# mandisplay.insert(END, "link", ("a", "href"+"http://hpc.nih.gov"))
# mandisplay.mark_set("insert", "1.0")
# mandisplay.pack( side = LEFT,fill=BOTH,padx=10,pady=10,expand=YES )
# mandisplay.configure(state="disabled")
# manscrollbar['command']=mandisplay.yview


# def searchmanual(*args):
#     global manquery
#     global mandisplay
#     query=manquery.get()
#     mandisplay.tag_delete("red")
#     mandisplay.tag_configure("red", foreground="#ff0000")
#     F="forwards" 
#     if F=="forwards":
#         start=mandisplay.index(INSERT)+"+1c"
#         end=mandisplay.index(END)
#         forwards=True
#         backwards=False

#     else:
#         start=mandisplay.index(INSERT)+"-1c"
#         end="1.0"
#         forwards=False
#         backwards=True

#     count=IntVar()
# #    tkinter.messagebox.showinfo("",start)

#     pos = mandisplay.search(query, start,stopindex=end,forwards=forwards,backwards=backwards,nocase=1,count=count)
#     if pos:
#         mandisplay.see(pos)
#         mandisplay.mark_set("insert", pos)
#         mandisplay.mark_set("matchStart", pos)
#         mandisplay.mark_set("matchEnd", "%s+%sc" % (pos, count.get()))
#         mandisplay.tag_add("red", "matchStart", "matchEnd")
#     return
  
# manquery.trace('w',searchmanual)
 
# Bu = Button(manualframe, text="Next", fg=widgetFgColor,bg=widgetBgColor,command=searchmanual)
# Bu.pack(side=TOP,padx=5,pady=10)

#backBu = Button(manualframe, text="back", fg=widgetFgColor,bg=widgetBgColor,command=lambda: searchmanual(F="backwards"))
#backBu.pack(side=TOP,padx=5,pady=10)
  
#####################
#The Run Frame
#####################

button = Button(runframe, text="Initialize Working Directory", fg="white",bg="red",command=initialize_results)
button.grid(row=0,column=0,padx=10,pady=10,sticky="w")
L=Label(runframe,text="The Working Directory must exist. Any data within it will be deleted.",anchor="ne",bg="red",fg="white")
L.grid(row=0,column=1,padx=10,pady=10,sticky="w")

initLock=StringVar()
initlock = Checkbutton(runframe, text="Unlock",variable=initLock,bg="red",fg="white",offvalue="locked",onvalue="unlocked",state="active")
initLock.set("locked")
#initlock.grid(row=0,column=2,padx=10,pady=10,sticky="w")


button = Button(runframe, text="Symbolically Link Data", fg=widgetFgColor,bg=widgetBgColor,command=lambda:symlink(datapath.get()))
button.grid(row=1,column=0,padx=10,pady=10,sticky="w")
L=Label(runframe,text="The starting data file will be linked from the Data Directory to the Working Directory.",anchor="ne",fg=widgetFgColor,bg=widgetBgColor)
L.grid(row=1,column=1,padx=10,pady=10,sticky="w")

button = Button(runframe, text="Perform Dry Run", fg=widgetFgColor,bg=widgetBgColor,command=lambda: cmd("--dryrun -s %s/Snakefile --rerun-incomplete -d %s" % (workpath.get(),workpath.get())))
button.grid(row=2,column=0,padx=10,pady=10,sticky="w")
L=Label(runframe,text="A test will be performed to verify that the pipeline will run as expected. No data is processed at this point.",anchor="ne",bg=widgetBgColor,fg=widgetFgColor)
L.grid(row=2,column=1,padx=10,pady=10,sticky="w")

button = Button(runframe, text="Save Project Profile", fg=widgetFgColor,bg=widgetBgColor,command=save_project)
button.grid(row=3,column=0,padx=10,pady=10,sticky="w")
L=Label(runframe,text="The pipeline parameters, project information and comments are saved.",anchor="ne",bg=widgetBgColor,fg=widgetFgColor)
L.grid(row=3,column=1,padx=10,pady=10,sticky="w")

button = Button(runframe, text="Submit Job to Biowulf2 via Slurm", fg=widgetFgColor,bg=widgetBgColor,command=runslurm)
button.grid(row=4,column=0,padx=10,pady=10,sticky="w")
L=Label(runframe,text="The pipeline job is submitted to the batch queuing system on the Biowulf2 cluster.",anchor="ne",bg=widgetBgColor,fg=widgetFgColor)
L.grid(row=4,column=1,padx=10,pady=10,sticky="w")


#####################
#The Files Pane
#####################

# filescrollbar = Scrollbar(filesframe)
# filescrollbar.pack( side = RIGHT, fill=Y )

# filedisplay = Text(filesframe,width=70,height=34,bg=statusBgColor,fg=statusFgColor,font=("nimbus mono","9"),yscrollcommand = filescrollbar.set)
# filedisplay.pack(side=RIGHT,expand=YES)

# filescrollbar['command']=filedisplay.yview

# errorbut = Button(filesframe, bd =2, width=12, bg="darkgreen",fg="white",text="No Errors",command=lambda: show(errorreport,"Errors",errorsFgColor,errorsBgColor,50,100))    
# errorbut.pack( side = BOTTOM,padx=5)


# def filestats():
# #def finderrors():
#     FNULL = open(os.devnull, 'w')
#     global errorreport
#     errorreport=""
#     cmd="grep 'Error' %s/Reports/snakemake.log %s/Reports/makeasnake.log"%(workpath.get(),workpath.get())
#     try:
#         p = Popen(cmd, stderr=FNULL,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
#         Out = p.communicate()
#         out=Out[0].decode("utf-8")
#         errorreport=errorreport+out
#     except:
#         pass
#     if errorreport !="":
#         errorbut.config(bg="indianred",text="There are Errors",command=lambda:show(errorreport,"Errors",errorsFgColor,errorsBgColor,50,100))
#     else:
#         errorbut.config(bg="darkgreen",text="No Errors",command=lambda: donothing())
# #    threading.Timer(60, finderrors).start()
    
# #def filestats():
#     report=""
#     total=1
#     auto_check=float(AC.get())
# #    tkinter.messagebox.showinfo("auto_check",str(auto_check))
#     cmd="date"
#     try:
#         p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=FNULL, close_fds=True)
#         Out = p.communicate()
#         out=Out[0].decode("utf-8")
#         report=report+"Pipeline: "+Pipeline.get()+"\n"+out
#     except:
#         pass

#     for FL in ['bam','sam','vcf','gvcf','fin.bam','recal.bam','realign.bam','dedup.bam','sorted.bam','pileup.bam','g.vcf']:
#         cmd="ls -l "+workpath.get()+"/*."+FL+"|awk '{print $5}'| paste -sd+ | bc -q"
# #        cmd="ls -l ../*."+FL 
#         try:
#             p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=FNULL,close_fds=True)
#             Out = p.communicate()
#             out=Out[0].decode("utf-8")
#         except:
#             pass
#         if str(out).find("No such file")>=0:    
#             out="0"
        
#         cmd="ls -l "+workpath.get()+"/*."+FL+"|wc|awk '{print $1}'"
#         try:
#             p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=FNULL,close_fds=True)
#             Out = p.communicate()
#             out2=Out[0].decode("utf-8")

#         except:
#             pass

        
#         report=report+FL+":\t"+out.strip()+"\tbytes in "+out2.strip()+" files\n"

#     cmd="du -s ../|cut -f1"
#     try:
#         p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=FNULL,close_fds=True)
#         Out = p.communicate()
#         total=float(Out[0].decode("utf-8"))
# #        tkinter.messagebox.showinfo("Disk",str(total))
#     except:
#         pass



#     color=["blue","blue","blue","blue","blue","blue","blue","blue","blue","green","orange","red","red","red"]
#     btn[1].config(text="Disk: {0:.2f}".format(total/10**6)+" G",fg="white",bg=color[int(math.log10(total))])
# #     if runmode==1:
# #         cmd="qstat|grep "+euser.get()
# #     else:
# #         cmd="sjobs|grep "+euser.get()
# #     
# #     try:
# #         p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
# #         Out = p.communicate()
# #         out=Out[0].decode("utf-8")
# #         report=report+"\n"+out
# #     except Exception as e:
# #         tkinter.messagebox.showinfo("Error",str(e))
# # 

# #    if runmode==1:
# #        cmd="jobload "+euser.get()
# #    else:
#     cmd="jobload -u "+euser.get()
#     try:
#         p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
#         Out = p.communicate()
#         out=Out[0].decode("utf-8")
#         try:
#             pat=re.compile("User Load Avg:.+?(\d+\.\d+)%")
#             jobload=float(pat.findall(out)[0])
#             color=["black","cyan4","green4","green4","green4","green4","green4","orange4","orange4","orange4","red4"]
#             btn[0].config(text="Load: "+str(jobload)+" %",fg="white",bg=color[int((jobload+9.9)/10)])
#         except:
#             jobload="No Load"
#             btn[0].config(text=str(jobload),fg="white",bg="black")

        
#         report=report+out
#     except Exception as e:
#         tkinter.messagebox.showinfo("Error",str(e))
        
#     filedisplay.delete("1.0",END)    
#     filedisplay.insert(INSERT, report)
#     F=open(workpath.get()+"/Reports/status.log","a")
#     F.write(report)
#     F.close()
# #    FNULL.close()

#     cmd="grep ') done' %s/Reports/snakemake.log"%(workpath.get())

#     p = Popen(cmd, stderr=FNULL,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
#     Out = p.communicate()
#     out=Out[0].decode("utf-8")
# #    tkinter.messagebox.showinfo("",out)        
#     try:
#         percent=out.split("\n")[-2]
#         percent_done.set(percent)
#     except:
#         pass
    
#     if auto_check>=15:
#         threading.Timer(auto_check, filestats).start()         

# def status_auto_check_on():
#     AC.set("60")
#     #auto_check=1
#     filestats()
#     #btn[9].config(bg="red",text="Auto Check ON")
    
    
# def status_auto_check_off():
#     #btn[9].config(bg="black",fg="white",text="Auto Check OFF")
#     #global auto_check
#     AC.set("0")

# percent_done=StringVar()
# percent_done.set("0% done")
# Entry(filesframe,textvariable=percent_done,width=15,relief=RAISED).pack(padx=5,side=TOP)

# #filestats()
# button = Button(filesframe, text="Update Now", width=12, fg=widgetFgColor,bg=widgetBgColor,command=filestats)
# button.pack( side = BOTTOM,padx=5)

# btn=[0 for x in range(9)]
# for i in range(9):
#     btn[i] = Label(filesframe, bd =2, width=15, bg=baseColor,text="")    
#     #btn[i] = Button(filesframe, text="", fg="gray80",bg="gray80")
#     btn[i].pack( side = BOTTOM,padx=5)

# AC=StringVar()


# label=Label(filesframe,text="Auto Update",fg=widgetFgColor,bg=widgetBgColor)
# label.pack(side=TOP,pady=5,padx=5,anchor=W)

# Radiobutton(filesframe, text="Auto On", bg=widgetBgColor, fg=widgetFgColor,variable=AC, value=60,command=lambda: status_auto_check_on()).pack(side=TOP,anchor=W,pady=2,padx=5)
# Radiobutton(filesframe, text="Auto Off", bg=widgetBgColor, fg=widgetFgColor, variable=AC, value=0,command=lambda: status_auto_check_off()).pack(side=TOP,anchor=W,pady=5,padx=5)

# L=Label(filesframe,text="Update rate (secs)",anchor="ne",bg=widgetBgColor,fg=widgetFgColor)
# E = Entry(filesframe, bd =2, width=5,bg=entryBgColor,fg=entryFgColor,textvariable=AC)
# AC.set("0")
# L.pack(side=TOP,pady=5,padx=5)
# E.pack(side=TOP,pady=5,padx=5)







######################
#The Project Entry Form
######################
BS=StringVar()

eprojectid= StringVar()
eorganism=StringVar()
eanalyst=StringVar()
epoc=StringVar()
epi=StringVar()
euser=StringVar()
eplatform=StringVar()
eplatform.set("Illumina")
technique=StringVar()

projpanel1 = LabelFrame(editframe,text="Project Information",fg=textLightColor,bg=baseColor)
projpanel1.pack( side = LEFT,fill=BOTH,padx=10,pady=10,expand=YES )
projpanel2 = LabelFrame(editframe,text="Project Description",fg=textLightColor,bg=baseColor)
projpanel2.pack( side = LEFT,fill=BOTH,padx=10,pady=10,expand=YES )




Dscrollbar = Scrollbar(projpanel2)
Dscrollbar.grid(row=0,column=4,rowspan=40)
description = Text(projpanel2,width=50,height=38,bg=commentBgColor,fg=commentFgColor,font=("nimbus mono bold","10"),yscrollcommand = Dscrollbar.set)
Dscrollbar['command']=description.yview
description.insert(INSERT, "Enter project Description and Notes here.")
description.bind('<FocusOut>',lambda _:makejson())
description.grid(row=2,column=3,sticky="e",padx=10,pady=10)



L=Label(projpanel1,text="Project Id",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20,bg=entryBgColor,fg=entryFgColor,textvariable=eprojectid)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=0,column=0,sticky=W,padx=10,pady=10)
E.grid(row=0,column=2,sticky=W,padx=10,pady=10)

eprojectid.trace('w', makejson)

L=Label(projpanel1,text="Biowulf Username",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=euser)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=1,column=0,sticky=W,padx=10,pady=10)
E.grid(row=1,column=2,sticky=W,padx=10,pady=10)

euser.trace('w', makejson)

L=Label(projpanel1,text="Organism",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=eorganism)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=2,column=0,sticky=W,padx=10,pady=10)
E.grid(row=2,column=2,sticky=W,padx=10,pady=10)
eorganism.trace('w', makejson)

L=Label(projpanel1,text="Analyst",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=eanalyst)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=3,column=0,sticky=W,padx=10,pady=10)
E.grid(row=3,column=2,sticky=W,padx=10,pady=10)

eanalyst.trace('w', makejson)

L=Label(projpanel1,text="Point of Contact",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=epoc)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=4,column=0,sticky=W,padx=10,pady=10)
E.grid(row=4,column=2,sticky=W,padx=10,pady=10)
epoc.trace('w', makejson)

L=Label(projpanel1,text="PI",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=epi)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=5,column=0,sticky=W,padx=10,pady=10)
E.grid(row=5,column=2,sticky=W,padx=10,pady=10)
epi.trace('w', makejson)


L=Label(projpanel1,text="Platform",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=eplatform)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=5,column=0,sticky=W,padx=10,pady=10)
E.grid(row=5,column=2,sticky=W,padx=10,pady=10)
eplatform.trace('w', makejson)

L=Label(projpanel1,text="Technique",anchor="ne",bg=baseColor,fg=textLightColor)
E = Entry(projpanel1, bd =2, width=20, bg=entryBgColor, fg=entryFgColor,textvariable=technique)
#L.pack(pady=5,padx=5)
#E.pack(pady=1,padx=5)
L.grid(row=6,column=0,sticky=W,padx=10,pady=10)
E.grid(row=6,column=2,sticky=W,padx=10,pady=10)
technique.trace('w', makejson)


#logo = Canvas(editframe)
#logo.config(width=200,height=300)
#img = PhotoImage(file="pipeliner-logo.gif")
#logo.create_image(0,0,image=img, anchor="nw")
#logo.image=img
#logo.grid(row=0,column=3)



Commentscrollbar = Scrollbar(pastewriteframe)
Commentscrollbar.grid(row=0,column=4,rowspan=30)
comments = Text(pastewriteframe,width=80,height=30,bg=commentBgColor,fg=commentFgColor,font=("nimbus mono bold","11"),yscrollcommand = Commentscrollbar.set)
Commentscrollbar['command']=comments.yview
comments.insert(INSERT, "Enter data here.  Then save in one of the standard formates listed.")
#comments.bind('<FocusOut>',lambda _:writepaste())
comments.grid(row=1,column=0,sticky="e")

L=Label(pastewriteframe, text="File Type",fg=textLightColor,bg=baseColor)
L.grid(row=11,column=1)

ftypes=['pairs','rg.tab','contrasts.tab']
ftype = StringVar()
ftype.set(ftypes[0])
om = OptionMenu(pastewriteframe, ftype, *ftypes, command=writeheader)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=11,column=2,sticky="e",padx=10,pady=10)

but1 = Button(pastewriteframe, text="Save", width=12, fg=widgetFgColor,bg=widgetBgColor,command=writepaste)
but1.grid(row=11,column=3)


Projscrollbar = Scrollbar(projframe)
Projscrollbar.grid(row=2,column=4,rowspan=30 )
jsonconf = Text(projframe,width=80,height=30,bg=projectBgColor,fg=projectFgColor,font=("nimbus mono bold","11"),yscrollcommand = Projscrollbar.set)
Projscrollbar['command']=jsonconf.yview
jsonconf.insert(INSERT, projectjson)
#jsonconf.pack(side=RIGHT,expand=YES)
jsonconf.grid(row=3,column=0,sticky="e",padx=10,pady=10)
#jsonconf.bind('<FocusOut>',lambda _:saveproject(jsonconf.get("1.0",END)))
jsonconf.bind('<FocusOut>',lambda _:makejson())

datapath=StringVar()
workpath=StringVar()

datapathL = Label(opts2, text="Full Path to Data Directory",fg=textLightColor,bg=baseColor)
#datapathL.pack(side=LEFT,padx=5)
datapathL.grid(row=5,column=0,sticky=W,padx=10,pady=10)

datapathE = Entry(opts2, bd =2, width=60, bg=entryBgColor,fg=entryFgColor,textvariable=datapath)
#datapathE.pack(side=LEFT)
datapathE.grid(row=5,column=1,sticky=W,padx=10,pady=10)
datapath.trace('w', makejson)

workL = Label(opts2, text="Full Path to Working Directory",fg=textLightColor,bg=baseColor)
#workL.pack(side=LEFT,padx=5,pady=5)
workL.grid(row=6,column=0,sticky=W,padx=10,pady=10)

workE = Entry(opts2, bd =2, width=60, bg=entryBgColor, fg=entryFgColor,textvariable=workpath)
workE.pack(side=LEFT,pady=5)
#workpath.set(defaultwork)
workE.grid(row=6,column=1,sticky=W,padx=10,pady=10)
workpath.trace('w', makejson)

#########################
# The RNASeq Pane
#########################

rframe = LabelFrame(rnaseqframe,text="Star",fg=textLightColor,bg=baseColor)
rframe.pack( side = TOP,fill=X,padx=10,pady=10,expand=YES )
# 
r2frame = LabelFrame(rnaseqframe,text="DEseq2",fg=textLightColor, bg=baseColor)
r2frame.pack( side = TOP,fill=X,padx=10,pady=10,expand=YES )

r3frame = LabelFrame(rnaseqframe,text="EdgeR",fg=textLightColor, bg=baseColor)
r3frame.pack( side = BOTTOM,fill=X,padx=10,pady=10,expand=YES )

r4frame = LabelFrame(rnaseqframe,text="LimmaVoom",fg=textLightColor, bg=baseColor)
r4frame.pack( side = BOTTOM,fill=X,padx=10,pady=10,expand=YES )



rnaseqopt1s=['hg19','mm10']
rnaseqopt1 = StringVar()
rnaseqopt1.set(rnaseqopt1s[0])
om = OptionMenu(rframe, rnaseqopt1, *rnaseqopt1s, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=2,column=1,sticky=W,padx=10,pady=10)

rnaseqopt2s=['hg19','mm10']
rnaseqopt2 = StringVar()
rnaseqopt2.set(rnaseqopt2s[0])
om = OptionMenu(r2frame, rnaseqopt2, *rnaseqopt2s, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=2,column=1,sticky=W,padx=10,pady=10)

rnaseqopt3s=['hg19','mm10']
rnaseqopt3 = StringVar()
rnaseqopt3.set(rnaseqopt2s[0])
om = OptionMenu(r3frame, rnaseqopt3, *rnaseqopt3s, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=2,column=1,sticky=W,padx=10,pady=10)

rnaseqopt4s=['hg19','mm10']
rnaseqopt4 = StringVar()
rnaseqopt4.set(rnaseqopt4s[0])
om = OptionMenu(r4frame, rnaseqopt4, *rnaseqopt4s, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=2,column=1,sticky=W,padx=10,pady=10)

rnaseqcheck1=StringVar()
rnaseqcheckb1 = Checkbutton(r4frame, text="check1",variable=rnaseqcheck1,bg=baseColor,fg=textLightColor,offvalue="no",onvalue="yes")
#bySample.pack(side=LEFT,padx=5,pady=5)
rnaseqcheckb1.grid(row=0,column=5,sticky=W,padx=10,pady=10)
rnaseqcheck1.set("no")
rnaseqcheck1.trace('w', makejson)


#########################
# The Option Menus
#########################

pipeline=["initialqc","bam2recal","wgslow","exomeseq-somatic","exomeseq-germline","exomeseq-germline-recal","exomeseq-germline-partial","custom"]

filetypes=['fastq','fastq.gz','R1.trimmed.fastq.gz','fastq.bz2','trimmed.fastq.bz2','sorted.bam','dedup.bam','fin.bam','recal.bam','realign.bam','bam','bai','sam','combined.gvcf','all.snp.dbnsfp.vcf','R1_fastqc.html','qualimapReport']
filetype = StringVar()
filetype.set(filetypes[1])


label=Label(customframe,text="Initial File Type",fg=textLightColor,bg=baseColor)
#label.pack(side=LEFT,padx=5)
label.grid(row=0,column=0,sticky=W,padx=10,pady=10)

om = OptionMenu(customframe, filetype, *filetypes, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=0,column=1,sticky=W,padx=10,pady=10)


label=Label(opts2,text="Software Set",fg=textLightColor,bg=baseColor)
#label.pack(side=LEFT,padx=5,pady=5)
label.grid(row=1,column=0,sticky=W,padx=10,pady=10)

binsets=['standard-bin']
binset = StringVar()
binset.set(binsets[0])
om = OptionMenu(opts2, binset, *binsets, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=1,column=1,sticky=W,padx=10,pady=10)

label=Label(opts2,text="Annotation Set",fg=textLightColor,bg=baseColor)
#label.pack(side=LEFT,padx=5,pady=5)
label.grid(row=2,column=0,sticky=W,padx=10,pady=10)

annotations=['hg19','mm10']
annotation = StringVar()
annotation.set(annotations[0])
om = OptionMenu(opts2, annotation, *annotations, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=2,column=1,sticky=W,padx=10,pady=10)

L=Label(opts2, text="Pipeline Group",fg=textLightColor,bg=baseColor)
L.grid(row=3,column=0)

pfamilys=['exomeseq','rnaseq','mirseq','chipseq','custom']
pfamily = StringVar()
pfamily.set(pfamilys[0])
om = OptionMenu(opts2, pfamily, *pfamilys, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=3,column=1,sticky=W,padx=10,pady=10)

# label=Label(opts2,text="Resource Use",fg=textLightColor,bg=baseColor)
# label.pack(side=LEFT,padx=5,pady=5)

# clusterparams=['medium','low','high']
# cluster = StringVar()
# cluster.set(clusterparams[0])
# om = OptionMenu(opts2, cluster, *clusterparams, command=makejson)
# om.config(bg = widgetBgColor,fg=widgetFgColor)
# om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
# om.pack(side=LEFT,padx=20,pady=5)
# 

label=Label(exomeseqframe,text="Pipeline",fg=textLightColor,bg=baseColor)
#label.pack(side=LEFT,padx=5,pady=5)
label.grid(row=3,column=0,sticky=W,padx=10,pady=10)

Pipeline = StringVar()
Pipeline.set(pipeline[0])
om = OptionMenu(exomeseqframe, Pipeline, *pipeline, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=3,column=1,sticky=W,padx=10,pady=10)



efiletype = StringVar()
efiletype.set(filetypes[0])



label=Label(customframe,text="Final File Type",fg=textLightColor,bg=baseColor)
#label.pack(side=LEFT,padx=5,pady=5)
label.grid(row=4,column=0,sticky=W,padx=10,pady=10)

om = OptionMenu(customframe, efiletype, *filetypes, command=makejson)
om.config(bg = widgetBgColor,fg=widgetFgColor)
om["menu"].config(bg = widgetBgColor,fg=widgetFgColor)
#om.pack(side=LEFT,padx=20,pady=5)
om.grid(row=4,column=1,sticky=W,padx=10,pady=10)

# freezeunits=StringVar()
# freezeUnits = Checkbutton(opts2, text="Freeze Units",variable=freezeunits,bg=baseColor,fg=textLightColor,offvalue="no",onvalue="yes")
# freezeUnits.pack(side=LEFT,padx=5,pady=5)
# freezeunits.set("no")
# 
# bysample=StringVar()
# bySample = Checkbutton(opts2, text="bysample",variable=bysample,bg=baseColor,fg=textLightColor,offvalue="no",onvalue="yes")
# #bySample.pack(side=LEFT,padx=5,pady=5)
# bySample.grid(row=0,column=5,sticky=W,padx=10,pady=10)
# bysample.set("no")
# bysample.trace('w', makejson)

makejson()
top.mainloop()


