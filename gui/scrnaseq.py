
import sys,os,math,time
import json,re
import contextlib
import webbrowser
import time,threading

# from subprocess import Popen, PIPE, STDOUT, check_output
from subprocess import Popen, PIPE, STDOUT, check_output, DEVNULL
from glob import glob
from shutil import copytree

import tkinter as tk
from tkinter import Toplevel
from tkinter import Tk, END, StringVar, LEFT, TOP, BOTTOM, X, BOTH, YES, INSERT, W, E
from tkinter import Text, Entry, OptionMenu

from tkinter import ttk
from tkinter.ttk import Label, Button, LabelFrame, Scrollbar, Frame, Notebook, Style

from tkinter.filedialog import askopenfilename, asksaveasfilename, askdirectory
from tkinter.messagebox import showerror, showinfo, showwarning

from pathlib import Path
from os.path import join, exists
from os import makedirs, listdir


from gui.frame import PipelineFrame
from gui.utils import ftypes, filetypes, filetype, efiletype, batchsize, customRules
from gui.utils import USER_HOME, PIPELINER_HOME, PIPELINER_CONF


class scRNASeqFrame( PipelineFrame ) :
    def __init__(self, pipepanel, pipeline_name, *args, **kwargs) :
        PipelineFrame.__init__(self, pipepanel, pipeline_name, *args, **kwargs)
        self.info = None
        
        eframe = self.eframe = LabelFrame(self,text="Options") 
        #,fg=textLightColor,bg=baseColor)
        eframe.grid( row=5, column=1, sticky=W, columnspan=7, padx=10, pady=5 )
        
        label = Label(eframe,text="Pipeline")#,fg=textLightColor,bg=baseColor)
        label.grid(row=3,column=0,sticky=W,padx=10,pady=5)
        PipelineLabels=["Initial QC","Differential Expression"]
        Pipelines=["scrnaQC","scrnaDE"]

        self.label2pipeline = { k:v for k,v in zip(PipelineLabels, Pipelines)}
        
        PipelineLabel = self.PipelineLabel = StringVar()
        self.Pipeline = StringVar()

        PipelineLabel.set(PipelineLabels[0])        
        om = OptionMenu(eframe, PipelineLabel, *PipelineLabels, command=self.option_controller)
        om.config()#bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=10,pady=5)


######## QUALITY CONTROL/FILTERING/CLUSTERING ###############        
        self.qcOptions = qcOptions = LabelFrame( eframe, text = "QC, Filtering, and Initial Clustering")

        #algorithm selection
        self.clustAlg = clustAlg = StringVar()
        clustAlgLabel = Label(qcOptions, text="Clustering Algorithm: ")
        clustAlgDropdown = ["SLM (Smart Local Moving)", "Louvain (Original)","Louvain (with Multilevel Refinement)"]
        clustAlg.set(clustAlgDropdown[0])
        clustAlgOptMenu = OptionMenu(qcOptions, clustAlg, *clustAlgDropdown)
        clustAlgLabel.grid(row=8,column=1,sticky=W,padx=10,pady=5)
        clustAlgOptMenu.grid(row=8,column=2,sticky=W,padx=0,pady=5)
        
        #clustering resolutions
        self.resolution = resolution = StringVar()
        resolution.set("0.4,0.6,0.8,1.0,1.2")
        resolutionLabel = Label(qcOptions,text="Clustering Resolution(s): \nSeparate with commas")
        resolutionEntry = Entry(qcOptions,bd=2,width=25,textvariable=resolution)
        resolutionLabel.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        resolutionEntry.grid(row=9,column=2,sticky=W,padx=0,pady=5)
        
        #annotation db: human: HPCA/BP Encode; mouse: immgen/mouseRNAseq
        self.annotDB = annotDB = StringVar()
        annotLabel = Label(qcOptions,text="Annotation database: ")
        annotHuman = ["Human Primary Cell Atlas","Blueprint/ENCODE","Monaco Immune","Database of Immune Cell Expression (DICE)"]
        annotMouse = ["ImmGen","Mouse RNASeq"]
        annotDropdown = []
        genome = self.global_info.annotation.get()
        if (genome == "GRCh38"):
            annotDropdown = annotHuman.copy()
        elif (genome == "mm10"):
            annotDropdown = annotMouse.copy()
        else:
            annotDropdown.append("No genome selected")
        annotDB.set(annotDropdown[0])
        annotDB_om = OptionMenu(qcOptions,annotDB,*annotDropdown)
        annotLabel.grid(row=10,column=1,sticky=W,padx =10, pady =5)
        annotDB_om.grid(row=10,column=2,sticky=W,padx=0,pady=5)
        
        self.citeseq = citeseq = StringVar()
        citeseqLabel = Label(qcOptions, text = "CITESeq Included: ")
        citeseqDropdown = ["Yes","No"]
        citeseq.set(citeseqDropdown[1])
        citeseqOptMenu = OptionMenu(qcOptions, citeseq, *citeseqDropdown)
        citeseqLabel.grid(row=11,column=1,sticky=W,padx=10,pady=5)
        citeseqOptMenu.grid(row=11,column=2,sticky=W,padx=0,pady=5)

        #set groups, mostly for relabeling purposes
        self.add_info( qcOptions ) #Position 13
        # self.option_controller()

        # groups_buttonL = Label(qcOptions, text="Sample Information: ")
        # groups_button = Button(qcOptions, 
        #                                     text="Set Groups", 
        #                                     command = self.popup_groups )
        # groups_buttonL.grid(row=12,column=1,sticky=W,padx=10,pady=5)
        # groups_button.grid(row=12,column=2,sticky=W,padx=0,pady=5)
        #
        #
######### DIFFERENTIAL EXPRESSION #########
        self.deOptions = deOptions = LabelFrame( eframe, text = "Differential Gene Expression")
        
        #option for which object to run DE on
        self.rdsObject = rdsObject = StringVar()
        rdsLabel = Label(deOptions, text = "Use the pre- or post-batch corrected data:")
        rdsDropdown = ["Merged (Pre-batch correction)","Integrated (Post-batch correction)","Both"]
        rdsObject.set(rdsDropdown[2])
        rdsOptMenu=OptionMenu(deOptions,rdsObject,*rdsDropdown)
        rdsLabel.grid(row=8,column=1,sticky=W,padx=10,pady=5)
        rdsOptMenu.grid(row=8,column=2,sticky=W,padx=0,pady=5)
        
        #option for cluster resolution for DE
        self.resolutionDE = resolutionDE = StringVar()
        resolutionDE.set("0.4,0.6,0.8,1.0,1.2")
        resolutionDELabel = Label(deOptions,text="Clustering Resolution: \nChoose a previous resolution or \nselect a new resolution to run")
        resolutionDEEntry = Entry(deOptions,bd=2,width=25,textvariable=resolutionDE)
        resolutionDELabel.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        resolutionDEEntry.grid(row=9,column=2,sticky=W,padx=0,pady=5)
        
        #options for difftest
        self.testDE = testDE = StringVar()
        testDELabel = Label(deOptions, text = "Statistical test for differential expression:")
        testDEDropdown = ["MAST","DESeq2","Likelihood Ratio","Logistic regression","Negative Binomial","Wilcoxon","Student's T"]
        testDE.set(testDEDropdown[0])
        testDEMenu = OptionMenu(deOptions,testDE,*testDEDropdown)
        testDELabel.grid(row=10,column=1,sticky=W,padx=10,pady=5)
        testDEMenu.grid(row=10,column=2,sticky=W,padx=0,pady=5)
        
        #options for filters
        self.deMinPct = deMinPct = StringVar()
        deMinPctLabel = Label(deOptions, text = "Minimum fraction of cells expressing DE genes:")
        deMinPct.set("0.1")
        deMinPctEntry = Entry(deOptions,bd=2,width = 10, textvariable=deMinPct)
        deMinPctLabel.grid(row=11,column=1,sticky=W,padx=10,pady=5)
        deMinPctEntry.grid(row=11,column=2,sticky=W,padx=0,pady=5)
        
        self.deMinFC = deMinFC = StringVar()
        deMinFCLabel = Label(deOptions, text = "Minimum fold change to report DE genes:")
        deMinFC.set("0.25")
        deMinFCEntry = Entry(deOptions,bd=2,width = 10, textvariable=deMinFC)
        deMinFCLabel.grid(row=12,column=1,sticky=W,padx=10,pady=5)
        deMinFCEntry.grid(row=12,column=2,sticky=W,padx=0,pady=5)
        
        #use groups and contrasts for differential expression, have options to create new groups
        self.om_groups = LabelFrame(deOptions, text="Sample Information")
        self.groups_button = Button(self.om_groups, text="Set Groups", command = self.popup_groups )
        self.groups_button.grid(row=5, column=5, padx=10, pady=5)
        self.contrasts_button = Button(self.om_groups, text="Set Contrasts", command = self.popup_groups )
        self.contrasts_button.grid(row=5, column=6, padx=10, pady=5)
        self.om_groups.grid(row=13,column=0,columnspan=6,sticky=W,padx=20,pady=10)

        #option for merged/integrated object
        # self.integrateOption = integrateOption = StringVar()
        # integrateLabel = Label(deOptions, text="Merged or Integrated (batch corrected): ")
        # integrateDropdown = ["Merged","Integrated"]
        # integrateOption.set(integrateDropdown[0])
        # integrateOptMenu = OptionMenu(deOptions, integrateOption, *integrateDropdown)
        # integrateLabel.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        # integrateOptMenu.grid(row=9,column=2,sticky=W,padx=0,pady=5)        

        self.option_controller()

### SUBFUNCTION FOR CREATING GROUPS.TAB AND CLUSTERS.TAB WINDOW ######    
    def add_info( self, parent ) :
        if not self.info :
            self.info = LabelFrame(parent, text="Sample Information")
            self.groups_button = Button(self.info, 
                                            text="Set Groups", 
                                            command = self.popup_groups )
            self.contrasts_button = Button(self.info,
                                            text="Set Contrasts",
                                            command = self.popup_contrasts )
          
            #self.pairs_load_button.pack( side=BOTTOM, padx=5, pady=5 )
            #self.pairs_save_button.pack( side=BOTTOM, padx=5, pady=5 )

            self.groups_button.grid( row=5, column=5, padx=10, pady=5 )
            self.contrasts_button.grid( row=5, column=6, padx=10, pady=5 )
            
        self.info.grid(row=13,column=0, columnspan=6, sticky=W, padx=20, pady=10)
                       
    def popup_groups( self ) :
        self.popup_window( "Groups Information", "groups.tab" )
    def popup_contrasts( self ) :
        self.popup_window( "Contrasts Information", "contrasts.tab" )
    def popup_window( self, text, filename ) :
        top = Toplevel()
  
        info = LabelFrame(top, text=text )#"Group Information")
        info_text = Text(info,
                              width=50,
                              height=8,
                              #bg=projectBgColor,
                              #fg=projectFgColor,
                              font=("nimbus mono bold","11")
                             )
        
        def savefunc() :
            self.writepaste( filename, info_text ) 
        
        def loadfunc() :
            self.readpaste( filename, info_text ) 
        
        info_save_button = Button(info, 
                                  text="Save", 
                                  command = savefunc )
        info_load_button = Button(info,
                                  text="Load",
                                  command = loadfunc )

        #self.pairs_load_button.pack( side=BOTTOM, padx=5, pady=5 )
        #self.pairs_save_button.pack( side=BOTTOM, padx=5, pady=5 )

        info_load_button.grid( row=5, column=5, padx=10, pady=5 )
        info_save_button.grid( row=5, column=6, padx=10, pady=5 )
        info_text.grid(row=1,  rowspan=3, 
                       column=1,  columnspan=7,
                       padx=5, pady=5 )

        info.grid(row=7,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
        top.focus_force()
    
    def set_count_directory( self ):
        fname = askdirectory( initialdir = USER_HOME, 
                             title="Select Data Directory")
        self.countpath.set(fname)

    def option_controller( self, *args, **kwargs ) :

        PipelineFrame.option_controller( self )

        self.Pipeline.set( self.label2pipeline[self.PipelineLabel.get()] )
        print( self.Pipeline.get() )

        if self.Pipeline.get() == 'scrnaQC' :
            self.deOptions.grid_forget()
            self.qcOptions.grid(row=8,column=0, columnspan=6, sticky=W, padx=20, pady=10 )

        elif self.Pipeline.get() == 'scrnaDE' :
            self.deOptions.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )
            self.qcOptions.grid_forget()

            
    def makejson_wrapper( self, *args, **kwargs ) :
        self.makejson(*args, **kwargs)


    def makejson(self, *args):
        ###################################################### FROM RNASEQ.PY, NOT SURE WHAT ALL THIS DOES #############################3
        #print(args[0])
        caller=args[0]
        #global PD
        #global UnitsBak
        #global RGbak
        D=dict()
        try:
            F=open(self.workpath.get()+"/samples","r")
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
            F=open(self.workpath.get()+"/pairs","r")
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
        try:
            F=open(self.workpath.get()+"/contrasts.tab","r")
    #        f=F.read().split('\n')
    #        F.close()
    #        D["rsamps"]=f[0].split()
    #        D["rgroups"]=f[1].split()
    #        D["rcontrasts"]=f[2].split()
    #        D["rlabels"]=f[3].split()        
            f=F.readlines()
            F.close()
    ##        sampl=[]
    ##        grp=[]
            cont=[]
    ##        lbl=[]
            for x in f:
                  if len(x.split()) == 2:
                     cont.append(x.split()[0])
                     cont.append(x.split()[1])
            D["rcontrasts"]=cont
    #        contrasts["rcontrasts"]=cont
            contrasts=D
        except:
            contrasts={"rcontrasts":"na"}
    ##        contrasts={"rsamps":"na","rgroups":"na","rcontrasts":"na","rlabels":"na"}   
    ##------
        D=dict()
        try:
            F=open(self.workpath.get()+"/groups.tab","r")
            f=F.readlines()
            F.close()
            sampl=[]
            grp=[]
    #        cont=[]
            lbl=[]
            for x in f:
    #           if len(x.split()) == 4 or len(x.split()) == 3:
                if len(x.split()) == 3:  
                    sampl.append(x.split()[0])
                    grp.append(x.split()[1])
                    lbl.append(x.split()[2])
            D["rsamps"]=sampl
            D["rgroups"]=grp
            D["rlabels"]=lbl
    #        D["rcontrasts"]="na"
    #        contrasts=D
            groups=D
        except:
    #        contrasts={"rsamps":"na","rgroups":"na","rcontrasts":"na","rlabels":"na"}
            groups={"rsamps":"na","rgroups":"na","rlabels":"na"}          
    ##------   
        D=dict() 
        FT = filetype#.get()
    #    p = Popen("ls "+workpath.get()+"/*."+FT, shell=True, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, close_fds=True)
        p = Popen("find "+self.workpath.get()+" -maxdepth 1 -type l -printf '%f\n' ", shell=True, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, close_fds=True)
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
            key=re.sub("\.R[12]","",key)
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
            F=open(self.workpath.get()+"/rg.tab","r")
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
    # ###########################################################################################################################

        PD=dict()

        smparams=[]
        
        algorithmDict = {"Louvain (Original)": 1, "Louvain (with Multilevel Refinement)": 2,"SLM (Smart Local Moving)":3}        
        algorithm = algorithmDict[self.clustAlg.get()]
        
        annotationDict = {"Human Primary Cell Atlas":"HPCA", "Blueprint/ENCODE":"BP_encode","Monaco Immune":"monaco","Database of Immune Cell Expression (DICE)":"dice","ImmGen":"immgen","Mouse RNASeq":"mouseRNAseq"}
        annotation = annotationDict[self.annotDB.get()]
        
        resolutionStr = self.resolution.get()
        resolutionStr = re.sub("\s",",",resolutionStr) #remove whitespaces
        resolutionStr = re.sub("[,]+",",",resolutionStr) #remove multiple commas
        
        rdsDict = {"Merged (Pre-batch correction)":"merged","Integrated (Post-batch correction)":"integrated","Both":"both"}
        rdsSelect = rdsDict[self.rdsObject.get()]

        deResolutionStr = self.resolutionDE.get()
        deResolutionStr = re.sub("\s",",",deResolutionStr) #remove whitespace
        deResolutionStr = re.sub("[,]+",",",deResolutionStr) #remove multiple commas       

        deTestDict = {"MAST":"MAST","DESeq2":"DESeq2","Likelihood Ratio":"bimod","Logistic regression":"LR","Negative Binomial":"negbinom","Wilcoxon":"wilcox","Student's T":"t"}
        testDESelect = deTestDict[self.testDE.get()]

        for i in range(len(self.parameters)):

            if cp[i].var.get()=="1":
                smparams.append(parameters[i])
                
        AD=eval( open( join(PIPELINER_HOME, 
                            self.annotation.get()+".json"
                           ), "r"
                     ).read()
               )
        gi = self.global_info
        species = ""
        if gi.annotation.get() == "GRCh38":
            species = "human"
        elif gi.annotation.get() == "mm10":
            species = "mouse"
        PD={
            'project': {
                'pfamily': gi.pfamily.get(),
                'units': units, 
                'samples': samples, 
                'pairs': pairs,
                'id': gi.eprojectid.get(), 
                'pi': gi.epi.get(), 
                'organism': gi.eorganism.get(), 
                'analyst': gi.eanalyst.get(), 
                'poc': gi.epoc.get(), 
                'pipeline': self.Pipeline.get(), 
                'version':"4.0", 
                'annotation': gi.annotation.get(), 
                'datapath': self.datapath.get(), 
                'targetspath': self.targetspath.get(), 
                'filetype': filetype , 
                'binset': "standard-bin", 
                'username': gi.euser.get(), 
                'flowcellid': gi.eflowcell.get(), 
                'platform': gi.eplatform.get(), 
                'custom': customRules, 
                'efiletype': efiletype, 
                'workpath': self.workpath.get(),
                'batchsize': batchsize, 
                "smparams": smparams,
                "rgid": RG,
                "cluster": "cluster_medium.json", 
                "description": gi.description.get('1.0',END), 
                "technique" : gi.technique.get(),
                
                "CLUSTALG" : algorithm,
                "ANNOTDB" : annotation,
                "CLUSTRESOLUTION": resolutionStr,
 #               "INTEGRATEDE": self.integrateOption.get(),
                "SPECIES": species,
                "CITESEQ": self.citeseq.get(),
                # "CRID": self.scrCRID.get(),
                # "EXPECTED": self.scrExpected.get(),
                # "COUNTSPATH": self.countpath.get(),
                # "MATTYPE": self.mattype.get(),
                # "DOCYCLEREGRESS": self.docycleregress.get(),
                # "RESOLUTION": self.scrRes.get(),
                # "PCS": self.scrPCs.get(),
                "FILEDE": rdsSelect,
                "CLUSTERDE": deResolutionStr,
                "TESTDE": testDESelect,
                "DEMINPCT" : self.deMinPct.get(),
                "DEMINFC": self.deMinFC.get(),
                "contrasts" : contrasts,
                "groups": groups
             }
        } 
        
        J=json.dumps(PD, sort_keys = True, indent = 4, ensure_ascii=True)
        gi.jsonconf.delete("1.0", END)    
        gi.jsonconf.insert(INSERT, J)
        self.saveproject(gi.jsonconf.get("1.0",END))
    
   
