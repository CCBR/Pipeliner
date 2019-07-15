
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
        PipelineLabels=["CellRanger","Initial/QC","Clustering","Multi-Sample Clustering" ]
        Pipelines=["cellranger","scrnaseqinit","scrnaseqcluster", "scrnaseqmulticluster"]

        self.label2pipeline = { k:v for k,v in zip(PipelineLabels, Pipelines)}
        
        PipelineLabel = self.PipelineLabel = StringVar()
        self.Pipeline = StringVar()

        PipelineLabel.set(PipelineLabels[0])        
        om = OptionMenu(eframe, PipelineLabel, *PipelineLabels, command=self.option_controller)
        om.config()#bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=10,pady=5)
        
        self.crOpts = crOpts = LabelFrame( eframe, 
                                              text="CellRanger Settings" )
        self.scrCRID = scrCRID = StringVar()
        scrCRID.set("SPECIFY_PREFIX_HERE")
        self.scrExpected = scrExpected = StringVar()
        scrExpected.set("3000")

        scrcridL = Label(crOpts, text="CellRanger Sample ID: ")
        scrcridE = Entry(crOpts, bd =2, width=25, textvariable=scrCRID)
        scrcridL.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        scrcridE.grid(row=9,column=2,sticky=W,padx=0,pady=5)
        
        screxpectedL = Label(crOpts, text="Expected number of cells: ")
        screxpectedE = Entry(crOpts, bd =2, width=8, textvariable=scrExpected)
        screxpectedL.grid(row=10,column=1,sticky=W,padx=10,pady=5)
        screxpectedE.grid(row=10,column=2,sticky=W,padx=0,pady=5)

        self.clusterOpts = clusterOpts = LabelFrame( eframe, 
                                              text="Clustering and tSNE Options" )

        self.scrPCs = scrPCs = StringVar()
        scrPCs.set("12")
        self.scrRes = scrRes = StringVar()
        scrRes.set("0.6")
        
        #scrPCs.trace('w', lambda a,b,c,x="scrPCs": makejson(x))

        #Filter out genes < [5] read counts in < [2] samples
        scrpcsL = Label(clusterOpts, text="Include principal components 1 through ")
        scrpcsE = Entry(clusterOpts, bd =2, width=3, textvariable=scrPCs)
        scrresL = Label(clusterOpts, text="with clustering resolution: ")
        scrresE = Entry(clusterOpts, bd =2, width=3, textvariable=scrRes)
        
        scrpcsL.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        scrpcsE.grid(row=9,column=2,sticky=W,padx=0,pady=5)
        scrresL.grid(row=9,column=3,sticky=W,padx=5,pady=5)
        scrresE.grid(row=9,column=4,sticky=W,padx=0,pady=5)
        #scrRes.trace('w', lambda a,b,c,x="scrPCs": makejson(x))
        
        clusterOpts.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )
        

        self.multiclusterOpts = multiclusterOpts = LabelFrame( eframe,
                                                text = "Multi-Sample Clustering and tSNE Options")

        scrccsL = Label(multiclusterOpts, text="Include canonical components 1 through ")
        scrccsE = Entry(multiclusterOpts, bd =2, width=3, textvariable=scrPCs)
        scrmcresL = Label(multiclusterOpts, text="with clustering resolution: ")
        scrmcresE = Entry(multiclusterOpts, bd =2, width=3, textvariable=scrRes)

        scrccsL.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        scrccsE.grid(row=9,column=2,sticky=W,padx=0,pady=5)
        scrmcresL.grid(row=9,column=3,sticky=W,padx=5,pady=5)
        scrmcresE.grid(row=9,column=4,sticky=W,padx=0,pady=5)


        self.qcOpts = qcOpts = LabelFrame( eframe, 
                                              text="Initial Settings" )
        countL = Label( qcOpts, text="Counts/Matrix Dir:" )
        countL.grid(row=9, column=1, sticky=W, padx=10, pady=5 )
        countpath=StringVar()  
        self.countpath = countpath
        count_entry = Entry(qcOpts, 
                           bd =2, 
                           width = 50, 
                           #bg = entryBgColor, 
                           #fg = entryFgColor, 
                           textvariable = countpath, state='normal'
                          )
        count_entry.grid( row=9, column=2, columnspan=3 )
        self.count_button = count_button = Button( qcOpts, 
                             text="Open Directory", 
                             command=self.set_count_directory )
        count_button.grid( row=9, column=5 )

        self.mattype = mattype = StringVar()
        mattypeL = Label(qcOpts, text="Count matrix format: ")
        scrMatTypeDropdown = ["cellranger", "cellranger_raw", "zumi", "biorad"]
        mattype.set(scrMatTypeDropdown[0])
        mattype_om = OptionMenu(qcOpts, mattype, *scrMatTypeDropdown)
        mattypeL.grid(row=10,column=1,sticky=W,padx=10,pady=5)
        mattype_om.grid(row=10,column=2,sticky=W,padx=0,pady=5)
        
        self.docycleregress = docycleregress = StringVar()
        docycleregressL = Label(qcOpts, text="Do cell cycle regression? ")
        scrCycleDropdown = ["TRUE", "FALSE"]
        docycleregress.set(scrCycleDropdown[0])
        cycle_om = OptionMenu(qcOpts, docycleregress, *scrCycleDropdown)
        docycleregressL.grid(row=11,column=1,sticky=W,padx=10,pady=5)
        cycle_om.grid(row=11,column=2,sticky=W,padx=0,pady=5)

        usecycleregressL_c = Label(clusterOpts, text="Use cell cycle regressed data? ")
        docycleregress.set(scrCycleDropdown[0])
        cycle_om_c = OptionMenu(clusterOpts, docycleregress, *scrCycleDropdown)
        usecycleregressL_c.grid(row=10,column=1,sticky=W,padx=10,pady=5)
        cycle_om_c.grid(row=10,column=2,sticky=W,padx=0,pady=5)

        usecycleregressL_mc = Label(multiclusterOpts, text="Use cell cycle regressed data? ")
        docycleregress.set(scrCycleDropdown[0])
        cycle_om_mc = OptionMenu(multiclusterOpts, docycleregress, *scrCycleDropdown)
        usecycleregressL_mc.grid(row=10,column=1,sticky=W,padx=10,pady=5)
        cycle_om_mc.grid(row=10,column=2,sticky=W,padx=0,pady=5)


        groups_buttonL = Label(qcOpts, text="SAMPLE INFORMATION: ")
        groups_button = Button(qcOpts, 
                                            text="Set Groups", 
                                            command = self.popup_groups )
        groups_buttonL.grid(row=12,column=1,sticky=W,padx=10,pady=5)
        groups_button.grid(row=12,column=2,sticky=W,padx=0,pady=5)
        #####################
        
        self.option_controller()

    def popup_groups( self ) :
        self.popup_window( "Groups Information", "groups.tab" )

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

        if self.Pipeline.get() == 'cellranger' :
            self.clusterOpts.grid_forget()
            self.crOpts.grid(row=8,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
            self.qcOpts.grid_forget()
            self.multiclusterOpts.grid_forget()
        elif self.Pipeline.get() == 'scrnaseqcluster' :
            self.clusterOpts.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )
            self.crOpts.grid_forget()
            self.qcOpts.grid_forget()
            self.multiclusterOpts.grid_forget()
        elif self.Pipeline.get() == 'scrnaseqinit' :
            self.clusterOpts.grid_forget()
            self.crOpts.grid_forget()
            self.qcOpts.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )
            self.multiclusterOpts.grid_forget()
        elif self.Pipeline.get() == 'scrnaseqmulticluster' :
            self.clusterOpts.grid_forget()
            self.crOpts.grid_forget()
            self.qcOpts.grid_forget()
            self.multiclusterOpts.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )


            
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

        for i in range(len(self.parameters)):

            if cp[i].var.get()=="1":
                smparams.append(parameters[i])
                
        AD=eval( open( join(PIPELINER_HOME, 
                            self.annotation.get()+".json"
                           ), "r"
                     ).read()
               )
        gi = self.global_info 
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
                "CRID": self.scrCRID.get(),
                "EXPECTED": self.scrExpected.get(),
                "COUNTSPATH": self.countpath.get(),
                "MATTYPE": self.mattype.get(),
                "DOCYCLEREGRESS": self.docycleregress.get(),
                "RESOLUTION": self.scrRes.get(),
                "PCS": self.scrPCs.get(),
                "groups": groups
             }
        } 
        
        J=json.dumps(PD, sort_keys = True, indent = 4, ensure_ascii=True)
        gi.jsonconf.delete("1.0", END)    
        gi.jsonconf.insert(INSERT, J)
        self.saveproject(gi.jsonconf.get("1.0",END))
    
   
