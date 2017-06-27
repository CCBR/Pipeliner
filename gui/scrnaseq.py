
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
        PipelineLabels=["CellRanger","Initial/QC","Clustering" ]
        Pipelines=["cellranger","scrnaseqinit","scrnaseqcluster"]

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
        self.scrExpected = scrExpected = StringVar()
        scrExpected.set("3000")

        screxpectedL = Label(crOpts, text="Expected number of cells: ")
        screxpectedE = Entry(crOpts, bd =2, width=8, textvariable=scrExpected)

        screxpectedL.grid(row=9,column=1,sticky=W,padx=10,pady=5)
        screxpectedE.grid(row=9,column=2,sticky=W,padx=0,pady=5)
        
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
        
        self.qcOpts = qcOpts = LabelFrame( eframe, 
                                              text="Initial Settings" )
        countL = Label( qcOpts, text="Data Directory:" )
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
        #####################
        
        self.option_controller()
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
        elif self.Pipeline.get() == 'scrnaseqcluster' :
            self.clusterOpts.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )
            self.crOpts.grid_forget()
            self.qcOpts.grid_forget()
        elif self.Pipeline.get() == 'scrnaseqinit' :
            self.clusterOpts.grid_forget()
            self.crOpts.grid_forget()
            self.qcOpts.grid( row=8, column=0, columnspan=4, sticky=W, padx=20, pady=10 )


            
    def makejson_wrapper( self, *args, **kwargs ) :
        self.makejson(*args, **kwargs)


    def makejson(self, *args):
        #print(args[0])
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
                'id': gi.eprojectid.get(), 
                'pi': gi.epi.get(), 
                'organism': gi.eorganism.get(), 
                'analyst': gi.eanalyst.get(), 
                'poc': gi.epoc.get(), 
                'pipeline': self.Pipeline.get(), 
                'version':"1.0", 
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
                "smparams": smparams,
                "cluster": "cluster_medium.json", 
                "description": gi.description.get('1.0',END), 
                "technique" : gi.technique.get(),
                "EXPECTED": self.scrExpected.get(),
                "COUNTSPATH": self.countpath.get(),
                "RESOLUTION": self.scrRes.get(),
                "PCS": self.scrPCs.get()
             }
        } 
        
        J=json.dumps(PD, sort_keys = True, indent = 4, ensure_ascii=True)
        gi.jsonconf.delete("1.0", END)    
        gi.jsonconf.insert(INSERT, J)
        self.saveproject(gi.jsonconf.get("1.0",END))
    
   
