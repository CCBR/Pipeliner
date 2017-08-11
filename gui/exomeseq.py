import sys,os,math,time
import json,re
import contextlib
import webbrowser
import time,threading

# from subprocess import Popen, PIPE, STDOUT, check_output
from subprocess import Popen, PIPE, STDOUT, check_output, DEVNULL
from glob import glob
from shutil import copytree
from os import listdir, makedirs
from os.path import join, exists

import tkinter as tk
from tkinter import Toplevel
from tkinter import Tk, END, StringVar, LEFT, TOP, BOTTOM, X, BOTH, YES, INSERT, W, E
from tkinter import Text, Entry, OptionMenu, Button

from tkinter import ttk
from tkinter.ttk import Label, LabelFrame, Scrollbar, Frame, Notebook, Style

from tkinter.filedialog import askopenfilename, asksaveasfilename, askdirectory
from tkinter.messagebox import showerror, showinfo, showwarning

from pathlib import Path

from gui.frame import PipelineFrame
from gui.utils import ftypes, filetypes, filetype, efiletype, batchsize, customRules
from gui.utils import USER_HOME, PIPELINER_HOME, PIPELINER_CONF


class ExomeSeqFrame( PipelineFrame ) :
    def __init__(self, pipepanel, pipeline_name, *args, **kwargs) :
        PipelineFrame.__init__(self, pipepanel, pipeline_name, *args, **kwargs)
        self.pairs = None
        
        eframe = self.eframe = LabelFrame(self,text="Options") 
        #,fg=textLightColor,bg=baseColor)
        eframe.grid( row=5, column=1, sticky=W, columnspan=7, padx=10, pady=5 )
        
        label = Label(eframe,text="Pipeline")#,fg=textLightColor,bg=baseColor)
        label.grid(row=3,column=0,sticky=W,padx=10,pady=5)
        PipelineLabels = ["Initial QC", "Germline", 'Somatic Tumor-Normal', 'Somatic Tumor-Only']
        Pipelines=["initialqc", "exomeseq-germline", "exomeseq-somatic", "exomeseq-somatic-tumoronly"]
        self.label2pipeline = { k:v for k,v in zip(PipelineLabels, Pipelines)}

        Pipeline = self.Pipeline = StringVar()
        PipelineLabel = self.PipelineLabel = StringVar()
        self.Pipeline = StringVar()

        PipelineLabel.set(PipelineLabels[0])
        
        om = OptionMenu(eframe, PipelineLabel, *PipelineLabels, command=self.option_controller)
        #om.config()#bg = widgetBgColor,fg=widgetFgColor)
        #om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=10,pady=5)

        targetsL=Label(eframe,
                       text="Target Capture Kit")
                       #,fg=textLightColor,bg=baseColor)
        targetsL.grid(row=5,column=0,sticky=W,padx=10,pady=5)
        targetsE = Entry(eframe,textvariable=self.targetspath, width=50)

        if self.genome=="hg19":
            self.targetspath.set( 
                "/data/CCBR_Pipeliner/db/PipeDB/lib/SS_v5_UTRs_hg19.bed" )
        else:
            self.targetspath.set(
                "/data/CCBR_Pipeliner/db/PipeDB/lib/SureSelect_mm10.bed")

        targetsE.grid(row=5,column=1,columnspan=6,sticky=W,padx=10,pady=5)
        self.targetspath.trace('w', lambda a,b,c,x="targetspath":self.makejson(x))
        label = Label (eframe, 
                       text = 
                       "By default, the path to the Agilent V5+UTR targets file is filled in here" ) 
        
        label.grid(row=6, column=0, columnspan=5, sticky=W, padx=10, pady=5)
        
        
            
    def option_controller( self, *args, **kwargs ) :
        PipelineFrame.option_controller( self )
        self.Pipeline.set( self.label2pipeline[self.PipelineLabel.get()] )
        print( self.Pipeline.get() )
        
        if self.Pipeline.get() == 'exomeseq-somatic' :
            self.add_pairs( self.eframe )
            #self.dry_button.config( state="disabled" )
        elif self.Pipeline.get() != 'exomeseq-somatic' :
            self.del_pairs( self.eframe )
            #if self.workpath.get() :
            #    self.dry_button.config( state='active' )
                
    
    def add_pairs( self, parent ) :
        if not self.pairs :
            self.pairs = LabelFrame(parent, text='Pairs')
            self.pairs_text = Text( self.pairs,
                             width=50,
                             height=8,
                             #bg=projectBgColor,
                             #fg=projectFgColor,
                             font=("nimbus mono bold","11")
                            )

            self.pairs_save_button = Button(self.pairs, 
                                            text="Save", 
                                            command = self.writepair )
            
            self.pairs_load_button = Button(self.pairs,
                                            text="Load",
                                            command = self.readpair )
            
            #self.pairs_load_button.pack( side=BOTTOM, padx=5, pady=5 )
            #self.pairs_save_button.pack( side=BOTTOM, padx=5, pady=5 )

            self.pairs_load_button.grid( row=5, column=5, padx=10, pady=5 )
            self.pairs_save_button.grid( row=5, column=6, padx=10, pady=5 )
            self.pairs_text.grid( row=1,  rowspan=3, 
                                 column=1,  columnspan=7,
                                 padx=5, pady=5 )

        self.pairs.grid(row=7,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
        
    def del_pairs( self, parent ) :
        if self.pairs :
            self.pairs.grid_forget()
            
    def writepair( self ) :
        self.writepaste( 'pairs', self.pairs_text )
        #self.dry_button.config( state='active' )
    
    def readpair( self ) :
        self.readpaste( 'pairs', self.pairs_text )
        #self.dry_button.config( state='active' )


