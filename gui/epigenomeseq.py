
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
from os import listdir, makedirs

from gui.frame import PipelineFrame
from gui.utils import ftypes, filetypes, filetype, efiletype, batchsize, customRules
from gui.utils import USER_HOME, PIPELINER_HOME, PIPELINER_CONF

class ChIPSeqFrame( PipelineFrame ) :
    def __init__(self, pipepanel, pipeline_name, *args, **kwargs) :
        PipelineFrame.__init__(self, pipepanel, pipeline_name, *args, **kwargs)
        self.info = None
        
        eframe = self.eframe = LabelFrame(self,text="Options") 
        #,fg=textLightColor,bg=baseColor)
        eframe.grid( row=5, column=1, sticky=W, columnspan=7, padx=10, pady=5 )
        
        label = Label(eframe,text="Pipeline")#,fg=textLightColor,bg=baseColor)
        label.grid(row=3,column=0,sticky=W,padx=10,pady=5)
        Pipelines=["ChIPSeq"]
        Pipeline = self.Pipeline = StringVar()
        Pipeline.set(Pipelines[0])
        
        om = OptionMenu(eframe, Pipeline, *Pipelines, command=self.makejson_wrapper)
        om.config()#bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=10,pady=5)

        label.grid(row=6, column=0, columnspan=5, sticky=W, padx=10, pady=5)
       
        
        self.add_info(eframe)
    
    def init_work_dir( self ) :
        #basic building!
        if PipelineFrame.init_work_dir( self ) :
            pass
        else :
            return
        
        fname = self.workpath.get()
        
        try :
            #need to be solved by making an empty dir in the Results-template
            #makedirs( join(fname, "QC") ) 
            #os.mknod can replace but OSX needs a super user prev.
            open( join(fname, "peakcallinfo.csv"), 'w' ).close() 
            #open( join(fname, "samples"), 'w' ).close()
            
            #print( "copying", 'template', "into", fname )
            #os.system( "cp -r %s/Results-template/* %s"%(PIPELINER_HOME, fname ) )
                
        except :
            showerror( "Initialization failed!", "Work directory data structure " )
            return
        
        if self.make_symlinks() :
            pass
        else :
            showerror( "Symlink failed", "" )
    
    
    def set_data_directory( self ):
        fname = askdirectory( initialdir = USER_HOME, 
                             title="Select Data Directory")
        self.datapath.set(fname)                                    
        #count number
        self.data_count['text'] = str(
            len([fn for fn in listdir(fname) if fn.endswith(filetype) or fn.endswith('.fastq')] )
        ) 
        print( "Found", self.data_count['text'], filetype, "files!" )
    
    
    def makejson_wrapper( self, *args, **kwargs ) :
        if self.Pipeline.get() == 'exomeseq-somatic' :
            self.add_pairs( self.eframe )
        elif self.Pipeline.get() != 'exomeseq-somatic' :
            self.del_pairs( self.eframe )
        self.makejson(*args, **kwargs)
    
    
    def add_info( self, parent ) :
        if not self.info :
            self.info = LabelFrame(parent, text='Peak Call Info')
            self.info_text = Text( self.info,
                             width=50,
                             height=8,
                             #bg=projectBgColor,
                             #fg=projectFgColor,
                             font=("nimbus mono bold","11")
                            )

            self.info_save_button = Button(self.info, 
                                            text="Save", 
                                            command = self.writeinfo )
            self.info_load_button = Button(self.info,
                                            text="Load",
                                            command = self.readinfo )
            
            #self.pairs_load_button.pack( side=BOTTOM, padx=5, pady=5 )
            #self.pairs_save_button.pack( side=BOTTOM, padx=5, pady=5 )

            self.info_load_button.grid( row=5, column=5, padx=10, pady=5 )
            self.info_save_button.grid( row=5, column=6, padx=10, pady=5 )
            self.info_text.grid( row=1,  rowspan=3, 
                                 column=1,  columnspan=7,
                                 padx=5, pady=5 )

        self.info.grid(row=7,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
        
            
    def writeinfo( self ) :
        self.writepaste( 'peakcallinfo.csv', self.info_text )
    
    def readinfo( self ) :
        self.readpaste( 'peakcallinfo.csv', self.info_text )
     