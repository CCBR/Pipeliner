
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

import tkinter as tk
from tkinter import Toplevel
from tkinter import Tk, END, StringVar, LEFT, TOP, BOTTOM, X, BOTH, YES, INSERT, W, E
from tkinter import Text, Entry, OptionMenu

from tkinter import ttk
from tkinter.ttk import Label, Button, LabelFrame, Scrollbar, Frame, Notebook, Style

from tkinter.filedialog import askopenfilename, asksaveasfilename, askdirectory
from tkinter.messagebox import showerror, showinfo, showwarning

from pathlib import Path

from gui.frame import PipelineFrame
from gui.utils import ftypes, filetypes, filetype, efiletype, batchsize, customRules
from gui.utils import USER_HOME, PIPELINER_HOME, PIPELINER_CONF


class GenomeSeqFrame( PipelineFrame ) :
    def __init__(self, pipepanel, pipeline_name, *args, **kwargs) :
        PipelineFrame.__init__(self, pipepanel, pipeline_name, *args, **kwargs)
        self.pairs = None
        
        eframe = self.eframe = LabelFrame(self,text="Options") 
        #,fg=textLightColor,bg=baseColor)
        eframe.grid( row=5, column=1, sticky=W, columnspan=7, padx=10, pady=5 )
        
        label = Label(eframe,text="Pipeline")#,fg=textLightColor,bg=baseColor)
        label.grid(row=3,column=0,sticky=W,padx=10,pady=5)
        Pipelines=["initialqcgenomeseq","wgslow"]
        Pipeline = self.Pipeline = StringVar()
        Pipeline.set(Pipelines[0])
        
        om = OptionMenu(eframe, Pipeline, *Pipelines, command=self.makejson_wrapper)
        om.config()#bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=10,pady=5)
        
        
     
    def init_work_dir( self ) :
        #basic building!
        if PipelineFrame.init_work_dir( self ) :
            pass
        else :
            showerror( "Initialization failed!", "Work directory could not be made." )
            return
        
        fname = self.workpath.get()
        try :
            #need to be solved by making an empty dir in the Results-template
            makedirs( join(fname, "QC") ) 
            #os.mknod can replace but OSX needs a super user prev.
            #open( join(fname, "pairs"), 'w' ).close() 
            #open( join(fname, "samples"), 'w' ).close()
            
            print( "copying", 'template', "into", fname )
            os.system( "cp -r %s/Results-template/* %s"%(PIPELINER_HOME, fname ) )
                
        except :
            showerror( "Initialization failed!", "Work directory data structure generation has failed." )
            return
        
        if self.make_symlinks() :
            showinfo( "Success", "The work directory has successfully initialized!")
        else :
            showerror( "Symlink failed", "" )
    
    def makejson_wrapper( self, *args, **kwargs ) :
        self.makejson(*args, **kwargs)
    
    