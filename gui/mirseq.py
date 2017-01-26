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
from os.path import exists, join

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


class MiRSeqFrame( PipelineFrame ) :
    def __init__(self, pipepanel, pipeline_name, *args, **kwargs) :
        PipelineFrame.__init__(self, pipepanel, pipeline_name, *args, **kwargs)
        self.info = None
        
        eframe = self.eframe = LabelFrame(self,text="Options") 
        #,fg=textLightColor,bg=baseColor)
        eframe.grid( row=5, column=1, sticky=W, columnspan=7, padx=10, pady=5 )
        
        label = Label(eframe,text="Pipeline")#,fg=textLightColor,bg=baseColor)
        label.grid(row=3,column=0,sticky=W,padx=10,pady=5)
        Pipelines=["CAPmirseq-plus"]
        Pipeline = self.Pipeline = StringVar()
        Pipeline.set(Pipelines[0])
        
        om = OptionMenu(eframe, Pipeline, *Pipelines, command=self.option_controller)
        om.config()#bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=10,pady=5)
     
        self.add_info( eframe )
        
    
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
            

        self.info.grid(row=10,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
   
  
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
 