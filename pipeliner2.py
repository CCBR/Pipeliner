#!/usr/bin/env python3

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


#get the environmental variables
from os import environ as ENV
from os import getcwd, listdir, makedirs, symlink
from os.path import join, exists, dirname, abspath

from gui.utils import PIPELINER_HOME, PIPELINER_CONF, USER_HOME
###############################################
#Legacy Utility Library 
#Mostly inherited from David's design
#Probably need to be checked for tidiness
###############################################
from gui.utils import CustomNotebook, Cfg, show, copy_dir
from gui.utils import ftypes, filetypes, filetype, efiletype, batchsize, customRules
from gui.utils import dryrunFgColor, dryrunBgColor
###############################################

PWD = getcwd()
VERBOSE=1

assert PIPELINER_HOME
assert PIPELINER_CONF
assert PWD

print( "USER_HOME:", USER_HOME )
print( "PIPELINER_HOME:", PIPELINER_HOME )
print( "PIPELINER_CONF:", PIPELINER_CONF )
print( "PWD:", PWD )

from gui.frame import PipelineFrame
from gui.rnaseq import RNASeqFrame
from gui.genomeseq import GenomeSeqFrame
from gui.mirseq import MiRSeqFrame
from gui.epigenomeseq import ChIPSeqFrame
from gui.exomeseq import ExomeSeqFrame

class PipelinerGUI(Tk):
    def __init__(self,parent=None):
        Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.grid()

        self.init_project_frame()
            
        self.grid_columnconfigure(0,weight=1)
        self.resizable(True,True)
        
        self.geometry("650x760")
        self.update()

    def init_project_frame( self ) :
        projframe = Frame(self)
        
        ######################
        #The Project Entry Form
        #
        #Current design inherited from others
        #is very bad. It should have been made as
        #json or dictionary as is in this point
        #then the subsequent coding will be much
        #eaiser.
        #BK
        ######################
        BS=StringVar()

        self.eprojectid = eprojectid= StringVar()
        self.eorganism = eorganism=StringVar()
        self.eanalyst = eanalyst=StringVar()
        self.epoc = epoc=StringVar()
        self.epi = epi=StringVar()
        self.euser = euser=StringVar()
        self.eflowcell = eflowcell=StringVar()
        self.eplatform = eplatform=StringVar()
        eplatform.set("Illumina")
        self.technique = technique=StringVar()

        editframe = Frame( self )
        editframe.grid( column=0, row=3, columnspan=3 )
        
        projpanel1 = LabelFrame(editframe,text="Project Information") 
        #,fg=textLightColor,bg=baseColor)
        projpanel1.pack( side = TOP, fill=X, padx=10, pady=5, expand=YES )
        
        pipeline_panel = LabelFrame(editframe, text="Global Settings")
        pipeline_panel.pack( side=TOP, fill=X, padx=10, pady=5, expand=YES )
        
        l=Label( pipeline_panel, text="Genome:" )
        l.grid(row=1,column=1,sticky=W,padx=0,pady=5)
        l = Label( pipeline_panel, text="Pipeline Family:" )
        l.grid(row=1,column=3,sticky=W,padx=0,pady=5)
        
        annotations=['hg19','mm10']
        self.annotation = annotation = StringVar()
        annotation.set(annotations[0])
        #annotation.trace('w', lambda *_ :settargets(annotation) )
        
        om = OptionMenu(pipeline_panel, annotation, *annotations)
        #, command=lambda _:makejson("refsets"))
        om.config() #bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config() #bg = widgetBgColor,fg=widgetFgColor)
        om.grid(row=1,column=2,sticky=W,padx=10,pady=10)
        
        pfamilys = ['exomeseq', 'rnaseq', 'genomeseq', 'mirseq', 'chipseq', 'custom']
        self.pfamily = pfamily = StringVar()
        pfamily.set(pfamilys[0])
        om = OptionMenu(pipeline_panel, pfamily, *pfamilys)
        #, command=lambda _:makejson(pfamily.get()))
        om.config() #bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config() #bg = widgetBgColor,fg=widgetFgColor)
        om.grid(row=1,column=4,sticky=W,padx=10,pady=5)
        add_button = Button( pipeline_panel, 
                            text="Set a pipeline", 
                            command=self.add_pipeline
                           )
        add_button.grid(row=1, column=5, sticky=W, padx=10, pady=5)
        
        self.notebook = notebook = Notebook( editframe )
        self.pipelineframe = None #the pipeline frame in the notebook frame!
        
        projpanel2 = Frame(notebook) #,fg=textLightColor,bg=baseColor)
        projpanel2.pack( side = LEFT, fill=BOTH, padx=10, pady=5, expand=YES )
        notebook.add( projpanel2, text="Project Description" )
        notebook.pack( side=LEFT, fill=BOTH, padx=10, pady=5, expand=YES )

        Dscrollbar = Scrollbar(projpanel2)
        Dscrollbar.grid(row=2,column=4,rowspan=40)
        
        self.description =description= Text(projpanel2,
                           width=75,
                           height=28,
                           #bg=commentBgColor,
                           #fg=commentFgColor,
                           font=("nimbus mono bold","10"),
                           yscrollcommand = Dscrollbar.set)  
        
        description.delete("1.0", END)
        description.insert(INSERT, "Enter CCBR Project Description and Notes here.")
        #description.bind('<FocusOut>',lambda _:makejson("none"))
        description.grid(row=2,column=3,sticky="e",padx=10,pady=5)
        Dscrollbar['command']=description.yview
        
        L=Label(projpanel1, text="Project Id", anchor="ne" )
        #,bg=baseColor,fg=textLightColor)
        E = Entry(projpanel1, 
                  bd =2, 
                  width=20, 
                  #bg=entryBgColor, 
                  #fg=entryFgColor, 
                  textvariable=eprojectid 
                 )
        #L.pack(pady=5,padx=5)
        #E.pack(pady=1,padx=5)
        L.grid(row=0, column=0, sticky=W, padx=10, pady=5)
        E.grid(row=0, column=2, sticky=W, padx=10, pady=5)
        eprojectid.set("project")
        #eprojectid.trace('w', makejson)
        L2=Label(projpanel1,
                 text = "(Examples: CCBR-nnn,Labname or short project name)",
                 anchor="ne"
                )#,bg="firebrick",fg="white")
        L2.grid(row=0, column=3, padx=10, pady=5, sticky="w")

        L=Label( projpanel1,
                text="Email address",
                anchor="ne")#,bg=baseColor,fg=textLightColor)
        
        E = Entry( projpanel1, 
                  bd =2, 
                  width=20, 
                  #bg=entryBgColor, 
                  #fg=entryFgColor,
                  textvariable=euser)
        
        L3 = Label(projpanel1,
                   text ="(Mandatory field: must use @nih.gov email address)",
                   anchor="ne")
        L3.grid(row=1,column=3,padx=10,pady=5,sticky="w")
        L.grid(row=1,column=0,sticky=W,padx=10,pady=5)
        E.grid(row=1,column=2,sticky=W,padx=10,pady=5)

        #euser.trace('w', makejson)

        L=Label(projpanel1,text="Flow Cell ID",anchor="ne")
        E = Entry(projpanel1, bd =2, width=20, textvariable=eflowcell )#, bg=entryBgColor)
        L4=Label( projpanel1,
                 text="(Examples: FlowCellID, Labname, date or short project name)",
                 anchor="ne" 
                )#,bg="firebrick",fg="white")
        L4.grid(row=2,column=3,padx=10,pady=5,sticky="w")
        L.grid(row=2,column=0,sticky=W,padx=10,pady=5)
        E.grid(row=2,column=2,sticky=W,padx=10,pady=5)
        eflowcell.set("stats")
        #eflowcell.trace('w', makejson)

        Projscrollbar = Scrollbar(projframe)
        Projscrollbar.grid(row=2,column=4,rowspan=30 )
        self.jsonconf = jsonconf = Text( projframe,
                        width=80,
                        height=30,
                        #bg=projectBgColor,
                        #fg=projectFgColor,
                        font=("nimbus mono bold","11"),
                        yscrollcommand = Projscrollbar.set)
        
        Projscrollbar['command']=jsonconf.yview
        jsonconf.grid(row=3,column=0,sticky="e",padx=10,pady=5)
        #jsonconf.bind('<FocusOut>',lambda _:makejson("none"))

        
    def init_global_options( self ) :
        pass
    
        
    def add_pipeline( self, *args ) :
        print( "set_pipeline", *args )
        print( self.pfamily )

        if self.pipelineframe :
            self.notebook.forget( self.pipelineframe )
            
        if self.pfamily.get() == 'exomeseq' :
            print( 'exomeseq' )
            self.pipelineframe = ExomeSeqFrame( self.notebook, 
                                               'ExomeSeq', #self.pfamily.get(), 
                                               self.annotation, global_info=self )

        elif self.pfamily.get() == 'genomeseq' :
            print( 'genomeseq' )
            self.pipelineframe = GenomeSeqFrame( self.notebook, 
                                               'GenomeSeq', #self.pfamily.get(), 
                                               self.annotation, global_info=self )

       	elif self.pfamily.get() == 'mirseq' :
            print( 'mirseq' )
            self.pipelineframe = MiRSeqFrame( self.notebook, 
                                               'miR-Seq', #self.pfamily.get(), 
                                               self.annotation, global_info=self )
     
        elif self.pfamily.get() == 'chipseq' :
            print( 'chipseq' )
            self.pipelineframe = ChIPSeqFrame( self.notebook, 
                                               'ChIPseq',#self.pfamily.get(), 
                                               self.annotation, global_info=self )
        elif self.pfamily.get() == 'rnaseq' :
            print( 'rnaseq' )
            self.pipelineframe = RNASeqFrame( self.notebook, 
                                               'RNAseq',#self.pfamily.get(), 
                                               self.annotation, global_info=self )
        else :
            pass

        self.notebook.select( self.pipelineframe )

if __name__ == "__main__":
    gui = PipelinerGUI()
    gui.title('CCBR Pipeliner')
    gui.mainloop()

    