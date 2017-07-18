import sys,os,math,time
import json,re
import contextlib
import webbrowser
import time,threading
from io import StringIO

from pathlib import Path
from os import listdir, makedirs
from os.path import join, exists

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
        
        label = Label(eframe,text="Pipeline:")#,fg=textLightColor,bg=baseColor)
        label.grid(row=3,column=0,sticky=W,padx=10,pady=5)
        Pipelines=["InitialChIPseqQC", "ChIPseq" ]
        Pipeline = self.Pipeline = StringVar()
        Pipeline.set(Pipelines[0])
        
        om = OptionMenu(eframe, Pipeline, *Pipelines, command=self.option_controller)
        om.config()#bg = widgetBgColor,fg=widgetFgColor)
        om["menu"].config()#bg = widgetBgColor,fg=widgetFgColor)
        #om.pack(side=LEFT,padx=20,pady=5)
        om.grid(row=3,column=1,sticky=W,padx=20,pady=5)

        readtypes = ['Single', 'Paired']
        self.readtype = readtype = StringVar()
        readtype.set(readtypes[0])
        readtype_menu = OptionMenu(eframe, readtype, *readtypes)
        readtype_menu.grid(row=3, column=3, sticky=E, pady=5)
        readtype_label = Label(eframe, text="-end   ")
        readtype_label.grid( row=3, column=4, stick=W, pady=5)

        self.add_info(eframe)
        self.option_controller()
        self.peakinfo_fn = 'peakcall.tab'
        self.contrast_fn = 'contrast.tab'

    def option_controller( self, *args, **kwargs ) :
        PipelineFrame.option_controller( self )
        if self.Pipeline.get() == 'InitialChIPseqQC' :
            self.info.grid_forget()
        else :
       	    self.info.grid(row=7,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
    
    def add_info( self, parent ) :
        if not self.info :
            self.info = LabelFrame(parent, text='Peak Call Info')
            self.peakinfo_button = Button(self.info, 
                                            text="Set Peak Infomation", 
                                            command = self.popup_peakinfo )
            
            self.contrast_button = Button(self.info, 
                                            text="Set Groups to compare peaks", 
                                            command = self.popup_window_contrast )
            
            self.peakinfo_button.grid( row=5, column=5, padx=10, pady=5 )
            self.contrast_button.grid( row=5, column=6, padx=10, pady=5 )
            
    def popup_peakinfo( self ) :
        self.popup_window( "Peak Information",
                          join(self.workpath.get(),self.peakinfo_fn), 
                          'peakinfo' )
        
    def popup_window( self, title, filename, kind='textbox' ) :
        if kind == 'textbox' :
            PipelineFrame.popup_window(self, title, filename)
        elif kind == 'peakinfo' :
            self.popup_window_peakinfo(title, filename)
            
    
    def popup_window_contrast(self) :
        text = "Comparing ChIP-seq Peaks"
        contrast_fn = join(self.workpath.get(),self.contrast_fn)
        peakinfo_fn  = join(self.workpath.get(),self.peakinfo_fn)
        NA='N/A'
        
        try :
            groups = list(set([l.split('\t')[-1].strip() for l in open(peakinfo_fn) if l.split('\t')[-1].strip()]))
            if len(groups) < 2 :
                showwarning( 'WARNING!', "More than one groups need to be fefined in Peak Information File!" )
                print( 'groups:', groups )
                return
            
        except IOError :
            showerror( 'Error', "Did you set peak information?" )
            print('Error:', 'Cannot process peakcall.tab file:')
            return
        
        top = Toplevel()
        info = LabelFrame(top, text=text )#"Group Information")
        
        print( groups )
        contrast_vars = []
        contrast_menus = []
        n = 0
        groups.insert(0, NA)
        for i in range( int((len(groups)-1)*(len(groups)-2)/2) ):
            n = n + 1
            v1, v2 = StringVar(), StringVar()
            contrast_vars.append( [v1, v2] )
            o1, o2 = OptionMenu(info, v1, *groups), OptionMenu(info, v2, *groups)
            contrast_menus.append( [o1, o2] )

            v1.set(NA)
            v2.set(NA)
            vslabel = Label(info, text="  VS  ")

            o1.grid( row=n, column=0, padx=4, pady=1 )
            o2.grid( row=n, column=2, padx=4, pady=1 )
            vslabel.grid( row=n, column=1, padx=4, pady=1 )

        def savefunc() :
            info_text = StringIO()
            for v1, v2 in contrast_vars :
                v1 = v1.get() if v1.get() != NA else ""
                v2 = v2.get() if v2.get() != NA else ""
                
                if v1 and v2 :
                    pass
                elif v1 or v2 :
                    showerror( 'Error', "None or Both columns should be selected!" )
                    return
                else:
                    continue
                    
                print( v1, v2, file=info_text, sep="\t" )
            
            fp = open( contrast_fn, 'w' )
            fp.write( info_text.getvalue() )
            fp.close()
                              
        
        def loadfunc() :
            #self.readpaste( filename, info_text ) 
            for i, l in enumerate(open( contrast_fn )) :
                v1, v2 = l.split('\t')
                v1 = v1.strip()
                v2 = v2.strip()
                
                if v1 :
                    try : assert v1 in groups
                    except :
                        showwarning('WARNING', 'Group name is not in the selection list!' )
                        print( 'v1:',v1 ) 
                        print( 'group:', groups )
                        continue
                
                if v2 :
                    try: assert v2 in groups
                    except :
                        showwarning('WARNING', 'Group name is not in the selection list!' )
                        print( 'v2:',v2 ) 
                        print( 'group:', groups )
                        continue
                
                contrast_vars[i][0].set(v1)
                contrast_vars[i][1].set(v2)
                
        
        def clearfunc() :
            for v1, v2 in contrast_vars :
                v1.set(NA)
                v2.set(NA)
                
        
        info_clear_button = Button(top, 
                                  text="Clear", 
                                  command = clearfunc )
        info_save_button = Button(top, 
                                  text="Save", 
                                  command = savefunc )
        info_load_button = Button(top,
                                  text="Load",
                                  command = loadfunc )

        #self.pairs_load_button.pack( side=BOTTOM, padx=5, pady=5 )
        #self.pairs_save_button.pack( side=BOTTOM, padx=5, pady=5 )

        info_clear_button.grid( row=5, column=3, padx=10, pady=5 )
        info_load_button.grid( row=5, column=4, padx=10, pady=5 )
        info_save_button.grid( row=5, column=5, padx=10, pady=5 )


        info.grid(row=7,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
        top.focus_force()
        
            
    def popup_window_peakinfo( self, text, filename ) :
        
        NA = 'N/A'
        selections = [fn.split('.R1.fastq.gz')[0] 
                     for fn in self.datafiles if fn.endswith( ".R1.fastq.gz" )]
        
        ################################
        #check availablity of the files before moving on!
        if not selections :
            showerror( "No FASTQ files available matching the pattern *.R1.fastq.gz" )
            return
        ################################
        
        selections.insert(0, NA) #Adding N/A for the groups where no input is available
        groups = ['Grp%d'%i for i in range(1, len(selections)+1)]
        
        ##checking for debugging purpose
        print( selections )
        assert len(selections) == len(set(selections))
        ##
        
  
        top = Toplevel()
        info = LabelFrame(top, text=text )#"Group Information")
        chip_vars = [ StringVar() for s in selections[1:] ]
        input_vars = [ StringVar() for s in selections[1:] ]
        group_vars = [ StringVar() for s in selections[1:] ]
        
        chip_menus = [OptionMenu(info,var,*selections) for var in chip_vars]
        input_menus = [OptionMenu(info,var,*selections) for var in input_vars]
        group_menus = [OptionMenu(info,var,*groups) for var in group_vars]
        group_entries = [Entry(info, bd=2, width=8, textvariable=var) for var in group_vars]
        
        chiplabel = Label(info, text= "ChIP Names")
        inputlabel = Label(info, text="Input Names")
        grouplabel = Label(info, text="Group Names")
        
        chiplabel.grid( row = 0, column = 1, padx=4, pady=1)
        inputlabel.grid( row = 0, column = 2, padx=4, pady=1)
        grouplabel.grid( row = 0, column = 3, padx=4, pady=1 )
        
        
        for i, (chvar, invar) in enumerate(zip(chip_vars, input_vars)) :
            chvar.set(selections[0])
            invar.set(selections[0])
            
            chip_menus[i].grid( row = i+1, column = 1, padx=4, pady=1 )
            input_menus[i].grid( row = i+1, column = 2, padx=4, pady=1 )
            group_entries[i].grid( row = i+1, column = 3, padx=4, pady=1 )
            group_menus[i].grid( row = i+1, column = 4, padx=4, pady=1 )
        
        
        def savefunc() :
            info_text = StringIO()
            for v1, v2, v3 in zip( chip_vars, input_vars, group_vars ) :
                v1 = v1.get().strip() if v1.get().strip() != NA else ""
                
                if not v1 :
                    continue
                    
                v2 = v2.get().strip() if v2.get().strip() != NA else ""
                v3 = v3.get().strip() if v3.get().strip() != NA else ""
                
                if not v3 :
                    showerror( "Error", "Missing Replicate group name detected.\nReplicate group names should be given!" )
                    print( "Error", "Missing Replicate group name detected.\nReplicate group names should be given!" )
                    return
                    
                print( v1, v2, v3, file=info_text, sep="\t" )
            
            fp = open( filename, 'w' )
            fp.write( info_text.getvalue() )
            fp.close()
        
        def loadfunc() :
            if not exists(filename) :
                print( filename, 'does not exists!' )
                return
            
            for i, l in enumerate(open( filename )) :
                v1, v2, v3 = l.split('\t')
                
                if v1 :
                    try : assert v1 in selections
                    except :
                        showwarning('WARNING', 'ChIP name is not in the selection list!' )
                        print( 'v1:',v1 ) 
                        print( 'selection:', selection )
                        continue
                
                if v2 :
                    try: assert v2 in selections
                    except :
                        showwarning('WARNING', 'Input name is not in the selection list!' )
                        print( 'v2:',v2 ) 
                        print( 'selection:', selection )
                        return
                    
                chip_vars[i].set(v1.strip())
                input_vars[i].set(v2.strip())
                group_vars[i].set(v3.strip())
            
        def clearfunc() :
            for i, (chvar, invar, grvar) in enumerate(zip(chip_vars, input_vars, group_vars)) :
                chvar.set(selections[0])
                invar.set(selections[0])
                grvar.set('')
        
        info_clear_button = Button(top, 
                                  text="Clear", 
                                  command = clearfunc )
        
        info_save_button = Button(top, 
                                  text="Save", 
                                  command = savefunc )
        
        info_load_button = Button(top,
                                  text="Load",
                                  command = loadfunc )

        #self.pairs_load_button.pack( side=BOTTOM, padx=5, pady=5 )
        #self.pairs_save_button.pack( side=BOTTOM, padx=5, pady=5 )

        info_clear_button.grid( row=5, column=3, padx=10, pady=5 )
        info_load_button.grid( row=5, column=4, padx=10, pady=5 )
        info_save_button.grid( row=5, column=5, padx=10, pady=5 )

        info.grid(row=7,column=0, columnspan=6, sticky=W, padx=20, pady=10 )
        top.focus_force()
        
        
    def makejson(self, *args):
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
    ##------   
        D=dict() 
        FT = filetype#.get()
    #    p = Popen("ls "+workpath.get()+"/*."+FT, shell=True, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, close_fds=True)
        p = Popen("find "+self.workpath.get()+" -maxdepth 1 -type l -printf '%f\n' ", shell=True, stdin=PIPE, stdout=PIPE, stderr=DEVNULL, close_fds=True)
        a = p.stdout.read().decode(encoding='UTF-8').split("\n")

        RG=dict()   
        b=a.pop()
        #tkinter.messagebox.showerror("",a)
        #if freezeunits.get()=="no":
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
            RG[key]={
                'rgsm':key,
                'rglb':'na',
                'rgpu':'na',
                'rgpl':'ILLUMINA',
                'rgcn':'na'
            }
        units=D
        UnitsBak=D

        try:
            F=open(self.workpath.get()+"/rg.tab","r")
            f=F.read().splitlines()
            F.close()
            for theLine in f:
                if not re.match("^ID",theLine):
                    (rgid,rgsm,rglb,rgpl,rgpu,rgcn)=theLine.split("\t")
                    RG[rgid]={
                        'rgsm': rgsm,
                        'rglb': rglb,
                        'rgpu': rgpu,
                        'rgpl': rgpl,
                        'rgcn': rgcn
                    }
        except:
            pass
        
        RGbak=RG
    #     else:
    #         units=UnitsBak
    #         RG=RGbak
    # 
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
        
        ################################
        #genome size for macs
        ################################
        gsize = 'hs' if self.genome[:2] == 'hg' else self.genome[:2] #hs or mm
        ################################
        peaks = {}
        contrast = []
        groups = {}
        if self.Pipeline.get() == 'ChIPseq' :
            ################################
            #peak calling info!
            ################################
            peaks['chips'] = []
            peaks['inputs'] = {}
            groups = {g:[] for g in list(set([l.split('\t')[-1].strip() for l in open(join(self.workpath.get(),self.peakinfo_fn)) if l.split('\t')[-1].strip()]))}
            
            for l in open(join(self.workpath.get(),self.peakinfo_fn)) :
                line = [el.strip() for el in l.split('\t')]
                
                peaks['chips'].append(line[0])
                if len(line) == 1 :
                    peaks['inputs'][line[0]] = ''
                elif len(line) > 1 :
                    peaks['inputs'][line[0]] = line[1]

                groups[line[2]].append( line[0] )
            ################################
            #Contrast info processing!
            ################################
            for l in open( join(self.workpath.get(),self.contrast_fn)) :
                c1, c2 = l.split()
                contrast.append( [c1, c2] )
    
        gi = self.global_info
        PD={
            'project': {
                'pfamily': gi.pfamily.get(),
                'units': units, 
                'samples': samples, 
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
                'batchsize': batchsize, 
                "smparams": smparams, 
                "rgid": RG, 
                "cluster": "cluster_medium.json", 
                "description": gi.description.get('1.0',END), 
                "technique" : gi.technique.get(), 
                'readtype': self.readtype.get(),
                'pairs': {'na':'na'},
                'peaks': peaks,
                'gsize': gsize,
                'groups': groups,
                'contrast': contrast,
             }
        } 
        
        J=json.dumps(PD, sort_keys = True, indent = 4, ensure_ascii=True)
        gi.jsonconf.delete("1.0", END)    
        gi.jsonconf.insert(INSERT, J)
        self.saveproject(gi.jsonconf.get("1.0",END))
    
   