import tkinter as tk
from tkinter import VERTICAL, HORIZONTAL, RIGHT, BOTTOM, INSERT, X, Y, BOTH, Text
from tkinter import ttk
from tkinter.ttk import Notebook


import os
import shutil

class CustomNotebook(Notebook):
    """A ttk Notebook with close buttons on each tab"""
    '''Donwloaded from stackoverflow and slightly modified'''

    __initialized = False

    def __init__(self, *args, **kwargs):
        if not self.__initialized:
            self.__initialize_custom_style()
            self.__inititialized = True

        kwargs["style"] = "CustomNotebook"
        Notebook.__init__(self, *args, **kwargs)

        self._active = None

        self.bind("<ButtonPress-1>", self.on_close_press, True)
        self.bind("<ButtonRelease-1>", self.on_close_release)

    def on_close_press(self, event):
        """Called when the button is pressed over the close button"""

        element = self.identify(event.x, event.y)

        if "close" in element:
            index = self.index("@%d,%d" % (event.x, event.y))
            self.state(['pressed'])
            self._active = index

    def on_close_release(self, event):
        """Called when the button is released over the close button"""
        if not self.instate(['pressed']):
            return

        element =  self.identify(event.x, event.y)
        index = self.index("@%d,%d" % (event.x, event.y))

        if "close" in element and self._active == index:
            self.forget(index)
            self.event_generate("<<NotebookTabClosed>>")

        self.state(["!pressed"])
        self._active = None

    def __initialize_custom_style(self):
        style = ttk.Style()
        self.images = (
            tk.PhotoImage("img_close", data='''
                R0lGODlhCAAIAMIBAAAAADs7O4+Pj9nZ2Ts7Ozs7Ozs7Ozs7OyH+EUNyZWF0ZWQg
                d2l0aCBHSU1QACH5BAEKAAQALAAAAAAIAAgAAAMVGDBEA0qNJyGw7AmxmuaZhWEU
                5kEJADs=
                '''),
            tk.PhotoImage("img_closeactive", data='''
                R0lGODlhCAAIAMIEAAAAAP/SAP/bNNnZ2cbGxsbGxsbGxsbGxiH5BAEKAAQALAAA
                AAAIAAgAAAMVGDBEA0qNJyGw7AmxmuaZhWEU5kEJADs=
                '''),
            tk.PhotoImage("img_closepressed", data='''
                R0lGODlhCAAIAMIEAAAAAOUqKv9mZtnZ2Ts7Ozs7Ozs7Ozs7OyH+EUNyZWF0ZWQg
                d2l0aCBHSU1QACH5BAEKAAQALAAAAAAIAAgAAAMVGDBEA0qNJyGw7AmxmuaZhWEU
                5kEJADs=
            ''')
        )

        style.element_create("close", "image", "img_close",
                            ("active", "pressed", "!disabled", "img_closepressed"),
                            ("active", "!disabled", "img_closeactive"), border=8, sticky='')
        style.layout("CustomNotebook", [("CustomNotebook.client", {"sticky": "nswe"})])
        style.layout("CustomNotebook.Tab", [
            ("CustomNotebook.tab", {
                "sticky": "nswe", 
                "children": [
                    ("CustomNotebook.padding", {
                        "side": "top", 
                        "sticky": "nswe",
                        "children": [
                            ("CustomNotebook.focus", {
                                "side": "top", 
                                "sticky": "nswe",
                                "children": [
                                    ("CustomNotebook.label", {"side": "left", "sticky": ''}),
                                    ("CustomNotebook.close", {"side": "left", "sticky": ''}),
                                ]
                        })
                    ]
                })
            ]
        })
    ])

        
def show(x,y,fg,bg,height,width):
    sm=tk.Toplevel()
    sm.title(y)
    smyscrollbar = ttk.Scrollbar(sm, orient=VERTICAL)
    smyscrollbar.pack( side = RIGHT, fill=Y )
    smxscrollbar = ttk.Scrollbar(sm, orient=HORIZONTAL)
    smxscrollbar.pack( side = BOTTOM, fill=X )
    smout = Text(sm,width=width,height=height,bg=bg,fg=fg,yscrollcommand = smyscrollbar.set, xscrollcommand = smxscrollbar.set)
#    smout = Text(sm,width=100,height=100,bg="gray30",fg="yellow",yscrollcommand = smyscrollbar.set, xscrollcommand = smxscrollbar.set)
    smout.insert(INSERT, x)
    smout.pack(expand=1,fill=BOTH,side=RIGHT)
    smyscrollbar.config(command=smout.yview)
    smxscrollbar.config(command=smout.xview)
    
def copy_dir(src, dst, *, follow_sym=True):
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    if os.path.isdir(src):
        shutil.copyfile(src, dst, follow_symlinks=follow_sym)
        shutil.copystat(src, dst, follow_symlinks=follow_sym)
    return dst
        
def load_config(cfile):
    Q=open(cfile,"r")
    Cfg=eval(Q.read())
    Q.close()
    return Cfg