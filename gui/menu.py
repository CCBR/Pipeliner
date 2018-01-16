
from tkinter import Menu
from tkinter.messagebox import showinfo

def build_in_progress( msg="" ) :
	if not msg :
		msg = 'Sorry! We are building this functionality'
	showinfo( 'Build In Progress:', msg )

def about():
    info="""
    CCBR Pipeliner
    Version 3.0
    December, 2017.
    """
    showinfo("CCBR Pipeliner\nVersion 3.0",info)

def add_menubar( root ):
	menubar = Menu(root)

	#file menu
	filemenu = Menu(menubar, tearoff=0)
	# filemenu.add_command( label="Load project", command=build_in_progress )
	# filemenu.add_command( label="Save project", command=build_in_progress )
	filemenu.add_separator()
	filemenu.add_command( label="Exit", command=root.quit )

	menubar.add_cascade( label="File", menu=filemenu )

	#view menu
	viewmenu = Menu(menubar, tearoff=0)
	viewmenu.add_command( label="Progress", command=root.progress )
	viewmenu.add_command( label="Workflow", command=root.workflow )

	menubar.add_cascade( label="View", menu=viewmenu )

	#help menu
	helpmenu = Menu(menubar, tearoff=0)
	helpmenu.add_command( label="Help", command=build_in_progress )
	helpmenu.add_separator()
	helpmenu.add_command( label="About", command=about )

	menubar.add_cascade(label="Help", menu=helpmenu )
	
	root.config( menu=menubar )
