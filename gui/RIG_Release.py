import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox

import numpy as np

#Define lists of options
RootMobilities = ["Cartesian", "Free", "Weld"]

#define variables:

default_text = "<None>"

G0 = 0 
T0 = 0  #These variables make sure we don't reinstantiate any object and lose the values the user assigned them
W0 = 0
X0 = 0
DidWorlds = 0 

def add2List(element, list0):
	list0.append(element)
#	print (list0)

def changeValue(variable, value):
	variable.set(value)

def NotYet(variable, value):
	messagebox.showwarning("Coming Soon", "This option isn't implemented yet.")
	variable.set(value)

def getParams(field):
	field.set(filedialog.askopenfilename(title = "Select Parameter file", filetypes = (("prmtop files","*.prmtop"),("all files","*.*"))))

def getCoords(field):
	field.set(filedialog.askopenfilename(title = "Select Parameter file",filetypes = (("inpcrd files","*.inpcrd"),("rst7 files","*.rst7"),("all files","*.*"))))

def getFlex(field):
	field.set(filedialog.askopenfilename(title = "Select Flexibility file",filetypes = (("flex files","*.flex"),("all files","*.*"))))

def getRB(field):
	field.set(filedialog.askopenfilename(title = "Select Rigid Bodies file",filetypes = (("rb files","*.rb"),("all files","*.*"))))

def disableEntry(field):
	field.configure(state = tk.DISABLED)

def enableEntry(field):
	field.configure(state = tk.NORMAL)

def add(variable, window, value):

	global G0, T0, W0, X0				#this is the most efficient way I've found of doing this. I can't find how to pass a global variable as a parameter to a function, and have that function change it.	
	if (variable == "G0"):
		TheGUI.setParams(value, "G0")
		TheGUI.updateStatus("G0")
	if (variable == "T0"):
		TheGUI.setParams(value, "T0")
		TheGUI.updateStatus("T0")
	if (variable == "W0"):
		TheGUI.setWorldParams(TheGUI.eExpvariable.get(), TheGUI.eW0variable.get())		#world uses a special function to set its params
		TheGUI.updateStatus("W0")
	if (variable == "X0"):
		TheGUI.setParams(value, "X0")
		TheGUI.updateStatus("X0")


def WorldParamsClose(window):

	global DidWorlds
	DidWorlds = 1
	window.destroy()

	
class MainWindow:

	def __init__(self, mainWin):
		self.mainWin = mainWin
		self.mainWin.title("Robosample Input Generator")
		self.eExpNumber=tk.IntVar()

		lNumberOfExperiments = tk.Label(self.mainWin, text="How many experiments do you wish to run?")
		lNumberOfExperiments.grid(row=0, column=0)
		eNumberOfExperiments = tk.Entry(self.mainWin, textvariable=self.eExpNumber)
		eNumberOfExperiments.grid(row=0, column=2)
		b1 = tk.Button(self.mainWin, text="Start", command=self.openParamsWindow, width=35)
		b1.grid(row=1,column=0,columnspan=3)

	def openParamsWindow(self):			

		#This is the main parameters window	
	
		if (self.eExpNumber.get() <= 0):
			messagebox.showwarning("Oops!", "Please select a positive, non-zero number of experiments!.")
			return
		self.ParamsWin = tk.Toplevel()
		self.ParamsWin.title("Parameters")
		experiments=[]

		for i in range(self.eExpNumber.get()):
			experiments.append(i)

		#we have to add a drop-down menu so the user can select what experiment he wants to choose parameters for
	
		lExp = tk.Label(self.ParamsWin, text="Experiment: ")
		lExp.grid(row=0,column=0)
		self.eExpvariable = tk.IntVar()
		self.eExpvariable.set(experiments[0])
		eExp = tk.OptionMenu(self.ParamsWin, self.eExpvariable,*experiments)
		eExp.grid(row=0,column=1,columnspan=1)
		

		bGen = tk.Button(self.ParamsWin, text = "General", command = self.openGeneralWindow)
		bThermo = tk.Button(self.ParamsWin, text= "Thermodynamics", command = self.openThermodynWindow)
		bWorld = tk.Button(self.ParamsWin, text="World Params", command = self.openWorldWindow)
		bOut = tk.Button(self.ParamsWin, text="Output Params", command = self.openOutputWindow)

		self.eFinishPath=tk.StringVar()
		self.eInputNumber=tk.IntVar()
		self.eInputNumber.set(1)
		bFinish = tk.Label(self.ParamsWin, text="Input Name Prefix")
		eFinish = tk.Entry(self.ParamsWin, textvariable=self.eFinishPath)
		bNumberOfInputs = tk.Label(self.ParamsWin, text="Number of inputs:")
		eNumberOfInputs = tk.Entry(self.ParamsWin, textvariable=self.eInputNumber)

		bFinishButton = tk.Button(self.ParamsWin, text="Write experiment to file!", command = lambda: self.writeAllToFile(self.eFinishPath.get(), self.eInputNumber.get(), self.eExpvariable.get()))


		bGen.grid(row=1, column=0)
		bThermo.grid(row=1, column=1)
		bWorld.grid(row=1, column=2)
		bOut.grid(row=1, column=3)
		l = tk.Label(self.ParamsWin, text = "")
		l.grid(row=3,column=0)
		bNumberOfInputs.grid(row=4, column=0, columnspan=2)
		eNumberOfInputs.grid(row=4, column=2, columnspan=2)
		bFinish.grid(row=5, column=0, columnspan=2)
		eFinish.grid(row=5, column=2, columnspan=2)
		bFinishButton.grid(row=6, column=1, columnspan=2)
		

		#This is the window that keeps track of all the parameters set thus far
		
		self.StatusWindow = tk.Toplevel()
		self.StatusWindow.title("Parameters Overview")
		DefaultText = tk.StringVar()
		DefaultText.set("<NOT SET>")

		#We add labelling that will remain the same throughout
			
		StatusLabelsDictG = {"StatusLabelG0":"Project Name", "StatusLabelG1":"Parameters File", "StatusLabelG2":"Coordinates File", \
										"StatusLabelG3":"Number of Worlds", "StatusLabelG4":"Number of Rounds", "StatusLabelG6":"Use OpenMM", \
										"StatusLabelG7":"Number of Threads", "StatusLabelG8":"Reproducible"}
		StatusLabelsG = ["StatusLabelG0", "StatusLabelG1", "StatusLabelG2", "StatusLabelG3", "StatusLabelG4", "StatusLabelG6", \
											"StatusLabelG7", "StatusLabelG8"]

		StatusCounter = 1
		for LabelObj in StatusLabelsG:
			LabelObj = tk.Label(self.StatusWindow, text=StatusLabelsDictG[LabelObj])
			LabelObj.grid(row=StatusCounter, column=0)
			StatusCounter = StatusCounter + 1

		StatusLabelSpace = tk.Label(self.StatusWindow, text ="")
		StatusLabelSpace.grid(row=StatusCounter,column=0)
		StatusCounter = StatusCounter + 1

		StatusLabelsDictT = {"StatusLabelT0":"Thermostat", "StatusLabelT1":"Initial Temp", "StatusLabelT2":"Final Temp", \
										"StatusLabelT3":"Boost Temp", "StatusLabelT4":"Forcefield Scale", "StatusLabelT5":"GBSA Scale Factor"}
		StatusLabelsT = ["StatusLabelT0", "StatusLabelT1", "StatusLabelT2", "StatusLabelT3", "StatusLabelT4", "StatusLabelT5"]
		
		for LabelObj in StatusLabelsT:
			LabelObj = tk.Label(self.StatusWindow, text=StatusLabelsDictT[LabelObj])
			LabelObj.grid(row=StatusCounter, column=0)
			StatusCounter = StatusCounter + 1

		StatusLabelSpace = tk.Label(self.StatusWindow, text ="")
		StatusLabelSpace.grid(row=StatusCounter,column=0)
		StatusCounter = StatusCounter + 1

		#We don't put the worlds in this window, jump straight to Output params
		
		StatusLabelsDictX = {"StatusLabelX0":"Data Print Frequency", "StatusLabelX1":"Coordinates Print Frequency", "StatusLabelX2":"Print Geometrical Features",	"StatusLabelX3":"Geom. Features to Print", "StatusLabelX4":"Use Visual?"}
		StatusLabelsX = ["StatusLabelX0", "StatusLabelX1", "StatusLabelX2", "StatusLabelX3", "StatusLabelX4"]
		
		for LabelObj in StatusLabelsX:
			LabelObj = tk.Label(self.StatusWindow, text=StatusLabelsDictX[LabelObj])
			LabelObj.grid(row=StatusCounter, column=0)
			StatusCounter = StatusCounter + 1

		StatusLabelSpace = tk.Label(self.StatusWindow, text ="")
		StatusLabelSpace.grid(row=StatusCounter,column=0)


		#column headers
		
		for i in range(self.eExpNumber.get()):
			experNumber = tk.Label(self.StatusWindow, text = i)
			experNumber.grid(row=0, column=i+1)
	
		
		#now we add the fields where the params go
		
		for i in range(self.eExpNumber.get()):
			for j in range(len(StatusLabelsDictG)):
		#General Fields
				GeneralField = tk.Entry(self.StatusWindow, textvariable = DefaultText, state = "readonly", justify=tk.CENTER, fg="red")
				GeneralField.grid(row=j+1, column=i+1)
		#Thermodynamic Fields
			space = tk.Label(self.StatusWindow, text="")
			space.grid(row=9,column=i+1)
			for j in range(len(StatusLabelsDictT)):
				ThermoField  = tk.Entry(self.StatusWindow, textvariable = DefaultText, state = "readonly", justify=tk.CENTER, fg="red")
				ThermoField.grid(row=j+10, column=i+1)
		#Output Fields
			ace = tk.Label(self.StatusWindow, text="")
			space.grid(row=16,column=i+1)
			for j in range(len(StatusLabelsDictX)):
				OutputField  = tk.Entry(self.StatusWindow, textvariable = DefaultText, state = "readonly", justify=tk.CENTER, fg="red")
				OutputField.grid(row=j+17, column=i+1)

	def openGeneralWindow(self):
		self.GeneralWin = tk.Toplevel()
		self.GeneralWin.title("")

		GeneralTitle = tk.Label(self.GeneralWin, text = "General Parameters", bd=2, relief=tk.SOLID, pady=5, padx=5)
		GeneralTitle.grid(row=0, column=0, columnspan=3)	

		space=tk.Label(self.GeneralWin, text=" ")
		space.grid(row=1,column=0,columnspan=3)

		#Make the labels
		StatusLabelsDictG = {"lG0":"Project Name:", "lG1":"Parameters", "lG2":"Coordinates",	"lG3":"Number of Worlds:", "lG4":"Number of Rounds:", "lG6":"Use OpenMM(GPU)", "lG7":"Number of Threads:", "lG8":"Reproducible Simulation"}
		StatusLabelslG = ["lG0", "lG1", "lG2", "lG3", "lG4", "lG6", "lG7", "lG8"]
	
		GeneralCounter = 2	
		for LabelObj in StatusLabelslG:
			LabelObj = tk.Label(self.GeneralWin, text=StatusLabelsDictG[LabelObj])
			LabelObj.grid(row=GeneralCounter, column=0)
			GeneralCounter = GeneralCounter + 1

		#make the fields where you actually add things

		global G0
		if (G0 == 1):
			pass
		else:
			#since we're working with multiple experiments, instead of variables, we use 1-D arrays of variables. Luckily for us, whereas the number of worlds can change, the number of experiments can not, so we don't have to program for that

			self.eG0variableArray = np.full(self.eExpNumber.get(),fill_value=default_text, dtype='object')	
			self.eG1variableArray = np.full(self.eExpNumber.get(),fill_value=default_text, dtype='object')	
			self.eG2variableArray = np.full(self.eExpNumber.get(),fill_value=default_text, dtype='object')	
			self.eG3variableArray = np.full(self.eExpNumber.get(),fill_value=-1, dtype='object')	
			self.eG4variableArray = np.full(self.eExpNumber.get(),fill_value=-1, dtype='object')	
			self.eG6variableArray = np.full(self.eExpNumber.get(),fill_value=2 , dtype='object')	
			self.eG7variableArray = np.full(self.eExpNumber.get(),fill_value=-1, dtype='object')	
			self.eG8variableArray = np.full(self.eExpNumber.get(),fill_value=2 , dtype='object')	
		
			self.eG0variable = tk.StringVar()	
			self.eG0variable.set("Experiment " + str(self.eExpvariable.get()))
			self.eG1variable = tk.StringVar()
			self.eG1variable.set(default_text)
			self.eG2variable = tk.StringVar()
			self.eG2variable.set(default_text)
			self.eG3variable = tk.IntVar()
			self.eG3variable.set(-1)
			self.eG4variable = tk.IntVar()
			self.eG4variable.set(-1)
			self.eG6variable = tk.BooleanVar()
			self.eG6variable.set(2) 
			self.eG7variable = tk.IntVar()
			self.eG7variable.set(0)
			self.eG8variable = tk.BooleanVar()
			self.eG8variable.set(2)
			G0 = 1 #so we know we set the variables in General and so we don't reinstantiate all fields every time

		self.eG0 = tk.Entry(self.GeneralWin, textvariable = self.eG0variable)
		self.eG0.grid(row=2, column=1, columnspan=2)


		self.eG1 = tk.Button(self.GeneralWin, text = "Browse...", command= lambda: getParams(self.eG1variable))
		self.eG1.grid(row=3, column=1, columnspan=1)
		self.eG1name = tk.Entry(self.GeneralWin, textvariable = self.eG1variable, state = "readonly")
		self.eG1name.grid(row=3,column=2, columnspan=1)


		self.eG2 = tk.Button(self.GeneralWin, text = "Browse...", command= lambda: getCoords(self.eG2variable))
		self.eG2.grid(row=4, column=1, columnspan=1)
		self.eG2name = tk.Entry(self.GeneralWin, textvariable = self.eG2variable, state = "readonly")
		self.eG2name.grid(row=4,column=2, columnspan=1)


		self.eG3 = tk.Entry(self.GeneralWin, textvariable = self.eG3variable)
		self.eG3.grid(row=5, column=1, columnspan=2)


		self.eG4 = tk.Entry(self.GeneralWin, textvariable = self.eG4variable)
		self.eG4.grid(row=6, column=1, columnspan=2)


		self.eG6a = tk.Radiobutton(self.GeneralWin, text="True", variable=self.eG6variable, value=True)
		self.eG6a.grid(row=8, column=1)
		self.eG6b = tk.Radiobutton(self.GeneralWin, text="False", variable=self.eG6variable, value=False)
		self.eG6b.grid(row=8, column=2)


		self.eG7 = tk.Entry(self.GeneralWin,textvariable = self.eG7variable)
		self.eG7.grid(row=9, column=1, columnspan=2)
		
		self.eG8a = tk.Radiobutton(self.GeneralWin, text="True", variable=self.eG8variable, value=True)
		self.eG8a.grid(row=10, column=1)
		self.eG8b = tk.Radiobutton(self.GeneralWin, text="False", variable=self.eG8variable, value=False)
		self.eG8b.grid(row=10, column=2)

		#Make the accept/cancel button, cancel forces a reinstancing of the values, and sets them back to the default
		#update: in the end it's just an add button

		fB1 = tk.Button(self.GeneralWin, text = "Set", command= lambda: add("G0", self.GeneralWin, self.eExpvariable.get()), width=10) 
		fB1.grid(row=12,column=0, columnspan=3)


	def openThermodynWindow(self):
		self.ThermoWin = tk.Toplevel()
		self.ThermoWin.title("")

		ThermoTitle = tk.Label(self.ThermoWin, text = "Thermodynamic Parameters", bd=2, relief=tk.SOLID, pady=5, padx=5)
		ThermoTitle.grid(row=0, column=0, columnspan=3)

		space=tk.Label(self.ThermoWin, text=" ")
		space.grid(row=1,column=0,columnspan=3)

		#Make the labels
		StatusLabelsDictT = {"lT0":"Thermostat:", "lT1":"Initial Temp:", "lT2":"Final Temp:",	"lT3":"Boost Temp:", "lT4":"Forcefield Scale:", "lT6":"Custom FF Scale Factor:", "lT7":"GBSA Scale Factor"} 
		StatusLabelslT = ["lT0", "lT1", "lT2", "lT3", "lT4", "lT6", "lT7"]
	
		ThermoCounter = 2	
		for LabelObj in StatusLabelslT:
			LabelObj = tk.Label(self.GeneralWin, text=StatusLabelsDictT[LabelObj])
			LabelObj.grid(row=ThermoCounter, column=0)
			ThermoCounter = ThermoCounter + 1

		#make the fields
		
		global T0
		if (T0 == 1):
			pass
		else:
			#again, we have to build arrays (see general ~313)
			self.eT0variableArray = np.full(self.eExpNumber.get(),fill_value=default_text ,dtype='object')	
			self.eT1variableArray = np.full(self.eExpNumber.get(),fill_value=-1 ,dtype='object')	
			self.eT2variableArray = np.full(self.eExpNumber.get(),fill_value=-1 ,dtype='object')	
			self.eT3variableArray = np.full(self.eExpNumber.get(),fill_value=-1 ,dtype='object')	
			self.eT4variableArray = np.full(self.eExpNumber.get(),fill_value=default_text ,dtype='object')	
			self.eT6variableArray = np.full(self.eExpNumber.get(),fill_value=-1 ,dtype='object')	
					
			self.eT0variable = tk.StringVar()
			self.eT0variable.set("Andersen")
			self.eT1variable = tk.IntVar()	
			self.eT1variable.set(-1)
			self.eT2variable = tk.IntVar()	
			self.eT2variable.set(-1)
			self.eT3variable = tk.IntVar()	
			self.eT3variable.set(-1)
			self.eT4variable = tk.StringVar()
			self.eT4variable.set("AMBER")
			self.eT6variable = tk.DoubleVar()	
			self.eT6variable.set(-1)
			T0 = 1	#so we don't reinstantiate every field every time
		
		self.eT0a = tk.Radiobutton(self.ThermoWin, text="None", variable=self.eT0variable, value="None")
		self.eT0a.grid(row=2, column=1)
		self.eT0b = tk.Radiobutton(self.ThermoWin, text="Andersen", variable=self.eT0variable, value="Andersen")
		self.eT0b.grid(row=2, column=2)


		self.eT1 = tk.Entry(self.ThermoWin, textvariable = self.eT1variable)
		self.eT1.grid(row=3, column=1, columnspan=2)
		

		self.eT2 = tk.Entry(self.ThermoWin, textvariable = self.eT2variable)
		self.eT2.grid(row=4, column=1, columnspan=2)


		self.eT3 = tk.Entry(self.ThermoWin, textvariable = self.eT3variable)
		self.eT3.grid(row=5, column=1, columnspan=2)
		

		self.eT4a = tk.Radiobutton(self.ThermoWin, text="AMBER", variable=self.eT4variable, value="AMBER", command = lambda: disableEntry(self.eT5))
		self.eT4a.grid(row=6, column=1)
		self.eT4b = tk.Radiobutton(self.ThermoWin, text="Custom", variable=self.eT4variable, value="Custom", command = lambda: enableEntry(self.eT5))
		self.eT4b.grid(row=6, column=2)


		self.eT5 = tk.Entry(self.ThermoWin, textvariable=self.eT4variable, state=tk.DISABLED)
		self.eT5.grid(row=7, column=1, columnspan=2)
		

		self.eT6 = tk.Entry(self.ThermoWin, textvariable = self.eT6variable)
		self.eT6.grid(row=8, column=1, columnspan=2)


		fB1 = tk.Button(self.ThermoWin, text = "Set", command= lambda: add("T0", self.ThermoWin, self.eExpvariable.get()), width=10) 
		fB1.grid(row=9,column=0, columnspan=3)


	def openWorldWindow(self):		
		
		if (G0 == 0):
			messagebox.showwarning("Oops!", "Please select the number of worlds before setting their parameters.")

		elif (any(i <= 0 for i in self.eG3variableArray)):
			messagebox.showwarning("Oops!", "There is at least one experiment with 0 or fewer worlds.")
		
		else:		

			global DidWorlds 
			DidWorlds = 0

			self.WorldWin = tk.Toplevel()
			self.WorldWin.title("")
		
			self.eW0variable = tk.IntVar() #variable that keeps track of selected world

			WorldTitle = tk.Label(self.WorldWin, text = "World Parameters", bd=2, relief=tk.SOLID, pady=5, padx=5)
			WorldTitle.grid(row=0, column=1, columnspan=4)

			space=tk.Label(self.WorldWin, text=" ")
			space.grid(row=1,column=0,columnspan=3)
		
			StatusLabelsDictW = {"lW0":"World", "lW1":"World Type", "lW2":"Run Type",	"lW3":"Flex Type", "lW4":"Flex File", "lW4a":"Rigid Bodies File", "lW5":"Reblock Frequency", "lW5a":"Roots", "lW6":"Sampler", "lW7":"Timestep", "lW8":"MD Steps", "lW8a":"MD Steps StDev", "lW9":"Boost MD Steps", "lW10":"Samples Per Round", "lW11":"Seed", "lW12":"Use Fixman Potential", "lW13":"Use Fixman Torque", "lW14":"Root Mobilities"} 
			StatusLabelslW = ["lW0", "lW1", "lW2", "lW3", "lW4", "lW4a", "lW5", "lW5a", "lW6", "lW7", "lW8", "lW8a", "lW9", "lW10", "lW11", "lW12", "lW13", "lW14"]
	
			WorldCounter = 2	
			for LabelObj in StatusLabelslW:
				LabelObj = tk.Label(self.GeneralWin, text=StatusLabelsDictW[LabelObj])
				LabelObj.grid(row=WorldCounter, column=0)
			
			self.setNumberofWorlds(self.eExpvariable.get()) #first time
			self.eExpvariable.trace("w", lambda name, index, mode: self.changeNumberofWorlds(self.eExpvariable.get())) #I couldn't figure out a way to remake the options menu without closing/opening it

			eW0 = tk.OptionMenu(self.WorldWin, self.eW0variable, *self.worlds)
			eW0.grid(row=2,column=2,columnspan=2)

			
			
			#rest of fields
		
			#making a 2-D array would be complicated, and it would be too confusing for the user to have to juggle 2 variables. As such, we will have to re-open the worlds window for each experiment. 
			#this does not mean we don't keep track of the worlds for each experiment, of course.	
			global W0 
			if (W0 == 1):
				#We need to create arrays that hold all the values. These arrays will be of dimension (np.max(self.eG3variableArray), self.eExpNumber.get()). At first it will be filled with default text
				# Then, as each experiment fills in the correct number of world parameters, the default_text will be replaced by actual values. When it comes to the update screen/ writing the input
				# we will ignore all fields for a certain experiment past the (self.eG3variableArray(ExperimentNumber)) for that particular experiment. It is *very* important we do not reinstate these ones, else the values get lost.
				pass	
						
			else:
		
				#TODO: REDO THIS OH MY GOD 	
				self.eW1variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW2variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW3variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW4variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW4avariableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW5variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW5avariableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW6variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW7variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW8variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW8avariableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW9variableGlobalArray  = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW10variableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW11variableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW12variableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW13variableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')
				self.eW14variableGlobalArray = np.full((np.max(self.eG3variableArray),self.eExpNumber.get()),fill_value=default_text, dtype='object')

				W0 = 1
				self.updateShow = 0



			self.eW1variable = tk.StringVar()
			self.eW1variable.set("IC")
			self.eW2variable = tk.StringVar()
			self.eW2variable.set("Normal")
			self.eW3variable = tk.StringVar()
			self.eW3variable.set("From File")
			self.eW4variable = tk.StringVar()
			self.eW4variable.set(default_text)
			self.eW4avariable = tk.StringVar()
			self.eW4avariable.set(default_text)
			self.eW5variable = tk.IntVar()
			self.eW5variable.set(0)
			self.eW5avariable = tk.IntVar()
			self.eW5avariable.set(0)
			self.eW6variable = tk.StringVar()
			self.eW6variable.set("HMC")
			self.eW7variable = tk.DoubleVar()
			self.eW7variable.set(0.002)
			self.eW8variable = tk.IntVar()
			self.eW8variable.set(0)
			self.eW8avariable = tk.IntVar()
			self.eW8avariable.set(0)
			self.eW9variable = tk.IntVar()
			self.eW9variable.set(0)
			self.eW10variable = tk.IntVar()
			self.eW10variable.set(0)
			self.eW11variable = tk.StringVar()
			self.eW11variable.set(self.eW0variable.get())		#the seed should be different for each world by default
			self.eW12variable = tk.BooleanVar()
			self.eW12variable.set(True)
			self.eW13variable = tk.BooleanVar()
			self.eW13variable.set(True)
			self.eW14variable = tk.StringVar()
			self.eW14variable.set(RootMobilities[2])
			
			W0 = 1 #so we don't reinstantiate all fields every time
	

			self.eW1a = tk.Radiobutton(self.WorldWin, text="IC", variable=self.eW1variable, value="IC")
			self.eW1a.grid(row=3, column=1)
			self.eW1b = tk.Radiobutton(self.WorldWin, text="RB", variable=self.eW1variable, value="RB")
			self.eW1b.grid(row=3, column=2)
			self.eW1c = tk.Radiobutton(self.WorldWin, text="TD", variable=self.eW1variable, value="TD")
			self.eW1c.grid(row=3, column=3)

		
			self.eW2a = tk.Radiobutton(self.WorldWin, text="Normal", variable=self.eW2variable, value="Normal")
			self.eW2a.grid(row=4, column=1)
			self.eW2b = tk.Radiobutton(self.WorldWin, text="Non-Eq", variable=self.eW2variable, value="Non-Eq", command = lambda: NotYet(self.eW2variable, "Normal"))
			self.eW2b.grid(row=4, column=3)
		
	
			self.eW3a = tk.Radiobutton(self.WorldWin, text="From File", variable=self.eW3variable, value="From File", command = lambda: disableEntry(self.eW5))
			self.eW3a.grid(row=5, column=1)
			self.eW3b = tk.Radiobutton(self.WorldWin, text="Auto Reblock", variable=self.eW3variable, value="Auto Reblock", command = lambda: NotYet(self.eW3variable, "From File"))
			self.eW3b.grid(row=5, column=3)


			self.eW4 = tk.Button(self.WorldWin, text = "Browse...", command= lambda: getFlex(self.eW4variable))
			self.eW4.grid(row=6, column=1,)
			self.eW4name = tk.Entry(self.WorldWin, textvariable = self.eW4variable, state = "readonly")
			self.eW4name.grid(row=6,column=2, columnspan=3)

			self.eW4a = tk.Button(self.WorldWin, text = "Browse...", command= lambda: getRB(self.eW4avariable))
			self.eW4a.grid(row=7, column=1,)
			self.eW4aname = tk.Entry(self.WorldWin, textvariable = self.eW4avariable, state = "readonly")
			self.eW4aname.grid(row=7,column=2, columnspan=3)
			
			self.eW5 = tk.Entry(self.WorldWin, textvariable = self.eW5variable, state= tk.DISABLED)
			self.eW5.grid(row=8, column=1, columnspan=3)
				
			self.eW5a = tk.Entry(self.WorldWin, textvariable = self.eW5avariable)
			self.eW5a.grid(row=9, column=1, columnspan=3)

			self.eW6a = tk.Radiobutton(self.WorldWin, text="HMC", variable=self.eW6variable, value="HMC")
			self.eW6a.grid(row=10, column=1)
			self.eW6b = tk.Radiobutton(self.WorldWin, text="MC", variable=self.eW6variable, value="MC")
			self.eW6b.grid(row=10, column=3)

			self.eW7 = tk.Entry(self.WorldWin, textvariable = self.eW7variable)
			self.eW7.grid(row=11, column=1, columnspan=3)

			self.eW8 = tk.Entry(self.WorldWin, textvariable = self.eW8variable)
			self.eW8.grid(row=12, column=1, columnspan=3)
			
			self.eW8a = tk.Entry(self.WorldWin, textvariable = self.eW8avariable)
			self.eW8a.grid(row=13, column=1, columnspan=3)
			
			self.eW9 = tk.Entry(self.WorldWin, textvariable = self.eW9variable)
			self.eW9.grid(row=14, column=1, columnspan=3)

			self.eW10 = tk.Entry(self.WorldWin, textvariable = self.eW10variable)
			self.eW10.grid(row=15, column=1, columnspan=3)

			#Define a Tracer object so the seed changes when the world changes (to streamline the process, if the user doesn't want to manually set a seed. This applies if reproducible is true, otherwise the seed is going to be a random number
	
#			print (self.eG8variableArray[self.eExpvariable.get()])
			
			if (self.eG8variableArray[self.eExpvariable.get()] == "True"):
#				print ("Reproducible is true")
				self.eW0variable.trace("w", lambda name, index, mode: changeValue(self.eW11variable, self.eW0variable.get()))
				self.eW11 = tk.Entry(self.WorldWin, textvariable = self.eW11variable)
			else:
				self.eW11variable.set("Seed will be random")
				self.eW11 = tk.Entry(self.WorldWin, textvariable = self.eW11variable, state= tk.DISABLED)
		
			self.eW11.grid(row=16, column=1, columnspan=3)

			self.eW12a = tk.Radiobutton(self.WorldWin, text="True", variable=self.eW12variable, value=True)
			self.eW12a.grid(row=17, column=1)
			self.eW12b = tk.Radiobutton(self.WorldWin, text="False", variable=self.eW12variable, value=False)
			self.eW12b.grid(row=17, column=3)

			self.eW13a = tk.Radiobutton(self.WorldWin, text="True", variable=self.eW13variable, value=True)
			self.eW13a.grid(row=18, column=1)
			self.eW13b = tk.Radiobutton(self.WorldWin, text="False", variable=self.eW13variable, value=False)
			self.eW13b.grid(row=18, column=3)

		
			self.eG14 = tk.OptionMenu(self.WorldWin, self.eW14variable, *RootMobilities, command = self.eW14variable.set)
			self.eG14.grid(row=19, column=1, columnspan=3)




			fB1 = tk.Button(self.WorldWin, text = "Set", command= lambda: add("W0", self.WorldWin, self.eExpvariable.get()), width=10) 
			fB1.grid(row=21,column=1, columnspan=2)


			self.WorldWin.protocol("WM_DELETE_WINDOW", lambda: WorldParamsClose(self.WorldWin)) #we do this so the Worlds window can refresh when we change the experiment, but so it doesn't reappear after the user closes it



			#We also want a "Worlds Status Window". each column  will represent a world, and each row a parameter. Then an empty row which has labels for the nextset of params, then the next experiment
			#there are 16 parameters/world and an unknown number of worlds per experiment. the rows should be 17 + i*17 for the i'th experiment

			if (self.updateShow == 0): #so this window doesn't keep appearing whenever we change the experiment
				self.WorldParams = tk.Toplevel()
				self.WorldParams.title("Worlds Overview")
				DefaultWorldText = tk.StringVar()
				DefaultWorldText.set("<PARAM NOT SET>")
			

			#Unlike the general parameter status window, here the number of labels or fields isn't known a priori, so it's easier to generate them in a loop

				WorldLabels = ["Experiment", "World Type", "Run Type","Flex Type", "Flex File", "Rigid Bodies File", "Reblock Frequency", "Roots", "Sampler", "Timestep", "MD Steps", "MD Steps StDev", "Boost MD Steps", "Samples Per Round", "Seed", "Use Fixman Potential", "Use Fixman Torque", "Root Mobilities"]

				for i in range(self.eExpNumber.get()):
					offset = i*17
					for j in range(self.eG3variableArray[i]):
							label = tk.Label(self.WorldParams, text = "Experiment " + str(i), bd=1, relief=tk.SOLID)
							label.grid(row=+offset,column=0)
					
						for LabelIx in range(1,len(WorldLabels):
							label = tk.Label(self.WorldParams, text = WorldLabels[LabelIx])
							label.grid(row=LabelIx+offset,column=0)
						
						##Fields (column > 0)
						entry   = tk.Label(self.WorldParams, text = " ")  #empty space that corresponds to the "Experiment i" row
						entry.grid(row=0+offset,column=j+1)
						for z in range (1, len(WorldLabels)):
							entry  = tk.Entry(self.WorldParams, textvariable  = DefaultWorldText, state="readonly", justify=tk.CENTER, fg="red")		#A much easier way to write entry fields for the world parameters
							entry.grid(row=z+offset,column=j+1)
							
						
				self.updateShow = 1
				
				
				self.WorldParams.protocol("WM_DELETE_WINDOW", self.WorldParams.iconify) #we do this so the user doesn't accidentally close this window.


	
			


	def openOutputWindow(self):
		self.OutputWin = tk.Toplevel()
		self.OutputWin.title("")
		
		WorldTitle = tk.Label(self.OutputWin, text = "Output Parameters", bd=2, relief=tk.SOLID, pady=5, padx=5)
		WorldTitle.grid(row=0, column=0, columnspan=3)

		space=tk.Label(self.OutputWin, text=" ")
		space.grid(row=1,column=0,columnspan=3)

		
		lX0 = tk.Label(self.OutputWin, text = "Data Print Frequency:")
		lX0.grid(row=2,column=0)
		lX1 = tk.Label(self.OutputWin, text = "Coordinates Print Frequency:")
		lX1.grid(row=3,column=0)
		lX2 = tk.Label(self.OutputWin, text = "Print Geometrical Features:")
		lX2.grid(row=4,column=0)
		lX3 = tk.Label(self.OutputWin, text = "Which Geometrical Features to print?")
		lX3.grid(row=5,column=0)
		lX4 = tk.Label(self.OutputWin, text = "Use Visual Mode?")
		lX4.grid(row=6,column=0)

		
		global X0
		if (X0 == 1):
			pass
		else:
			self.eX0variableArray = np.full(self.eExpNumber.get(), fill_value=-1,    dtype='object')	
			self.eX1variableArray = np.full(self.eExpNumber.get(), fill_value=-1,    dtype='object')	
			self.eX2variableArray = np.full(self.eExpNumber.get(), fill_value=False, dtype='object')	
			self.eX3variableArray = [] 	#this ones a list cause we don't know a priori how many elements each experiment will have.
			self.eX4variableArray = np.full(self.eExpNumber.get(), fill_value=False, dtype='object')	
		
			self.eX0variable = tk.IntVar()
			self.eX0variable.set(-1)     
			self.eX1variable = tk.IntVar()	
			self.eX1variable.set(-1)
			self.eX2variable = tk.BooleanVar()	
			self.eX2variable.set(False)
			self.eX3variable = tk.StringVar()	
			self.eX3variable.set("") #no defaults
			self.eX4variable = tk.BooleanVar()
			self.eX4variable.set(True)    
	 
			for i in range(self.eExpNumber.get()):
				self.eX3variableArray.append([])
			X0 = 1


		self.eX0 = tk.Entry(self.OutputWin, textvariable = self.eX0variable)
		self.eX0.grid(row=2, column=1, columnspan=3)
		
		self.eX1 = tk.Entry(self.OutputWin, textvariable = self.eX1variable)
		self.eX1.grid(row=3, column=1, columnspan=3)

		self.eX2a = tk.Radiobutton(self.OutputWin, text="True", variable=self.eX2variable, value=True, command = lambda: enableEntry(self.eX3))
		self.eX2a.grid(row=4, column=1)
		self.eX2b = tk.Radiobutton(self.OutputWin, text="False", variable=self.eX2variable, value=False, command = lambda: disableEntry(self.eX3))
		self.eX2b.grid(row=4, column=3)

		self.eX3 = tk.Button(self.OutputWin, text = "Add", command = lambda: add2List(self.eX3variable.get(), self.eX3variableArray[self.eExpvariable.get()]), state=tk.DISABLED)
		self.eX3.grid(row=5, column=3)
		self.eX3name = tk.Entry(self.OutputWin, textvariable = self.eX3variable)
		self.eX3name.grid(row=5,column=1, columnspan=2)
		
		self.eX4a = tk.Radiobutton(self.OutputWin, text="True", variable=self.eX4variable, value=True)
		self.eX4a.grid(row=6, column=1)
		self.eX4b = tk.Radiobutton(self.OutputWin, text="False", variable=self.eX4variable, value=False)
		self.eX4b.grid(row=6, column=3)




		fB1 = tk.Button(self.OutputWin, text = "Set", command = lambda: add("X0", self.OutputWin, self.eExpvariable.get()), width=10) 
		fB1.grid(row=8,column=0, columnspan=3)


	def writeAllToFile(self, outFile, numberOfFiles, experimentNumber):		#Given that now all of the variables are stored in arrays, we will use "self.eExpvariable" to figure out which values to write in what file and the corresponding arrays

		global G0, T0, W0, X0

		if (G0 == 0):											#we check if the user opened all windows before trying to write to file
			messagebox.showwarning("Oops!", "General parameters not set. Please set them before trying to write.")
		elif (T0 == 0):
			messagebox.showwarning("Oops!", "Thermodynamic parameters not set. Please set them before trying to write.")
		elif (W0 == 0):
			messagebox.showwarning("Oops!", "World parameters not set. Please set them before trying to write.")
		elif (X0 == 0):
			messagebox.showwarning("Oops!", "Output parameters not set. Please set them before trying to write.")
		
		else:	
			for j in range(numberOfFiles):
				self.eG0variableArray[experimentNumber] = self.eG0variableArray[experimentNumber].replace(" ", "")
				outF = open(outFile + "." + str(j) + ".RS.inp", "w")
				self.distances = []
				self.angles = []
				self.dihedrals = []
				outF.write("#This script has been automatically generated by RIG\n\n")	
				outF.write("#General Info\n\n")
				outF.write("MOLECULES " + str(self.eG0variableArray[experimentNumber]) + "\n")
				outF.write("PRMTOP " + str(self.eG1variableArray[experimentNumber] + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("INPCRD " + str(self.eG2variableArray[experimentNumber] + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("RBFILE ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW4avariableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("FLEXFILE ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW4variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("OUTPUT_DIR " + str(self.eG0variableArray[experimentNumber]) + ".outDir" + "\n\n")
			
				outF.write("#Simulation Parameters\n")
				outF.write("RUN_TYPE ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW2variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("ROOT_MOBILITY ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW14variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("ROUNDS " + str(self.eG4variableArray[experimentNumber]) + "\n")
				outF.write("ROUNDS_TILL_REBLOCK ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW5variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("WORLDS ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW1variableGlobalArray[i][experimentNumber]) + str(i) + " ")
				outF.write("\n")
				outF.write("ROOTS ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW5avariableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("SAMPLER ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW6variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("TIMESTEPS ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW7variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("MDSTEPS ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW8variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("MDSTEPS_STD ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW8avariableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("BOOST_MDSTEPS ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW9variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("SAMPLES_PER_ROUND ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW10variableGlobalArray[i][experimentNumber]) + " ")
				outF.write("\n")
				outF.write("REPRODUCIBLE " + (str(self.eG8variableArray[experimentNumber]).upper() + " ") * self.eG3variableArray[experimentNumber] + "\n")
				if (self.eG8variableArray[experimentNumber] == "True"):
					outF.write("SEED ")
					for i in range(self.eG3variableArray[experimentNumber]):
						outF.write(str(int(self.eW11variableGlobalArray[i][experimentNumber]) + j) + " ") #we add "j" in order to have different seed values for the same world across multiple repetitions
				else:
					pass
				
	
				outF.write("\n\n")
			
				outF.write("#Thermodynamic Parameters\n")
				outF.write("THERMOSTAT " + str(self.eT0variableArray[experimentNumber] + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("TEMPERATURE_INI " + (str(self.eT1variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("TEMPERATURE_FIN " + (str(self.eT2variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("BOOST_TEMPERATURE " + (str(self.eT3variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("FFSCALE " + (str(self.eT4variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("GBSA " + (str(self.eT6variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n\n")
	
				outF.write("#Generalized Coordinates Parameters\n")
				outF.write("FIXMAN_POTENTIAL ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW12variableGlobalArray[i][experimentNumber]).upper() + " ")
				outF.write("\n")
				outF.write("FIXMAN_TORQUE ")
				for i in range(self.eG3variableArray[experimentNumber]):
					outF.write(str(self.eW13variableGlobalArray[i][experimentNumber]).upper() + " ")
				outF.write("\n\n")
			
				outF.write("#Output Parameters\n")
				outF.write("VISUAL " + (str(self.eX4variableArray[experimentNumber]).upper() + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("PRINT_FREQ " + (str(self.eX0variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("WRITEPDBS " + (str(self.eX1variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("GEOMETRY " + (str(self.eX2variableArray[experimentNumber]).upper() + " ") * self.eG3variableArray[experimentNumber] + "\n")
				if (self.eX2variableArray[experimentNumber] == "True"):
					for i in range (len(self.eX3variableArray[experimentNumber])): 				#We get which geometric feature is what kind (distance, angle, dihedral)
						if (len(self.eX3variableArray[experimentNumber][i].split()) == 2):
							self.distances.append(self.eX3variableArray[experimentNumber][i].split())		
						elif (len(self.eX3variableArray[experimentNumber][i].split()) == 3):
							self.angles.append(self.eX3variableArray[experimentNumber][i].split())		
						elif (len(self.eX3variableArray[experimentNumber][i].split()) == 4):
							self.dihedrals.append(self.eX3variableArray[experimentNumber][i].split())
	
	
#					print (self.distances)
#					print (self.angles)
#					print (self.dihedrals)
					if (len(self.distances) > 0):
						outF.write("DISTANCE ")
#						dist = [str(i) for i in str(self.distances).split()]
						for i in range(len(self.distances)):
							for j in range(2):
								outF.write(self.distances[i][j] + " ")
						outF.write("\n")
					if (len(self.angles) > 0):
						outF.write("ANGLE ")	
#						angl = [str(i) for i in str(self.angles).split()]
						for i in range(len(self.angles)):
							for j in range(3):
								outF.write(self.angles[i][j] + " ")
						outF.write("\n")
					if (len(self.dihedrals) > 0):
						outF.write("DIHEDRAL ")
#						dihe = [str(i) for i in str(self.dihedrals).split()]
						for i in range(len(self.dihedrals)):
							for j in range(4):
								outF.write(self.dihedrals[i][j] + " ")
						outF.write("\n\n")

#					print (dist)
#					print (angl)
#					print (dihe)
			
				outF.write("#Software Parameters\n")
				outF.write("THREADS " + (str(self.eG7variableArray[experimentNumber]) + " ") * self.eG3variableArray[experimentNumber] + "\n")
				outF.write("OPENMM " + (str(self.eG6variableArray[experimentNumber]).upper() + " ") * self.eG3variableArray[experimentNumber] + "\n\n")
		
				outF.write("#End of input file.\n")
			
				outF.close()	
		
		
		
	def updateStatus(self, variable):
		global G0,T0,W0,X0			
		
		#assign general values (ONLY if the G0 value is > 0, so that we don't call on variables that haven't been set yet 
		if (variable == "G0" and G0 == 1):		
			for i in range(self.eExpNumber.get()):
				G0value = tk.StringVar()
				if (self.eG0variableArray[i] != default_text):
					G0value.set(self.eG0variableArray[i])
					G0valueSet =  tk.Entry(self.StatusWindow, textvariable = G0value, state = "readonly", fg="black", justify=tk.CENTER)
					G0valueSet.grid(row=1, column=i+1)
				
				G1value = tk.StringVar()
				if (self.eG1variableArray[i] != default_text):
					G1value.set(self.eG1variableArray[i])
					G1valueSet =  tk.Entry(self.StatusWindow, textvariable = G1value, state = "readonly", fg="black", justify=tk.CENTER)
					G1valueSet.grid(row=2, column=i+1)
				
				G2value = tk.StringVar()
				if (self.eG2variableArray[i] != default_text):
					G2value.set(self.eG2variableArray[i])
					G2valueSet =  tk.Entry(self.StatusWindow, textvariable = G2value, state = "readonly", fg="black", justify=tk.CENTER)
					G2valueSet.grid(row=3, column=i+1)
				
				G3value = tk.IntVar()
				if (self.eG3variableArray[i] != -1):
					G3value.set(self.eG3variableArray[i])
					G3valueSet =  tk.Entry(self.StatusWindow, textvariable = G3value, state = "readonly", fg="black", justify=tk.CENTER)
					G3valueSet.grid(row=4, column=i+1)
				
				G4value = tk.IntVar()
				if (self.eG4variableArray[i] != -1):
					G4value.set(self.eG4variableArray[i])
					G4valueSet =  tk.Entry(self.StatusWindow, textvariable = G4value, state = "readonly", fg="black", justify=tk.CENTER)
					G4valueSet.grid(row=5, column=i+1)
				
				G6value = tk.StringVar()
				if (self.eG6variableArray[i] == "True" or self.eG6variableArray[i] == "False"):
					G6value.set(self.eG6variableArray[i])
					G6valueSet =  tk.Entry(self.StatusWindow, textvariable = G6value, state = "readonly", fg="black", justify=tk.CENTER)
					G6valueSet.grid(row=6, column=i+1)
				
				G7value = tk.IntVar()
				if (self.eG7variableArray[i] != -1):
					G7value.set(self.eG7variableArray[i])
					G7valueSet =  tk.Entry(self.StatusWindow, textvariable = G7value, state = "readonly", fg="black", justify=tk.CENTER)
					G7valueSet.grid(row=7, column=i+1)
				
				G8value = tk.StringVar() 
				if (self.eG8variableArray[i] == "True" or self.eG8variableArray[i] == "False"):
					G8value.set(self.eG8variableArray[i])
					G8valueSet =  tk.Entry(self.StatusWindow, textvariable = G8value, state = "readonly", fg="black", justify=tk.CENTER)
					G8valueSet.grid(row=8, column=i+1)

		#assign Thermo values
		if (variable == "T0" and T0 == 1):
			for i in range(self.eExpNumber.get()):
				T0value = tk.StringVar()
				if (self.eT0variableArray[i] != default_text):
					T0value.set(self.eT0variableArray[i])
					T0valueSet =  tk.Entry(self.StatusWindow, textvariable = T0value, state = "readonly", fg="black", justify=tk.CENTER)
					T0valueSet.grid(row=10, column=i+1)
				T1value = tk.IntVar()
				if (self.eT1variableArray[i] != -1):
					T1value.set(self.eT1variableArray[i])
					T1valueSet =  tk.Entry(self.StatusWindow, textvariable = T1value, state = "readonly", fg="black", justify=tk.CENTER)
					T1valueSet.grid(row=11, column=i+1)
				T2value = tk.IntVar()
				if (self.eT2variableArray[i] != -1):
					T2value.set(self.eT2variableArray[i])
					T2valueSet =  tk.Entry(self.StatusWindow, textvariable = T2value, state = "readonly", fg="black", justify=tk.CENTER)
					T2valueSet.grid(row=12, column=i+1)
				T3value = tk.IntVar()
				if (self.eT3variableArray[i] != -1):
					T3value.set(self.eT3variableArray[i])
					T3valueSet =  tk.Entry(self.StatusWindow, textvariable = T3value, state = "readonly", fg="black", justify=tk.CENTER)
					T3valueSet.grid(row=13, column=i+1)
				T4value = tk.StringVar()
				if (self.eT4variableArray[i] != default_text):
					T4value.set(self.eT4variableArray[i])
					T4valueSet =  tk.Entry(self.StatusWindow, textvariable = T4value, state = "readonly", fg="black", justify=tk.CENTER)
					T4valueSet.grid(row=14, column=i+1)
				T6value = tk.IntVar()
				if (self.eT6variableArray[i] != -1):
					T6value.set(self.eT6variableArray[i])
					T6valueSet =  tk.Entry(self.StatusWindow, textvariable = T6value, state = "readonly", fg="black", justify=tk.CENTER)
					T6valueSet.grid(row=15, column=i+1)

		#assign World values
		if (variable == "W0" and W0 == 1):
			for i in range(self.eExpNumber.get()):
				for j in range(self.eG3variableArray[i]): 
					offset = i*17 #same as when we generate the empty window (line ~830)
	
					W1value = tk.StringVar()	#start from W1 because W0 is used for the menus
					if (self.eW1variableGlobalArray[j][i] != default_text):
						W1value.set(self.eW1variableGlobalArray[j][i])
						W1valueSet = tk.Entry(self.WorldParams, textvariable = W1value, state = "readonly", fg="black", justify=tk.CENTER)
						W1valueSet.grid(row=1+offset,column=j+1)
					W2value = tk.StringVar()
					if (self.eW2variableGlobalArray[j][i] != default_text):
						W2value.set(self.eW2variableGlobalArray[j][i])
						W2valueSet = tk.Entry(self.WorldParams, textvariable = W2value, state = "readonly", fg="black", justify=tk.CENTER)
						W2valueSet.grid(row=2+offset,column=j+1)
					W3value = tk.StringVar()
					if (self.eW3variableGlobalArray[j][i] != default_text):
						W3value.set(self.eW3variableGlobalArray[j][i])
						W3valueSet = tk.Entry(self.WorldParams, textvariable = W3value, state = "readonly", fg="black", justify=tk.CENTER)
						W3valueSet.grid(row=3+offset,column=j+1)
					W4value = tk.StringVar()
					if (self.eW4variableGlobalArray[j][i] != default_text):
						W4value.set(self.eW4variableGlobalArray[j][i])
						W4valueSet = tk.Entry(self.WorldParams, textvariable = W4value, state = "readonly", fg="black", justify=tk.CENTER)
						W4valueSet.grid(row=4+offset,column=j+1)
					W4avalue = tk.StringVar()
					if (self.eW4avariableGlobalArray[j][i] != default_text):
						W4avalue.set(self.eW4avariableGlobalArray[j][i])
						W4avalueSet = tk.Entry(self.WorldParams, textvariable = W4avalue, state = "readonly", fg="black", justify=tk.CENTER)
						W4avalueSet.grid(row=5+offset,column=j+1)
					W5value = tk.StringVar()
					if (self.eW5variableGlobalArray[j][i] != default_text):
						W5value.set(self.eW5variableGlobalArray[j][i])
						W5valueSet = tk.Entry(self.WorldParams, textvariable = W5value, state = "readonly", fg="black", justify=tk.CENTER)
						W5valueSet.grid(row=6+offset,column=j+1)
					W5avalue = tk.StringVar()
					if (self.eW5avariableGlobalArray[j][i] != default_text):
						W5avalue.set(self.eW5avariableGlobalArray[j][i])
						W5avalueSet = tk.Entry(self.WorldParams, textvariable = W5avalue, state = "readonly", fg="black", justify=tk.CENTER)
						W5avalueSet.grid(row=7+offset,column=j+1)
					W6value = tk.StringVar()
					if (self.eW6variableGlobalArray[j][i] != default_text):
						W6value.set(self.eW6variableGlobalArray[j][i])
						W6valueSet = tk.Entry(self.WorldParams, textvariable = W6value, state = "readonly", fg="black", justify=tk.CENTER)
						W6valueSet.grid(row=8+offset,column=j+1)
					W7value = tk.StringVar()
					if (self.eW7variableGlobalArray[j][i] != default_text):
						W7value.set(self.eW7variableGlobalArray[j][i])
						W7valueSet = tk.Entry(self.WorldParams, textvariable = W7value, state = "readonly", fg="black", justify=tk.CENTER)
						W7valueSet.grid(row=9+offset,column=j+1)
					W8value = tk.StringVar()
					if (self.eW8variableGlobalArray[j][i] != default_text):
						W8value.set(self.eW8variableGlobalArray[j][i])
						W8valueSet = tk.Entry(self.WorldParams, textvariable = W8value, state = "readonly", fg="black", justify=tk.CENTER)
						W8valueSet.grid(row=10+offset,column=j+1)
					W8avalue = tk.StringVar()
					if (self.eW8avariableGlobalArray[j][i] != default_text):
						W8avalue.set(self.eW8avariableGlobalArray[j][i])
						W8avalueSet = tk.Entry(self.WorldParams, textvariable = W8avalue, state = "readonly", fg="black", justify=tk.CENTER)
						W8avalueSet.grid(row=11+offset,column=j+1)
					W9value = tk.StringVar()
					if (self.eW9variableGlobalArray[j][i] != default_text):
						W9value.set(self.eW9variableGlobalArray[j][i])
						W9valueSet = tk.Entry(self.WorldParams, textvariable = W9value, state = "readonly", fg="black", justify=tk.CENTER)
						W9valueSet.grid(row=12+offset,column=j+1)
					W10value = tk.StringVar()
					if (self.eW10variableGlobalArray[j][i] != default_text):
						W10value.set(self.eW10variableGlobalArray[j][i])
						W10valueSet = tk.Entry(self.WorldParams, textvariable = W10value, state = "readonly", fg="black", justify=tk.CENTER)
						W10valueSet.grid(row=13+offset,column=j+1)
					W11value = tk.StringVar()
					if (self.eW11variableGlobalArray[j][i] != default_text):
						W11value.set(self.eW11variableGlobalArray[j][i])
						W11valueSet = tk.Entry(self.WorldParams, textvariable = W11value, state = "readonly", fg="black", justify=tk.CENTER)
						W11valueSet.grid(row=14+offset,column=j+1)
					W12value = tk.StringVar()
					if (self.eW12variableGlobalArray[j][i] != default_text):
						W12value.set(self.eW12variableGlobalArray[j][i])
						W12valueSet = tk.Entry(self.WorldParams, textvariable = W12value, state = "readonly", fg="black", justify=tk.CENTER)
						W12valueSet.grid(row=15+offset,column=j+1)
					W13value = tk.StringVar()
					if (self.eW13variableGlobalArray[j][i] != default_text):
						W13value.set(self.eW13variableGlobalArray[j][i])
						W13valueSet = tk.Entry(self.WorldParams, textvariable = W13value, state = "readonly", fg="black", justify=tk.CENTER)
						W13valueSet.grid(row=16+offset,column=j+1)
					W14value = tk.StringVar()
					if (self.eW14variableGlobalArray[j][i] != default_text):
						W14value.set(self.eW14variableGlobalArray[j][i])
						W14valueSet = tk.Entry(self.WorldParams, textvariable = W14value, state = "readonly", fg="black", justify=tk.CENTER)
						W14valueSet.grid(row=17+offset,column=j+1)


		#assign Output values
		if (variable == "X0" and X0 == 1):
			for i in range(self.eExpNumber.get()):
				X0value = tk.IntVar()
				if (self.eX0variableArray[i] != -1):
					X0value.set(self.eX0variableArray[i])
					X0valueSet =  tk.Entry(self.StatusWindow, textvariable = X0value, state = "readonly", fg="black", justify=tk.CENTER)
					X0valueSet.grid(row=17, column=i+1)
				X1value = tk.IntVar()
				if (self.eX1variableArray[i] != -1):
					X1value.set(self.eX1variableArray[i])
					X1valueSet =  tk.Entry(self.StatusWindow, textvariable = X1value, state = "readonly", fg="black", justify=tk.CENTER)
					X1valueSet.grid(row=18, column=i+1)
				X2value = tk.StringVar()
				if (self.eX2variableArray[i] == "True" or self.eX2variableArray[i] == "False"):
					X2value.set(self.eX2variableArray[i])
					X2valueSet =  tk.Entry(self.StatusWindow, textvariable = X2value, state = "readonly", fg="black", justify=tk.CENTER)
					X2valueSet.grid(row=19, column=i+1)
				X3value = tk.StringVar()
				if (self.eX3variableArray[i] != ""):
					X3valueSet =  tk.Listbox(self.StatusWindow)
					X3valueSet.insert(tk.END, *self.eX3variableArray[i])
					X3valueSet.grid(row=20, column=i+1)
				X4value = tk.StringVar()
				if (self.eX4variableArray[i] == "True" or self.eX4variableArray[i] == "False"):
					X4value.set(self.eX4variableArray[i])
					X4valueSet =  tk.Entry(self.StatusWindow, textvariable = X4value, state = "readonly", fg="black", justify=tk.CENTER)
					X4valueSet.grid(row=21, column=i+1)



				


	def setParams(self, value, variable):
		global G0, T0, W0, X0
			
		if (G0 == 1 and variable == "G0"):
			self.eG0variableArray[value] = self.eG0variable.get()
			self.eG1variableArray[value] = self.eG1variable.get()
			self.eG2variableArray[value] = self.eG2variable.get()
			self.eG3variableArray[value] = self.eG3variable.get()
			self.eG4variableArray[value] = self.eG4variable.get()
			self.eG6variableArray[value] = str(self.eG6variable.get())
			self.eG7variableArray[value] = self.eG7variable.get()
			self.eG8variableArray[value] = str(self.eG8variable.get())
		
		if (T0 == 1 and variable == "T0"):
			self.eT0variableArray[value] = str(self.eT0variable.get())
			self.eT1variableArray[value] = self.eT1variable.get()
			self.eT2variableArray[value] = self.eT2variable.get()
			self.eT3variableArray[value] = self.eT3variable.get()
			self.eT4variableArray[value] = str(self.eT4variable.get())
			self.eT6variableArray[value] = self.eT6variable.get()
			
		if (X0 == 1 and variable == "X0"):
			self.eX0variableArray[value] = (self.eX0variable.get())
			self.eX1variableArray[value] = self.eX1variable.get()
			self.eX2variableArray[value] = str(self.eX2variable.get())
			self.eX4variableArray[value] = str(self.eX4variable.get())

				
	def setWorldParams(self, ExpNumber, WorldNumber):
		self.eW1variableGlobalArray[WorldNumber][ExpNumber]  = self.eW1variable.get()
		self.eW2variableGlobalArray[WorldNumber][ExpNumber]  = self.eW2variable.get()
		self.eW3variableGlobalArray[WorldNumber][ExpNumber]  = self.eW3variable.get()
		self.eW4variableGlobalArray[WorldNumber][ExpNumber]  = self.eW4variable.get()
		self.eW4avariableGlobalArray[WorldNumber][ExpNumber] = self.eW4avariable.get()
		self.eW5variableGlobalArray[WorldNumber][ExpNumber]  = self.eW5variable.get()
		self.eW5avariableGlobalArray[WorldNumber][ExpNumber] = self.eW5avariable.get()
		self.eW6variableGlobalArray[WorldNumber][ExpNumber]  = self.eW6variable.get()
		self.eW7variableGlobalArray[WorldNumber][ExpNumber]  = self.eW7variable.get()
		self.eW8variableGlobalArray[WorldNumber][ExpNumber]  = self.eW8variable.get()
		self.eW8avariableGlobalArray[WorldNumber][ExpNumber] = self.eW8avariable.get()
		self.eW9variableGlobalArray[WorldNumber][ExpNumber]  = self.eW9variable.get()
		self.eW10variableGlobalArray[WorldNumber][ExpNumber] = self.eW10variable.get()
		self.eW11variableGlobalArray[WorldNumber][ExpNumber] = self.eW11variable.get()
		self.eW12variableGlobalArray[WorldNumber][ExpNumber] = str(self.eW12variable.get())
		self.eW13variableGlobalArray[WorldNumber][ExpNumber] = str(self.eW13variable.get())
		self.eW14variableGlobalArray[WorldNumber][ExpNumber] = str(self.eW14variable.get())



	def setNumberofWorlds(self, value):
		self.worlds = []
		for i in range(self.eG3variableArray[value]):
			self.worlds.append(i)
		


	def changeNumberofWorlds(self, value):
		global DidWorlds
		print (DidWorlds)
		if (DidWorlds == 0):

			self.worlds = []
			for i in range(self.eG3variableArray[value]):
				self.worlds.append(i)
		
			self.eW0variable.set(self.worlds[0])
			self.WorldWin.destroy()
			self.openWorldWindow()
		


TheGUI = MainWindow(tk.Tk())
tk.mainloop()
