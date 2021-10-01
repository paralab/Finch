from enum import Enum, unique
from state_constants import *
from domain import *
class State:
	def __init__(self, configDict):
		self.configDict = {}
		for param in RequiredParameters:
			if param in configDict:
				self.setParameter(param, configDict[param])
				del configDict[param]
			else:
				self.setParameter(param)

		#merge in options with string keys
		self.configDict = {**self.configDict, **configDict}

	def setParameter(self, id, value=None):
		key = id
		if isinstance(id, RequiredParameters):
			key = id.name
			if value is None:
				#use default
				value = id.value

		self.configDict[key] = value

	def getParameter(self, id):
		key = id
		if isinstance(id, RequiredParameters):
			key = id.name

		return self.configDict[key]

	def generateCode(self, type=GenerationType.CPP):
		print("Generating Code")

#enum value is the default if unset by user.
#can be simple value, enum, or full class.
class RequiredParameters(Enum):
	Dimensions = 2
	Solver = Solvers.CG
	Domain = Domain([(0,1),(0,1)], {}, Geometry.Square, Decomposition.Unstructured)
	TimeStepper = TimeSteppers.RK4
	OutputFormat = OutputFormats.VTK
	FunctionSpace = FunctionSpaces.LEGENDRE
	NonlinearSolver = NonlinearSolvers.NEWTON_RAPHSON
	#... and many more ...  