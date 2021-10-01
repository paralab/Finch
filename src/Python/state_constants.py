from enum import Enum, unique

@unique
class Solvers(Enum):
	CG = 1
	DG = 2
	HDG = 3
    
@unique
class TimeSteppers(Enum):
    EULER_EXPLICIT = 1
    EULER_IMPLICIT = 2
    RK4 = 3

@unique
class OutputFormats(Enum):
    VTK = 1
    RAW = 2
    CUSTOM = 3

@unique
class FunctionSpaces(Enum):
	LEGENDRE = 1

@unique
class NonlinearSolvers(Enum):
	NEWTON_RAPHSON = 1

class GenerationType(Enum):
    CPP = 1
    MATLAB = 2