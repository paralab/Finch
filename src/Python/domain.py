from enum import Enum, unique

@unique
class Geometry(Enum):
	Irregular = 1
	Square = 2
	Triangle = 3

@unique
class BoundaryConditions(Enum):
	Dirichlet = 1
	Neumann = 2

@unique
class Decomposition(Enum):
	Octree = 1
	Unstructured = 2

class Domain:
	#configure boundaries, boundary conditions, geometry, mesh
	def __init__(self, extent, boundaries, geometry=Geometry.Irregular, decomposition=Decomposition.Unstructured):
		self.extent = extent
		self.boundaries = boundaries
		self.geometry = geometry
		self.decomposition = decomposition
		#...and so on...

class Boundary:
	def __init__(self, id, extent, conditionType=BoundaryConditions.Dirichlet, condition = 0):
		self.id = id
		self.extent = extent
		self.conditionType = conditionType
		self.condition = condition