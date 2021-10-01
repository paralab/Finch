###########################################
# A structure for storing finch comfiguration
###########################################

from Finch_constants import Finch_constants as FS

class Finch_config:
    # How should we define the states?
    # As integers assigned to constant names in an enum?

    # Initializes a default config.
    def __init__(self):
        # Domain 
        self.dimension = 1
        self.geometry = FS.SQUARE
        self.mesh_type = FS.TREE

        # FEM details
        self.solver_type = FS.DG
        self.basis_type = FS.NODAL
        self.trial_function = FS.LEGENDRE
        self.test_function = FS.LEGENDRE
        self.basis_nodes = FS.LOBATTO
        self.p_adaptive = True
        self.basis_order_min = 4
        self.basis_order_max = 8

        # Variables
        # Maybe info about variables and boundaries should be
        # in problem specification rather than here.
        
        # Other solver details
        # These are some basic examples, there are more needed.
        self.linear = True
        self.t_adaptive = True
        self.stepper = FS.RK4
        self.linalg_matrixfree = False
        self.linalg_backend = FS.PETSC
        self.output_format = FS.VTK


