# Is this a reasonable way to define states?
# These are just some examples, clearly many more would be needed

from enum import Enum

class Finch_constants(Enum):
    # Domain geometry
    SQUARE = 1
    IRREGULAR = 2

    # Domain decomposition
    TREE = 11
    UNSTRUCTURED = 12

    # Solver type
    CG = 21
    DG = 22
    HDG = 23

    NODAL = 31
    MODAL = 32

    # Function space
    LEGENDRE = 41

    # Element node positions
    UNIFORM = 51
    GAUSS = 52
    LOBATTO = 53

    # Nonlinear solver methods
    NONLINEAR_NEWTON = 62
    NONLINEAR_SOMETHING_ELSE = 63

    # Time steppers
    EULER_EXPLICIT = 71
    EULER_IMPLICIT = 72
    RK4 = 73
    LSRK4 = 74
    ABM4 = 75

    # Linear system solvers/structures
    OURS = 81
    PETSC = 82

    # Output format
    VTK = 91
    RAW_OUTPUT = 92
    CUSTOM_OUTPUT = 93

