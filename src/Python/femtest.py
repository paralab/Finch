from state import *
from state import RequiredParameters as RP
boundaries = [
    Boundary("A", [(0,0), (1,0)], BoundaryConditions.Dirichlet, 0),
    Boundary("B", [(1,0), (1,1)], BoundaryConditions.Dirichlet, 0),
    Boundary("C", [(1,1), (0,1)], BoundaryConditions.Dirichlet, 0),
    Boundary("D", [(0,1), (0,0)], BoundaryConditions.Dirichlet, 0)
]
domain = Domain([(0,1),(0,1)], boundaries, Geometry.Square, Decomposition.Unstructured)

#basic usage. Initialize with a dict of options. Everything left out will be initialized to default value
femstate = State({
    RP.Dimensions : 2,
    RP.Domain: domain,
    #can also just use a magic string instead of enums
    'Solver': Solvers.CG,
    #RequiredParameters.SolutionMethod: SolutionMethods.CG,
    RP.TimeStepper: TimeSteppers.RK4,
    RP.OutputFormat: OutputFormats.VTK

    ############ TODO ##########
    # Some sort of specification of the problem using sympy (weak form, variables, initial conditions, etc)

})

femstate.generateCode(GenerationType.CPP)