#############################################################################################################
# Python minimal implementation of the function solveCobraLP from Matlab COBRA Toolbox
# Only implements Gurobi case
# Daniel Ribeiro 2021
#############################################################################################################
import numpy as np
from scipy import sparse
from libsmop import *
from solveCobraLP_helper import *

# COBRA - cobrapy
import cobra
# GUROBI
import gurobipy as gp
from gurobipy import GRB

import time


# define GLOBAL variables
# Pytohn str
global CBT_LP_SOLVER
global CBT_MILP_SOLVER
global CBT_QP_SOLVER
global CBT_MIQP_SOLVER
global CBT_NLP_SOLVER
global GUROBI_PATH

# Python dict of (parameter: values)
global CBT_LP_PARAMS
global CBT_MILP_PARAMS
global CBT_QP_PARAMS
global CBT_MIQP_PARAMS
global CBT_NLP_PARAMS

#DEBUG
#comment 'model' only for debug
@function
def solveCobraLP(LPproblem, model, *args, **kwargs):
    # Daniel Ribeiro 2021: Python minimal implementation of the function solveCobraLP from Matlab COBRA Toolbox
    # Only implements Gurobi case
    #
    # Solves constraint-based LP problems
    #
    # USAGE:
    #
    #    solveCobraLP(LPproblem, varargin)
    #
    # Daniel 2021: Use libsmop struct class
    # INPUT:
    #    LPproblem:     Structure containing the following fields describing the LP problem to be solved
    #
    #                     * .A - LHS matrix
    #                     * .b - RHS vector
    #                     * .c - Objective coeff vector
    #                     * .lb - Lower bound vector
    #                     * .ub - Upper bound vector
    #                     * .osense - Objective sense (-1 means maximise (default), 1 means minimise)
    #                     * .csense - Constraint senses, a string containting the constraint sense for
    #                       each row in A ('E', equality, 'G' greater than, 'L' less than).

    # Daniel 2021: Partially implemented
    # OPTIONAL INPUTS:
    #    varargin:      Additional parameters either as parameter struct, or as
    #                   parameter/value pairs. A combination is possible, if
    #                   the parameter struct is either at the beginning or the
    #                   end of the optional input.
    #                   All fields of the struct which are not COBRA parameters
    #                   (see `getCobraSolverParamsOptionsForType`) for this
    #                   problem type will be passed on to the solver in a
    #                   solver specific manner. Some optional parameters which
    #                   can be passed to the function as parameter value pairs,
    #                   or as part of the options struct are listed below:
    #
    #    printLevel:    Printing level
    #
    #                     * 0 - Silent (Default)
    #                     * 1 - Warnings and Errors
    #                     * 2 - Summary information
    #                     * 3 - More detailed information
    #                     * > 10 - Pause statements, and maximal printing (debug mode)
    #
    #    saveInput:     Saves LPproblem to filename specified in field.
    #                   i.e. parameters.saveInput = 'LPproblem.mat';
    #
    #    minNorm:       {(0), scalar , `n x 1` vector}, where `[m, n] = size(S)`;
    #                   If not zero then, minimise the Euclidean length
    #                   of the solution to the LP problem. minNorm ~1e-6 should be
    #                   high enough for regularisation yet maintain the same value for
    #                   the linear part of the objective. However, this should be
    #                   checked on a case by case basis, by optimization with and
    #                   without regularisation.
    #
    #    primalOnly:    {(0), 1}; 1 = only return the primal vector (lindo solvers)
    #
    #    solverParams:  solver-specific parameter structure. Formats supported
    #                   are ILOG cplex and Tomlab parameter syntax. see example
    #                   for details.

    # Daniel 2021: Use libsmop struct class
    # OUTPUT:
    #    solution:      Structure containing the following fields describing a LP solution:
    #                     * .full:         Full LP solution vector
    #                     * .obj:          Objective value
    #                     * .rcost:        Reduced costs, dual solution to :math:`lb <= v <= ub`
    #                     * .dual:         dual solution to `A*v ('E' | 'G' | 'L') b`
    #                     * .solver:       Solver used to solve LP problem
    #                     * .algorithm:    Algorithm used by solver to solve LP problem
    #                     * .stat:         Solver status in standardized form
    #
    #                       * 1 - Optimal solution
    #                       * 2 - Unbounded solution
    #                       * 3 - Partial success (OPTI-csdp) - will not give desired
    #                         result from OptimizeCbModel
    #                       * 0 - Infeasible
    #                       * -1 - No solution reported (timelimit, numerical problem etc)
    #                     * .origStat:     Original status returned by the specific solver
    #                     * .time:         Solve time in seconds
    #                     * .basis:        (optional) LP basis corresponding to solution
    #
    # NOTE:
    #           Optional parameters can also be set through the
    #           solver can be set through `changeCobraSolver('LP', value)`;
    #           `changeCobraSolverParams('LP', 'parameter', value)` function. This
    #           includes the minNorm and the `printLevel` flags.
    #
    # EXAMPLE:
    #
    #    #Optional parameters can be entered in three different ways {A,B,C}
    #
    #    #A) as a generic solver parameter followed by parameter value:
    #    [solution] = solveCobraLP(LPCoupled, 'printLevel', 1);
    #    [solution] = solveCobraLP(LPCoupled, 'printLevel', 1, 'feasTol', 1e-8);
    #
    #    #B) parameters structure with field names specific to a particular solvers
    #    # internal parameter fields
    #    [solution] = solveCobraLP(LPCoupled, parameters);
    #
    #    #C) as parameter followed by parameter value, with a parameter structure
    #    #with field names specific to a particular solvers internal parameter,
    #    #fields as the LAST argument
    #    [solution] = solveCobraLP(LPCoupled, 'printLevel', 1, 'feasTol', 1e-6, parameters);
    #
    # .. Authors:
    #       - Markus Herrgard, 08/29/06
    #       - Ronan Fleming, 11/12/08 'cplex_direct' allows for more refined control
    #       of cplex than tomlab tomrun
    #       - Ronan Fleming, 04/25/09 Option to minimise the Euclidean Norm of internal
    #       fluxes using either 'cplex_direct' solver or 'pdco'
    #       - Jan Schellenberger, 09/28/09 Changed header to be much simpler.  All parameters
    #       now accessed through changeCobraSolverParams(LP, parameter,value)
    #       - Richard Que, 11/30/09 Changed handling of optional parameters to use
    #       getCobraSolverParams().
    #       - Ronan Fleming, 12/07/09 Commenting of input/output
    #       - Ronan Fleming, 21/01/10 Not having second input, means use the parameters as specified in the
    #       global paramerer variable, rather than 'default' parameters
    #       - Steinn Gudmundsson, 03/03/10 Added support for the Gurobi solver
    #       - Ronan Fleming, 01/24/01 Now accepts an optional parameter structure with nonstandard
    #       solver specific parameter options
    #       - Tim Harrington, 05/18/12 Added support for the Gurobi 5.0 solver
    #       - Ronan Fleming, 07/04/13 Reinstalled support for optional parameter structure

    print("Solving LP problem...")
    varargin = solveCobraLP.varargin
    nargin = solveCobraLP.nargin

    # comment only if debugging
    solver = cobra.util.solver.interface_to_str(model.solver.interface)

    # Check that solver is gurobi
    if solver != 'gurobi':
        print("ERORR: Gurobi is not the solver")
        print("\tsolveCobraLPpy only implemented for Gurobi solver")
    else:
        #Initialize global variables for SOLVER
        global CBT_LP_SOLVER
        global CBT_MILP_SOLVER
        global CBT_QP_SOLVER
        global CBT_MIQP_SOLVER
        global CBT_NLP_SOLVER

        CBT_LP_SOLVER   = solver
        CBT_MILP_SOLVER = solver
        CBT_QP_SOLVER   = solver
        CBT_MIQP_SOLVER = solver
        CBT_NLP_SOLVER  = solver

    #Initialize global variables for ProbemType
    global CBT_LP_PARAMS
    global CBT_MILP_PARAMS
    global CBT_QP_PARAMS
    global CBT_MIQP_PARAMS
    global CBT_NLP_PARAMS

    CBT_LP_PARAMS   = getCobraSolverParamsOptionsForType('LP')
    CBT_MILP_PARAMS = getCobraSolverParamsOptionsForType('MILP')
    CBT_QP_PARAMS   = getCobraSolverParamsOptionsForType('QP')
    CBT_MIQP_PARAMS = getCobraSolverParamsOptionsForType('MIQP')
    CBT_NLP_PARAMS  = getCobraSolverParamsOptionsForType('NLP')

    CBT_LP_PARAMS   = getCobraSolverParams('LP', CBT_LP_PARAMS, 'default')
    CBT_MILP_PARAMS = getCobraSolverParams('MILP', CBT_MILP_PARAMS, 'default')
    CBT_QP_PARAMS   = getCobraSolverParams('QP', CBT_QP_PARAMS, 'default')
    CBT_MIQP_PARAMS = getCobraSolverParams('MIQP', CBT_MIQP_PARAMS, 'default')
    CBT_NLP_PARAMS  = getCobraSolverParams('NLP', CBT_NLP_PARAMS, 'default')

    problemTypeParams, solverParams = parseSolverParameters('LP', varargin[1:][:])
    #Solver was set above
    #FIX
    # set the solver
    #solver = problemTypeParams.solver

    # check solver compatibility with minNorm option
    if max(problemTypeParams.minNorm) != None and not any(strcmp(solver, np.asarray(['cplex_direct', 'cplex']))):
        error('minNorm only works for LP solver ''cplex_direct'' from this interface, use optimizeCbModel for other solvers.')

    # save Input if selected
    try:
        fileName = getattr(problemTypeParams, 'saveInput')
    except AttributeError:
        fileName = None
    if fileName != None:    
        if regexp(fileName, '.mat', 'match') == None:
            fileName = ''.join([fileName, '.mat'])
        disp(f'Saving LPproblem in {fileName}')
        with open(fileName, 'w') as f:
            f.write(str(LPproblem))

    # support for lifting of ill-scaled models
    if problemTypeParams.lifting == 1:
        print("Support for lifting of ill-scaled models is not implemented")
        raise NotImplementedError

    # assume constraint matrix is S if no A provided.
    if (not hasattr(LPproblem, 'A')) and hasattr(LPproblem, 'S'):
        print("Assume constraint matrix is S if no A provided")
        setattr(LPproblem, 'A', LPproblem.S)

    # assume constraint A*v = b if csense not provided
    if not hasattr(LPproblem, 'csense'):
        print("Assume constraint A*v = b if csense not provided")
        # if csense is not declared in the model, assume that all
        # constraints are equalities.
        setattr(LPproblem, 'csense', np.full((LPproblem.A.shape[0], 1), 'E'))

    # assume constraint S*v = 0 if b not provided
    if not hasattr(LPproblem, 'b'):
        print("Assume constraint S*v = 0 if b not provided")
        setattr(LPproblem, 'b', np.zeros((LPproblem.A.shape[0], 1)))

    # assume max c'v s.t. S v = b if osense not provided
    if not hasattr(LPproblem, 'osense'):
        print("Assume max c'v s.t. S v = b if osense not provided")
        setattr(LPproblem, 'osense', -1)

    if not hasattr(LPproblem, 'modelID'):
        print("Set modelID")
        setattr(LPproblem, 'modelID', 'aModelID')

    # extract the problem from the structure
    A, b, c, lb, ub, csense, osense, modelID = deal((sparse.csr_matrix(LPproblem.A),
                                                    LPproblem.b,
                                                    LPproblem.c,
                                                    LPproblem.lb,
                                                    LPproblem.ub,
                                                    LPproblem.csense,
                                                    LPproblem.osense,
                                                    LPproblem.modelID))

    if hasattr(LPproblem, 'basis'):
        if not isempty(LPproblem.basis):
            basis = getattr(LPproblem, 'basis')
    else:
        basis = struct('vbasis', None, 'cbasis', None)

    # defaults in case the solver does not return anything
    f = None
    x = None
    y = None
    w = None
    stat = 0
    origStat = None
    origStatText = None
    algorithm = 'default'
    t_start = time.perf_counter()

    ######
    # switch solver
    #   case 'gurobi' in original MATLAB file
    ######
    # Free academic licenses for the Gurobi solver can be obtained from
    # http://www.gurobi.com/html/academic.html
    # resultgurobi = struct('x',[],'objval',[],'pi',[]);

    #  The params struct contains Gurobi parameters. A full list may be
    #  found on the Parameter page of the reference manual:
    #     http://www.gurobi.com/documentation/5.5/reference-manual/node798#sec:Parameters
    #  For example:
    #   params.outputflag = 0;          # Silence gurobi
    #   params.resultfile = 'test.mps'; # Write out problem to MPS file

    # params.method gives the algorithm used to solve continuous models
    # -1=automatic,
    #  0=primal simplex,
    #  1=dual simplex,
    #  2=barrier,
    #  3=concurrent,
    #  4=deterministic concurrent
    # i.e. params.method     = 1;          # use dual simplex method
    if solver == 'gurobi':

        param = solverParams
        if not hasattr(param, 'OutputFlag'):
            if problemTypeParams.printLevel == 0:
                setattr(param, 'OutputFlag', 0)
                setattr(param, 'DisplayInterval', 1)
            elif problemTypeParams.printLevel == 1:
                setattr(param, 'OutputFlag', 0)
                setattr(param, 'DisplayInterval', 1)
            else:
                # silent
                setattr(param, 'OutputFlag', 0)
                setattr(param, 'DisplayInterval', 1)

        if hasattr(param, 'FeasibilityTol'):
            # update tolerance according to actual setting
            setattr(problemTypeParams, 'feasTol', param.FeasibilityTol)
        else:
            setattr(param, 'FeasibilityTol', problemTypeParams.feasTol)

        if hasattr(param, 'OptimalityTol'):
            # update tolerance according to actual setting
            setattr(problemTypeParams, 'optTol', param.OptimalityTol)
        else:
            setattr(param, 'OptimalityTol', problemTypeParams.optTol)

        gurobiLP = {}
        gurobiLP['sense'] = np.full(b.shape, '=').reshape(-1)
        gurobiLP['sense'][csense == 'L'] = '<'
        gurobiLP['sense'][csense == 'G'] = '>'

        # modelsense (optional)
        # The optimization sense. Allowed values are 'min' (minimize) or 'max' (maximize). When absent, the default optimization sense is minimization.
        # gurobipy requires int, not str
        if osense == -1 or osense == 'max':
            gurobiLP['modelsense'] = -1
        else:
            gurobiLP['modelsense'] = 1

        gurobiLP['A'  ] = A   # matix
        gurobiLP['rhs'] = b   # vector
        gurobiLP['lb' ] = lb  # vector
        gurobiLP['ub' ] = ub  # vector
        #gurobi wants a dense double vector as an objective
        gurobiLP['obj'] = c # full

        # basis reuse - Ronan
        if basis.cbasis != None and basis.vbasis !=None:
            gurobiLP['cbasis'] = basis.cbasis.toarray() if sparse.issparse(basis.cbasis) else basis.cbasis
            gurobiLP['vbasis'] = basis.vbasis.toarray() if sparse.issparse(basis.vbasis) else basis.vbasis

        # set the solver specific parameters
        param = updateStructData(param, solverParams)

        # >>>Daniel Ribeiro 2021: Add extra code to reflect gurobipy interface>>>
        gurobi_model = gp.Model("metabolism")
        # Add Vars
        x = gurobi_model.addMVar(shape=gurobiLP['obj'].shape[0],
                                 lb=gurobiLP['lb'],
                                 ub=gurobiLP['ub'],
                                 obj=gurobiLP['obj'])
        gurobi_model.update()
        # Add Contraints
        gurobi_model.addMConstr(A=gurobiLP['A'],
                                x=x,
                                sense=gurobiLP['sense'],
                                b=gurobiLP['rhs'])
        gurobi_model.update()
        # Optimization sense
        gurobi_model.setAttr('modelsense', gurobiLP['modelsense'])
        # Set parameter
        for p in param.__dict__.keys():
            gurobi_model.setParam(p, getattr(param, p))
        gurobi_model.update()
        gurobi_model.optimize()  
        # <<<Daniel Ribeiro 2021: back at Cobra Toolbox<<<

        # see the solvers original status -Ronan
        origStat = gurobi_model.getAttr('status')
        if gurobi_model.getAttr('status') == GRB.OPTIMAL:
            stat = 1; # optimal solution found
            if stat == 1 and isempty(gurobi_model.x):
                    error('solveCobraLP: gurobi reporting OPTIMAL but no solution')
            x,f,y,w,s = deal((gurobi_model.x,
                              gurobi_model.objval,
                              osense * gurobi_model.pi,
                              osense * gurobi_model.rc,
                              gurobi_model.slack))

            if getattr(problemTypeParams, 'printLevel') > 2:
                res1 = A * x + s - b
                print(f"res1 shape: {res1.shape}")
                disp(np.linalg.norm(res1, np.inf))
                res2 = osense * c - A.T * y - w
                disp(np.linalg.norm(res2, np.inf))
                disp("Check osense * c - A.T * lam - w = 0 (stationarity):")
                res22 = gurobiLP['obj'] - gurobiLP['A'].T * gurobi_model.pi - gurobi_model.rc
                disp(res22)
                if not all(res22 < 1e-8):
                    time.sleep(0.1)

                time.sleep(0.1)

                # save the basis
            setattr(basis, 'vbasis', gurobi_model.vbasis)
            setattr(basis, 'cbasis', gurobi_model.cbasis)
        elif gurobi_model.getAttr('status') == GRB.INFEASIBLE:
            stat = 0 # infeasible
        elif gurobi_model.getAttr('status') == GRB.UNBOUNDED:
            stat = 2 # unbounded
        elif gurobi_model.getAttr('status') == GRB.INF_OR_UNBD:
            # we simply remove the objective and solve again.
            # if the status becomes 'OPTIMAL', it is unbounded, otherwise it is infeasible.
            gurobiLP['obj'][:] = 0
            x.setAttr('obj', gurobiLP['obj'])
            gurobi_model.optimize()
            if gurobi_model.getAttr('status') == GRB.OPTIMAL:
                stat = 2
            else:
                stat = 0 # Gurobi reports infeasible *or* unbounded
            end
        else:
            stat = -1 # Solution not optimal or solver problem

        if hasattr(param, 'Method'):
            # -1=automatic,
            #  0=primal simplex,
            #  1=dual simplex,
            #  2=barrier,
            #  3=concurrent,
            #  4=deterministic concurrent
            # i.e. params.method     = 1;          # use dual simplex method
            method = getattr(param, 'Method')
            if method == -1:
                algorithm = 'automatic'
            elif method == 1:
                algorithm = 'primal simplex'
            elif method == 2:
                algorithm = 'dual simplex'
            elif method == 3:
                algorithm = 'barrier'
            elif method == 4:
                algorithm = 'concurrent'
            else:
                algorithm = 'deterministic concurrent'

    else:
        if isempty(solver):
            error("There is no solver for LP problems available")
        else:
            error(''.join(['Unknown solver: ', solver]))

    if stat == -1:
        # this is slow, so only check it if there is a problem
        if any(any(not np.isfinite(A))):
            error("Cannot perform LP on a stoichiometric matrix with NaN of Inf coefficents.")

    if stat == 1 and not strcmp(solver, 'mps'):
        # TODO: pull out slack variable from every solver interface (see list of solvers below)
        if 's' not in globals() or 's' not in locals():
            # slack variables required for optimality condition check, if they are
            # not already provided
            s = b - A * x
            # optimality condition check should still check for satisfaction of the
            # optimality conditions
            s[csense == 'E'] = 0
        else:
            # optimality condition check should still check for satisfaction of the
            # optimality conditions
            s[csense == 'E'] = 0
    else:
        s = None

    if not strcmp(solver, 'cplex_direct') and not strcmp(solver, 'mps'):
    # assign solution
        t = time.perf_counter() - t_start
        solution = struct()
        if 'basis' not in globals() or 'basis' not in locals():
            basis = struct('vbasis', None, 'cbasis', None)
            setattr(solution, 'full', None)
            setattr(solution, 'obj', None)
            setattr(solution, 'rcost', None)
            setattr(solution, 'dual', None)
            setattr(solution, 'slack', None)
            setattr(solution, 'solver', None)
            setattr(solution, 'algorithm', None)
            setattr(solution, 'stat', None)
            setattr(solution, 'origStat', None)
            setattr(solution, 'origStatText', None)
            setattr(solution, 'time', None)
            setattr(solution, 'basis', struct('vbasis', None, 'cbasis', None))
            solution.full, solution.obj, solution.rcost, solution.dual, solution.slack,\
            solution.solver, solution.algorithm, solution.stat, solution.origStat,\
            solution.origStatText, solution.time, solution.basis = deal((x, f, w, y, s, solver, algorithm, stat, origStat, origStatText, t, basis))
    elif strcmp(solver,'mps'):
        solution = struct()

    # check the optimality conditions for various solvers
    if not any(strcmp(solver, 'mps')):
        if solution.stat == 1:
            if not isempty(solution.slack) and not isempty(solution.full):
                # determine the residual 1
                res1 = A * solution.full + solution.slack - b
                res1[logical_not(np.isfinite(res1))] = 0
                tmp1 = np.linalg.norm(res1, inf)
                # evaluate the optimality condition 1
                if tmp1 > (problemTypeParams.feasTol * 1e2):
                    disp(solution.origStat)
                    error(f"[{solver}] Primal optimality condition in solveCobraLP not satisfied, residual = {tmp1}, while feasTol = {problemTypeParams.feasTol}")
                else:
                    if problemTypeParams.printLevel > 0:
                        fprintf(f"\n > [{solver}] Primal optimality condition in solveCobraLP satisfied.")
            if not isempty(solution.rcost) and not isempty(solution.dual) and not any([strcmp(solver, 'glpk'), strcmp(solver, 'matlab')]):
                # determine the residual 2
                res2 = osense * c  - A.T * solution.dual - solution.rcost
                tmp2 = np.linalg.norm(res2, inf)    #TODO matlab linprog still does not pass Testing testDifferentLPSolvers using matlab
                # evaluate the optimality condition 2
                if tmp2 > (problemTypeParams.optTol * 1e2):
                    disp(solution.origStat)
                    if not (length(A) == 1 and strcmp(solver,'pdco')): #todo, why does pdco choke on small A?
                       error(f"[{solver}] Dual optimality condition in solveCobraLP not satisfied, residual = {tmp2}, while optTol = {problemTypeParams.optTol}")
                else:
                    if problemTypeParams.printLevel > 0:
                        fprintf(f"\n > [{solver}] Dual optimality condition in solveCobraLP satisfied.\n")

    return solution


############
# Helper functions
############
# Daniel Ribeiro 2021:
# Bellow are the needed helper functions to implement solveCobraLP in a minimal fashion

@function
def parseSolverParameters(problemType, *args, **kwargs):
    # Daniel Ribeiro 2021: Python implementation of the function parseSolverParameters from Matlab COBRA Toolbox
    # Parse the solver parameters for a specified problem
    #
    # USAGE:
    #    [cobraParams, solverParams] = parseSolverParameters(problemType,varargin)
    #
    # INPUT:
    #    problemType:       The type of the problem to get parameters for
    #                       ('LP','MILP','QP','MIQP','NLP')
    #
    # OPTIONAL INPUTS:
    #    varargin:          Additional parameters either as parameter struct, or as
    #                       parameter/value pairs. A combination is possible, if
    #                       the parameter struct is either at the beginning or the
    #                       end of the optional input.
    #                       All fields of the struct which are not COBRA parameters
    #                       (see `getCobraSolverParamslverParamsOptionsForType`) for this
    #                       problem type will be passed on to the solver in a
    #                       solver specific manner.
    #
    # OUTPUTS:
    #    problemTypeParams: The COBRA Toolbox specific parameters for this
    #                       problem type given the provided parameters
    #
    #    solverParams:      Additional parameters provided which are not part
    #                       of the COBRA parameters and are assumed to be part
    #                       of direct solver input structs.
    varargin = parseSolverParameters.varargin
    nargin = parseSolverParameters.nargin
    varargin = varargin[1:][:]

    cobraSolverParameters = getCobraSolverParamsOptionsForType(problemType); # build the default Parameter Structure

    # set the solver Type
    exec(''.join(['global CBT_', problemType, '_SOLVER']))
    exec(''.join(['global defaultSolver; defaultSolver = CBT_', problemType, '_SOLVER']))

    # get the default variables for the correct solver.
    solverVars = getCobraSolverParams(problemType, cobraSolverParameters, struct('solver', defaultSolver))
    defaultParams = {}
    for key in cobraSolverParameters.keys():
        try:
            defaultParams[key] = solverVars[key]
        except KeyError:
            print(f"Key error (non-existing) when trying to set defaultParams[{key}] = solverVars[{key}]")

    # parse the supplied parameters
    if numel(varargin) > 0:
        # we should have a struct at the end (1 object). #Daniel 2021: The struct comming form easyLP
        if mod(numel(varargin),2) == 1: #>1 args, ( struct)
            optParamStruct = varargin[-1][0]  #get last struct 
            if type(optParamStruct) != struct:
                # but it could also be at the first position, so test that as well.
                optParamStruct = varargin[0]
                varargin = np.delete(varargin, 0, axis=0)
                if type(optParamStruct) != struct:
                    error(''.join(['Invalid Parameters supplied.\n',
                           'Parameters have to be supplied either as parameter/Value pairs, or as struct.\n',
                           'A combination is possible, if the last or first input argument is a struct, and all other arguments are parameter/value pairs']))
            else:
                varargin = np.delete(varargin, -1, axis=1)
        else:
            # no parameter struct. so initialize an empty one.
            optParamStruct = struct()
        # now, loop through all parameter/value pairs.
        for i in range(0, numel(varargin), 2):
            cparam = varargin[i]
            if type(cparam) != str:
                print(f"Type: {type(cparam)}, content: {cparam}")
                error("Parameters have to be supplied as 'parameterName'/Value pairs")
            # the param struct overrides the only use the parameter if it is
            # not a field of the parameter struct.
            if not hasattr(optParamStruct, cparam):
                try:
                    setattr(optParamStruct, cparam, varargin[i+1])
                except:
                    error('All parameters have to be valid matlab field names. %s is not a valid field name' % cparam)
            else:
                print('WARNING: Duplicate parameter %s, both supplied as a field name and a parameter/value pair!' % cparam)
    else:
        # no optional parameters.
        optParamStruct = struct()

    # set up the cobra parameters
    problemTypeParams = struct()

    for key in defaultParams.keys():
        # if the field is part of the optional parameters (i.e. explicitly provided) use it.
        if hasattr(optParamStruct, key):
            setattr(problemTypeParams, key, optParamStruct[key])
            # and remove the field from the struct for the solver specific parameters.
            delattr(optParamStruct, key)
        else:
            # otherwise use the default parameter
            setattr(problemTypeParams, key, defaultParams[key])

    # assign all remaining parameters to the solver parameter struct.
    solverParams = optParamStruct
    
    return problemTypeParams, solverParams

@function
def getCobraSolverParams(solverType, paramNames, parameters, *args, **kwargs):
    # Daniel Ribeiro 2021: Python minimal implementation of the function getCobraSolverParams from Matlab COBRA Toolbox
    # This function gets the specified parameters in `paramNames` from
    # parameters, the global cobra paramters variable or default values set within
    # this script. It will use values with the following priority
    #
    # parameters > global parameters > default
    #
    # The specified parameters will be delt to the specified output arguements.
    # See examples below.
    #
    # USAGE:
    #
    #    varargout = getCobraSolverParams(solverType, paramNames, parameters)
    #
    # INPUTS:
    #    solverType:    Type of solver used: 'LP', 'MILP', 'QP', 'MIQP'
    #                   Daniel 2021: Python str
    #    paramNames:    Cell array of strings containing parameter names OR one
    #                   parameter name as string
    #                   Daniel 2021: Python list
    #
    # OPTIONAL INPUTS:
    #    parameters:    Structure with fields pertaining to parameter values that
    #                   should be used in place of global or default parameters.
    #                   parameters can be set to 'default' to use the default
    #                   values set within this script.
    #                   Daniel 2021: Python dict with paramNames: vals or string 'default' to use default values                  
    #
    # OUTPUTS:
    #    varargout:     Variables which each value corresponding to paramNames
    #                   is outputted to.
    #                   Daniel 2021: Python dictionary
    #
    # EXAMPLE:
    #    parameters.saveInput = 'LPproblem.mat';
    #    parameters.printLevel = 1;
    #    [printLevel, saveInput] = getCobraSolverParams('LP', {'printLevel', 'saveInput'}, parameters);
    #
    #    #Example using default values
    #    [printLevel, saveInput] = getCobraSolverParams('LP', {'printLevel','saveInput'}, 'default');
    #
    # .. Authors:
    #       - Richard Que (12/01/2009)
    #       - Ronan (16/07/2013) default MPS parameters are no longer global variables
    varargin = getCobraSolverParams.varargin
    nargin = getCobraSolverParams.nargin

    if nargin < 2:
        error('getCobraSolverParams: No parameters specified')

    if nargin < 3:
        parameters = ''

    # Persistence will make names specific to one type of solver.
    # Default Values
    # For descriptions of the different settings please have a look at 
    # getCobraSolverParamsOptionsForType
    #Daniel Ribeiro 2021: Use dict
    valDef = {
        'minNorm': None,
        'objTol': 1e-6,     # this should be used only when comparing the values of two objective functions
        'optTol': 1e-6,     # (dual) optimality tolerance
        'feasTol': 1e-6,    # (primal) feasibility tolerance
        'printLevel': 0,
        'verify': 0,
        'primalOnly': 0,
        'timeLimit': 1e36,
        'iterationLimit': 1000,
        'logFile': ''.join(['Cobra', solverType, 'Solver.log']),
        'saveInput': None,
        'PbName': ''.join([solverType, 'problem']),
        'debug': 0,
        'lifting': 0,

        'method': -1,
        # CPLEX parameters
        'DATACHECK': 1,
        'DEPIND': 1,
        'checkNaN': 0,
        'warning': 0,
        # tolerances
        'intTol': 1e-12,
        'relMipGapTol': 1e-12,
        'absMipGapTol': 1e-12,
        'NUMERICALEMPHASIS': 1
    }

    # Daniel Ribeiro 2021: Its a python list
    #if ~iscell(paramNames):
    #    paramNames = {paramNames}

    if solverType == 'LP':
        global CBT_LP_PARAMS
        parametersGlobal = CBT_LP_PARAMS
    elif solverType == 'QP':
        global CBT_QP_PARAMS
        parametersGlobal = CBT_QP_PARAMS
    elif solverType == 'MILP':
        global CBT_MILP_PARAMS
        parametersGlobal = CBT_MILP_PARAMS
    elif solverType == 'MIQP':
        global CBT_MIQP_PARAMS
        parametersGlobal = CBT_MIQP_PARAMS
    elif solverType == 'NLP':
        global CBT_NLP_PARAMS
        parametersGlobal = CBT_NLP_PARAMS
    else:
        fprintf('Unrecognized solver type')
        return

    # Initialize vararg out with None
    varargout = {}
    for p in paramNames:
        # set values to default
        if p in valDef.keys():
            varargout[p] = valDef[p]

        if not strcmp(parameters, 'default'): # skip if using default values
            # set values to global values
            if p in parametersGlobal.keys():
                varargout[p] = parametersGlobal[p]

            # set values to specified values
            if hasattr(parameters, p):
                varargout[p] = getattr(parameters, p)

    return varargout

@function
def getCobraSolverParamsOptionsForType(solverType, *args, **kwargs):
    # Daniel Ribeiro 2021: Python minimal implementation of the function getCobraSolverParamsOptionsForType from Matlab COBRA Toolbox
    # This function returns the available optional parameters for the specified
    # solver type.
    #
    # USAGE:
    #    paramnames = getCobraSolverParamsOptionsForType(solverType)
    #
    # INPUT:
    #    solverType:        One of the solver types available in the cobra
    #                       Toolbox ('LP','QP','MILP','MIQP','NLP')
    # OUPTUT:
    #    paramNames:        The possible parameters that can be set for the
    #                       given solver Type (depends on the solver Type
    varargin = getCobraSolverParamsOptionsForType.varargin
    nargin = getCobraSolverParamsOptionsForType.nargin
    
    if solverType == 'LP':
        paramNames = {'verify':None,         # verify that it is a suitable  LP problem
                      'minNorm':None,        # type of normalization used.
                      'printLevel':None,     # print Level
                      'primalOnly':None,     # only solve for primal
                      'saveInput':None,      # save the input to a file (specified)
                      'feasTol':None,        # feasibility tolerance
                      'optTol':None,         # optimality tolerance
                      'solver':None,         # solver to use (overriding set solver)
                      'debug':None,          # run debgugging code
                      'logFile':None,        # file (location) to write logs to
                      'lifting':None,        # whether to lift a problem
                      'method':None}              # solver method: -1 = automatic, 0 = primal simplex, 1 = dual simplex, 2 = barrier, 3 = concurrent, 4 = deterministic concurrent, 5 = Network Solver(if supported by the solver)

    elif solverType == 'QP':
        paramNames = {'verify':None,          # verify that it is a suitable  QP problem
                      'method':None,          # solver method: -1 = automatic, 0 = primal simplex, 1 = dual simplex, 2 = barrier, 3 = concurrent, 4 = deterministic concurrent, 5 = Network Solver(if supported by the solver)
                      'printLevel':None,      # print level
                      'saveInput':None,       # save the input to a file (specified)
                      'debug':None,           # run debgugging code
                      'feasTol':None,         # feasibility tolerance
                      'optTol':None,          # optimality tolerance
                      'logFile':None,         # file (location) to write logs to
                      'solver':None}          # the solver to use


    elif solverType == 'MILP':
        paramNames = {'intTol':None,          # integer tolerance (accepted derivation from integer numbers)
                      'relMipGapTol':None,    # relative MIP Gap tolerance
                      'absMipGapTol':None,    # absolute MIP Gap tolerance
                      'timeLimit':None,       # maximum time before stopping computation (if supported by the solver)
                      'logFile':None,         # file (location) to write logs to
                      'printLevel':None,      # print level
                      'saveInput':None,       # save the input to a file (specified)
                      'feasTol':None,         # feasibility tolerance
                      'optTol':None,          # optimality tolerance
                      'solver':None,          # solver to use (overriding set solver)
                      'debug':None}           # run debgugging code

    elif solverType == 'MIQP':
        paramNames = {'timeLimit':None,       # maximum time before stopping computation (if supported by the solver)
                      'method':None,          # solver method: -1 = automatic, 0 = primal simplex, 1 = dual simplex, 2 = barrier, 3 = concurrent, 4 = deterministic concurrent, 5 = Network Solver(if supported by the solver)
                      'feasTol':None,         # feasibility tolerance
                      'optTol':None,          # optimality tolerance
                      'intTol':None,          # integer tolerance (accepted derivation from integer numbers)
                      'relMipGapTol':None,    # relative MIP Gap tolerance
                      'absMipGapTol':None,    # absolute MIP Gap tolerance
                      'printLevel':None,      # print level
                      'saveInput':None,       # save the input to a file (specified)
                      'logFile':None,         # file (location) to write logs to
                      'solver':None}          # the solver to use

    elif solverType == 'NLP':
        paramNames = {'warning':None,         # whether to display warnings
                      'checkNaN':None,        # check for NaN solutions
                      'PbName':None,          # name of the problem
                      'iterationLimit':None,  # maximum number of iterations before stopping computation (if supported by the solver)
                      'timeLimit':None,       # time limit for the calculation
                      'logFile':None,         # file (location) to write logs to
                      'printLevel':None,      # print level
                      'saveInput':None,       # save the input to a file (specified)
                      'solver':None}          # the solver to use
    else:
        error(f'Solver type {solverType} is not supported by the Toolbox')

    return paramNames