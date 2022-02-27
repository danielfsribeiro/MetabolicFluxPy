# Generated with SMOP  0.41-beta
# Modified by Daniel Ribeiro 2020
from numpy.core.arrayprint import format_float_positional
from numpy.core.defchararray import asarray
from libsmop import *

#COBRA - cobrapy
import cobra
from solveCobraLP import solveCobraLP

#@function
#def call_Lee12(model=None,gene_names=None,gene_exp=None,gene_exp_sd=None,gene_scale_rxn=None,
#                flux_scale_rxn=None,flux=None,*args,**kwargs):
#    varargin = call_Lee12.varargin
#    nargin = call_Lee12.nargin

# Calls the implementation of Lee-12 provided in [Lee et al, BMC Sys Bio, 2012].
#
# Note: fixes to the original code are tagged with FIX

# INPUTS
#       model - cobra model
#       gene_names - genes ids
#       gene_exp - gene expression
#       gene_exp_sd - gene expression std
#       gene_scale_rxn - reaction to scale gene expression
#       flux_scale_rxn - reaction to scale flux distribution
#       flux - flux of scaling reaction

# OUTPUTS
#       fluxes - flux distribution

# Author: Daniel Machado, 2013

#Daniel Ribeiro: Adapt function call to COBRApy
global MODEL_ID
@function
def call_Lee12(model=None,
                gene_names=None,
                gene_exp=None,
                gene_exp_sd=None,
                flux=None,
                model_id= None,
                *args, **kwargs):
    varargin = call_Lee12.varargin
    nargin = call_Lee12.nargin
# INPUTS
#       model - cobra model
#       gene_names - genes ids (cobra gene format)
#       gene_exp - gene expression
#       gene_exp_sd - gene expression std
#       gene_names_cobra - dict of key(gene names) mapped to value(Cobra Genes)
#       gene_scale_rxn - reaction to scale gene expression
#       flux_scale_rxn - reaction to scale flux distribution
#       flux - flux of scaling reaction (the measured experimental fluxes)
#       model_id - name for the model
# OUTPUTS
#       fluxes - flux distribution

    global MODEL_ID
    MODEL_ID = model_id
    #Daniel Ribeiro: Change idx finding method to cobrapy
    unmeasured = setdiff(model.genes, gene_names)
    gene_names = concat([[gene_names],[unmeasured]], axis=None)
    gene_exp = concat([[gene_exp],[zeros(length(unmeasured),1)]])
    if isempty(gene_exp_sd):
        gene_exp_sd = zeros(length(gene_exp))
        disp('warning: no gene expression std given')
    else:
        gene_exp_sd = concat([[gene_exp_sd],[zeros(length(unmeasured),1)]])
    print("\nMapping gene to reactions...")
    rxn_exp, rxn_exp_sd = geneToReaction(model, gene_names, gene_exp, gene_exp_sd, nargout=2)
    rxn_exp = np.asarray(rxn_exp.flatten(order='F').tolist(), order='F')
    rxn_exp_sd = np.asarray(rxn_exp_sd.flatten(order='F').tolist(), order='F')

# Reaction to scale gene expression (comes from call_lee12.m form Machado, D. et al., 2014)
# Granata et al., did not perform these normalizations - comment out
    #if rxn_exp[scale_gene_idx] == 0:
        #error('Error: gene expression for scaling enzyme is zero.')
    #Scale all gene expression to GLUT1 gene (normalizing step)
    #scale_exp = rxn_exp[scale_gene_idx]
    #rxn_exp = rxn_exp / scale_exp
    #if any(rxn_exp_sd > 0):
    #    rxn_exp_sd[rxn_exp_sd == 0] = min(rxn_exp_sd[rxn_exp_sd > 0]) / 2
    #    rxn_exp_sd = rxn_exp_sd / scale_exp
    # FIX: scale all bounds accordingly, not just glc uptake

    # Daniel 2021: Use the exprimental fluxes as constriants instead of scaling all bounds (as Granata et al. did)
    # Constrain model according to experimental fluxes
    for key, val in flux.items():
        reaction_name = key
        flux = val[1]
        model.reactions.get_by_id(reaction_name).lower_bound = flux
        model.reactions.get_by_id(reaction_name).upper_bound = flux

    print("\nPerform flux analysis...")
    # Gene expression constraint FBA
    fluxes = dataToFlux(model,rxn_exp,rxn_exp_sd)

    return fluxes

if __name__ == '__main__':
    pass


@function
def geneToReaction(m=None,g=None,t=None,t_sd=None,*args,**kwargs):
    varargin = geneToReaction.varargin
    nargin = geneToReaction.nargin

    # kieran: 16 sep 11
    # Daniel Ribeiro 2020:
    # This function will recreate the reaction tha can be modeled in our system, since not all genes are expressed
    # m - cobra model
    # g - gene_names
    # t - expression
    # t_sd - expression sd

    # Daniel Ribeiro 2020: adapt to cobrapy
    r = zeros(*size(m.reactions))
    r_sd = zeros(*size(m.reactions))
    # Substitute nay gene_name containing '-' bt '_' in DictList
    for k in arange(0, length(g) - 1).reshape(-1):
        g[k] = strrep(g[k], '-', '_')
        #print(k, g[k])
    # Daniel Ribeiro 2020: gene_reaction_rule is accessible only by reaction
    # (in matlab you can access it from the model obj)
    # adapt code to call gene_reaction_rule from reaction obj
    print("Process GPRs for the model...")
    print("\tMatch [norm_counts]pm[SD] for each gene in GPR")
    print("\tSolve the GPR 'AND/OR' logic for each GPR")
    for k in range(len(m.reactions)):
        ga = m.reactions[k].gene_reaction_rule
        ga = strrep(ga, '-', '_')
        # Substitute [and or] by [AND, OR]. 
        # Later eval() will be used and the lowercase tokens are considered operators in Python
        ga = strrep(ga, 'and', 'AND')
        ga = strrep(ga, 'or', 'OR')
        # findall in Python is equivalent to match in Matlab
        # Get all gene names associated with GPR
        w = regexp(ga, '\\b\w*\\b', 'findall')
        # Clear empty strings after regex
        #print("w before '' clearing: ", w)
        indexes = []
        for i in arange(0, length(w) - 1).reshape(-1):
            if w[i] == '':
                indexes.append(i)
        w = cellarray(np.delete(w, indexes))
        # Collect only gene)names (discard AND, OR)
        w = setdiff(np.asarray(w), np.asarray(['AND','OR']))
        #Get expression data for each gene in w array
        for kk in arange(0, length(w) - 1).reshape(-1):
            # Convert Dictlist to numpy arrays
            # Get the index of gene_name in w, from DictList g
            j = find(strcmp(w[kk], np.asarray(g, dtype = np.str_)))
            # Get expression from normCounts and stdev list
            n = t[j]
            n_sd = t_sd[j]
            # Daniel Ribeiro 2021: join str in Python format
            ga = regexprep(ga, 
                            ''.join(['\\b', w[kk], '\\b']),
                            ''.join(np.asarray([num2str(n), 'pm', num2str(n_sd)]).flatten(order='F').tolist())
                            )
        #print("\tSubstituted ApmB GPR: ", ga)

        n, n_sd = addGeneData(ga, nargout=2)
        r[k,0] = n
        r_sd[k,0] = n_sd
    return r, r_sd

if __name__ == '__main__':
    pass


@function
def AandB(str1=None,str2=None,*args,**kwargs):
    varargin = AandB.varargin
    nargin = AandB.nargin
    
    #FIX Daniel
    #Daniel Ribeiro 2021: '\.' '\-' is needed, as opposed to '.' '-'
    #ApmB='([0-9\.\-e])+pm([0-9\.\-e]+)'
    ApmB = '([0-9\.\-e]+)pm([0-9\.\-e]+)'
    match_expr = ApmB
    #Daniel Ribeiro 2021: adapt gegex for Python
    m1 = eval(regexprep(str1, match_expr, '\g<1>'))
    s1 = eval(regexprep(str1, match_expr, '\g<2>'))
    m2 = eval(regexprep(str2, match_expr, '\g<1>'))
    s2 = eval(regexprep(str2, match_expr, '\g<2>'))
    m,j = min(concat([m1,m2]), nargout=2)
    if j == 0:
        s = s1
    else:
        s = s2
    string = ''.join([num2str(m), 'pm', num2str(s), ' '])
    return string
    
if __name__ == '__main__':
    pass
    
    
@function
def AorB(str1=None,str2=None,*args,**kwargs):
    varargin = AorB.varargin
    nargin = AorB.nargin
    
    #FIX Daniel
    #Daniel Ribeiro: '\.' '\-' is needed, as opposed to '.' '-'
    #ApmB='([0-9\.\-e])+pm([0-9\.\-e]+)'
    ApmB = '([0-9\.\-e]+)pm([0-9\.\-e]+)'
    match_expr = ApmB
#Daniel Ribeiro 2021: adapt gegex for Python
    m1 = eval(regexprep(str1,match_expr,'\g<1>'))
    s1 = eval(regexprep(str1,match_expr,'\g<2>'))
    m2 = eval(regexprep(str2,match_expr,'\g<1>'))
    s2 = eval(regexprep(str2,match_expr,'\g<2>'))
    m = m1 + m2
    s = sqrt(s1 ** 2 + s2 ** 2)
    string = ''.join([num2str(m),'pm',num2str(s), ' '])
    return string
    
if __name__ == '__main__':
    pass
    
    
@function
def addGeneData(g=None,*args,**kwargs):
    varargin = addGeneData.varargin # args tuple (cell array)
    nargin = addGeneData.nargin     # len args tuple
    #g is arg=0, access by g[0]

    # kieran: 22 july 11
    
    n = np.nan
    n_sd = np.nan
    
    #FIX Daniel
    # Daniel Ribeiro: '\.' '\-' is needed, as opposed to '.' '-'
    #ApmB='([0-9\.\-e])+pm([0-9\.\-e]+)'
    ApmB = '([0-9\.\-e]+)pm([0-9\.\-e]+)'
    tries = 0
    # Lambda function: Adapted for Python regex, where function can only take match object as argument
    f_and = lambda match_obj=None: AandB(match_obj.group(1), match_obj.group(4))
    # Lambda function:
    f_or = lambda match_obj=None: AorB(match_obj.group(1), match_obj.group(4))
    if logical_not(isempty(g) or g == ''):
        while isnan(n):
            tries = tries + 1
            if tries > 10000:
                fprintf(1,'\tWarning: stuck at loop evaluating... %s\n', g)
                break
            #Daniel Ribeiro 2021: eval() cannot evaluate empty strings or invalid tokens, exception is raised
            try:
                match_expr = ApmB
                #Daniel Ribeiro 2021: adapt to Python regex
                g_av = regexprep(g, match_expr, '\g<1>')
                g_sd = regexprep(g, match_expr, '\g<2>')
                n = eval(g_av)
                n_sd = eval(g_sd)
            #TODO there is code missing for exception hanling. Check call_Lee12.m
            except:
                # Daniel Ribeiro 2021: adpat match_expr for Python regex expression
                # replace parenthesis
                # each ApmB may have parenthesis of the style '( ApmB )'. Remove those
                #match_expr      = ['\(\s*(',ApmB,')\s*\)'];
                if '(' in g or ')' in g:
                    parenthesis = True
                else:
                    parenthesis = False
                match_expr      = '\(\s*(' + ApmB + ')\s*\)'
                replace_expr    = '\g<1>'
                g = regexprep(g, match_expr, replace_expr)
                # replace and
                #match_expr = ['(',ApmB,')\s+and\s+(',ApmB,')'];
                match_expr = '(' + ApmB + ')\s+AND\s+(' + ApmB + ')'
                # The match_expr reads: (([0-9\.\-e]+)pm([0-9\.\-e]+))\s+AND\s+(([0-9\.\-e]+)pm([0-9\.\-e]+))\s+
                # Groups are read form left to right, from outer-level to inner-lever before changing to the next group at the same level
                # (g<1>(g<2>)pm(g<3>)) AND (g<4>(g<5>)pm(g<6>))
                #replace_expr = '${AandB($1,$2)}';
                replace_expr = f_and     #Call AandB func
                g = regexprep(g, match_expr, replace_expr, 'once')
                # replace or
                #match_expr      = ['(',ApmB,')\s+or\s+(',ApmB,')']
                match_expr = '(' + ApmB + ')\s+OR\s+(' + ApmB + ')'
                replace_expr = f_or     #Call AorB func
                g = regexprep(g, match_expr, replace_expr, 'once')

    #print(f"\tAfter all AND OR replacements: {g}. Took {tries} tries.")
    return n, n_sd
    
if __name__ == '__main__':
    pass
    
    
@function
def dataToFlux(m=None,r=None,r_sd=None,*args,**kwargs):
    varargin = dataToFlux.varargin
    nargin = dataToFlux.nargin

    # kieran: 21 sep 11
    
    # Daniel Ribeiro 2021: Adapt to Cobrapy
    print("Working on flux analysis...")
    rev = np.zeros((len(m.reactions),1), dtype=bool).reshape(-1)
    
    nR_old = 0 # number of reversible reactions in previous state
    v_sol = np.zeros(len(m.reactions))
    
    tries = 0
    # rev – 0 = irrev, 1 = rev #Reversibility of reactions
    while sum(logical_not(m.reactions.list_attr('reversibility'))) > nR_old:
        tries = tries + 1
        print(f"Performing {tries}/100 tries...\tirreversible old {nR_old}, irreversible new {sum(logical_not(m.reactions.list_attr('reversibility')))}")
        if tries > 100:
            fprintf(1, 'warning: stuck at loop. irrev old %d, irrev new %d \n', nR_old, sum(logical_not(m.reactions.list_attr('reversibility'))))
            break

        nR_old = sum(logical_not(m.reactions.list_attr('reversibility')))
        ### Step 1 - fit to data
        N = cobra.util.array.create_stoichiometric_matrix(m)
        L = np.asarray(m.reactions.list_attr('lower_bound'))
        U = np.asarray(m.reactions.list_attr('upper_bound'))
        f = np.zeros(len(m.reactions)).T
        b = np.zeros(len(m.metabolites))
        # numpy array resize and copy takes too long
        # Resize to a pre allocated matrix
        # Find initial array size
        s1, s2 = size(N, nargout=2)
        new_rows = 0
        new_cols = 0
        for k in arange(0, length(m.reactions)-1).reshape(-1):
            d = r[k]
            s = r_sd[k]
            if logical_not(m.reactions[k].reversibility) and logical_not(isnan(d)) and s > 0: 
                new_rows = new_rows + 1
                new_cols = new_cols + 2
        new_rows = new_rows + s1
        new_cols = new_cols + s2
        # Copy N array into N_temp
        N_temp = np.zeros((new_rows, new_cols))
        N_temp[0:N.shape[0], 0:N.shape[1]] = N
        
        # keep track of the N_temp row, col expansion
        s1, s2 = size(N, nargout=2)
        for k in arange(0, length(m.reactions)-1).reshape(-1):
            d = r[k]
            s = r_sd[k]
            if logical_not(m.reactions[k].reversibility) and logical_not(isnan(d)) and s > 0:  
                # In MATLAB, you can add one or more elements to a matrix by placing them outside of the existing row and column index boundaries. 
                # The matrix is padded with zeros
                # However this is not possible with Python so easily (index error raised). Adapt code to reflect this behaviour.
                N_temp[s1-1 + 1, k] = 1
                N_temp[s1-1 + 1, s2-1 + 1] = -1
                N_temp[s1-1 + 1, s2-1 + 2] = 1
                col_expand = np.zeros(2)
                L = np.concatenate((L, col_expand))
                L[s2-1 + 1] = 0
                L[s2-1 + 2] = 0
                col_expand = np.zeros(2)
                U = np.concatenate((U, col_expand))
                U[s2-1 + 1] = np.inf
                U[s2-1 + 2] = np.inf
                col_expand = np.zeros(1)
                b = np.concatenate((b, col_expand))
                b[s1-1 + 1] = d
                col_expand = np.zeros(2)
                f = np.concatenate((f, col_expand))
                f[s2-1 + 1] = -1 / s
                f[s2-1 + 2] = -1 / s
                #Update positions
                s1 = s1 + 1
                s2 = s2 + 2
        N = N_temp  # copy N_temp back to N

        v, fOpt, conv = easyLP(f, N, b, L, U, m, nargout=3)
        if conv:
            v_sol = v[0: length(m.reactions)].copy()
            for k in arange(0,length(m.reactions)-1).reshape(-1):
                if rev[k]:
                    v_sol[k] = -v_sol[k]
           
            ### 2. run FVA           
            N = np.concatenate((N, f.reshape(1,-1)), axis=0) # add f row
            b= np.concatenate((b,np.asarray(fOpt).reshape(-1))) # extent b with fOpt (b size should match N.row_size)
            for k in arange(0, length(m.reactions)-1).reshape(-1):
                if m.reactions[k].reversibility:
                    f = np.zeros(L.size).reshape(-1)
                    f[k] = -1
                    __,fOpt,conv=easyLP(f,N,b,L,U,m,nargout=3)
                    print(f"\tFor k={k} -> fOpt={fOpt}, conv={conv}")
                    if conv and (-fOpt >= 0):   # irreversibly forward
                        print(f"\t-fOpt >=0: {-fOpt >= 0}, irreversibly forward")
                        
                        m.reactions[k].lower_bound = builtins.max(m.reactions[k].lower_bound, 0)
                    else:
                        f[k]=1
                        __,fOpt,conv=easyLP(f,N,b,L,U,m,nargout=3)
                        if conv and builtins.abs(fOpt) <= 0:    # irreversibly backward
                            # cobrapy does not seem to let you edit the S matrix as CobraToolbox does
                            # The workaround is to iterate over reactions and metabolites specifically, and manually set values
                            # Look at cobra.util.array.create_stoichiometric_matrix for hints on how to do this
                            # Here we want to multiply by -1 all values in the reaction (reverses the reaction) 
                            #m.S[:,k] = -m.S[arange(),k]    #original matlab code
                            print(f"\tabs(fOpt) <=0: {builtins.abs(fOpt) <= 0}, irreversibly backward")
                            metabolites = m.reactions[k].metabolites # Dictlist of metabolites and their stoichiometric coeficients (a copy)
                            for metab in metabolites:
                                #Multiply coefficient by -1
                                metabolites[metab] *= -1
                            # Substiture coeficients. combine=False ensures that coefficient is replaced
                            m.reactions[k].add_metabolites(metabolites, combine=False)
                            m.reactions[k].upper_bound = -m.reactions[k].upper_bound
                            m.reactions[k].lower_bound = -m.reactions[k].lower_bound
                            ub = m.reactions[k].upper_bound
                            m.reactions[k].upper_bound = m.reactions[k].lower_bound
                            m.reactions[k].lower_bound = ub
                            m.reactions[k].lower_bound = builtins.max(m.reactions[k].lower_bound, 0)
                            #m.reactions[k].reversibility = 0    # This does nothing and issues a warning in cobrapy
                            rev[k] = logical_not(rev[k])
    return v_sol
    
if __name__ == '__main__':
    pass
    
    
@function
def easyLP(f=None,a=None,b=None,vlb=None,vub=None, model=None,*args,**kwargs):
    varargin = easyLP.varargin
    nargin = easyLP.nargin
# easyLP
    
# solves the linear programming problem: 
#   max f'x subject to 
#   a x = b
#   vlb <= x <= vub.
    
# Usage: [v,fOpt,conv] = easyLP(f,a,b,vlb,vub)
    
#   f           objective coefficient vector (f,c)
#   a           LHS matrix (N,A)
#   b           RHS vector (b)
#   vlb         lower bound vector (L)
#   vub         upper bound vector (U)
    
#   v           solution vector
#   fOpt        objective value
#   conv        convergence of algorithm [0/1]
    
# the function is a wrapper for on the "solveCobraLP" script provided with
# the COBRA (COnstraint-Based Reconstruction and Analysis) toolbox 
# http://opencobra.s.net/
    
    #kieran, 20 april 2010
    
    # matlab can crash if inputs nan
    if any(isnan(f)) or any(any(isnan(a))) or any(isnan(b)) or any(isnan(vlb)) or any(isnan(vub)):
        error('nan inputs not allowed')
    
    # initialize
    v = np.zeros(vlb.shape)
    v = ravel(v).reshape(-1)
    f = ravel(f).reshape(-1)
    vlb = ravel(vlb).reshape(-1)
    vub = ravel(vub).reshape(-1)

    # remove any tight contstraints as some solvers require volume > 0
    j1 = (vlb != vub)
    j2 = (vlb == vub)
    v[j2] = vlb[j2]
    b = (ravel(b).reshape(-1) - a @ v).reshape(-1)
    #a(:,j2) = []; remove columns where j2 is true 
    # logically this is the same as keeping the columns where (not j2) is true
    a = a[:,np.logical_not(j2)]
    vlb = vlb[np.logical_not(j2)]
    vub = vub[np.logical_not(j2)]
    f0 = f
    f = f[np.logical_not(j2)]
    fOpt = np.nan
    csense = np.full(b.shape, 'E').reshape(-1)  #constraint sense. Default is 'E' equality
    solution = solveCobraLP(struct('A',a, 'b',b, 'c',f, 'lb',vlb, 'ub',vub, 'osense',-1, 'csense',csense, 'modelID',f'model_{MODEL_ID}'), model)
    conv = (solution.stat == 1)
    if conv:
        v0 = solution.full
        v[j1] = v0
        fOpt = f0.T @ v
    
    return v, fOpt, conv
    
if __name__ == '__main__':
    pass
    