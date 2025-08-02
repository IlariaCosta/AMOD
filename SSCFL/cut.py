
from typing import Tuple
import numpy as np
import fractions
import logging
import cplex #####
import io

columns=["name", "cluster_type", "nvar","nconstraints","optimal_sol","sol","sol_is_integer","status","ncuts","elapsed_time","gap","relative_gap","iterations"]

logging.basicConfig(filename='resolution.log', format='%(asctime)s - %(message)s',level=logging.INFO, datefmt='%d-%b-%y %H:%M:%S')

def getProblemData(f, c_matrix, demands, capacity) -> Tuple:
    """
    Costruisce vettore c, matrice A e vettore b per la formulazione LP del problema UFL.

    Args:
        f: array costi apertura facilities (shape m,)
        c_matrix: matrice costi assegnamento (shape m x n)

    Returns:
        c: vettore funzione obiettivo (shape m + m*n,)
        A: matrice vincoli (shape (n + m*n) x (m + m*n))
        b: vettore termini noti (shape n + m*n,)
    """
    m = len(f)             # numero facilities
    c_matrix = np.transpose(c_matrix)
    c_matriX = np.array(c_matrix)
    n = c_matrix.shape[1]  # numero clienti
    print("numero clienti ->",n,"; numero facilities ->", m)
    # calcolo vettore coefficienti funzione obiettivo
    # [f_1..f_m, c_11, c_12, ..., c_mn]
    c = np.concatenate((f, c_matrix.flatten()))
    
    # Numero variabili = m + m*n (+ slack (facility*cliente + facility))
    num_vars = m + m * n + (m * n + m)
    print(f"numero variabili = {num_vars}")
   
    # Numero vincoli = n + m*n + m
    num_constraints = n + m * n + m
    print(f"numero vincoli = {num_constraints}")
    # inizializza matrice di zeri
    A = np.zeros((num_constraints, num_vars))
    b = np.zeros(num_constraints)
    #print("\n qui")
    # Cliente j assegnato ad una sola facility
    # le prime n righe descrivono i vincoli
    for j in range(n):
        for i in range(m):
            A[j][m + i * n + j] = 1
        b[j] = 1
    #print("\n qui")
    # Vincoli x_ij <= y_i (equivalente a -y_i + x_ij <= 0)
    # da riga n a riga n + m*n -1
    # per ogni i,j
    for i in range(m):
        for j in range(n):
            row_idx = n + i * n + j
            y_idx = i
            x_idx = m + i * n + j
            s_idx = m * n + m + (i * n + j)
            A[row_idx, y_idx] = -1
            A[row_idx, x_idx] = 1
            A[row_idx, s_idx] = 1
            b[row_idx] = 0
            #print(row_idx, y_idx,x_idx,s_idx)
    #print("\n qui")
    # d[j] * x[j,i] - capacity[i] * y[i] + s_c[i]<= 0
    # da riga n+n*m fino a n + n*m + m
    start_capacity_row = n + (n*m)
    for i in range(m):
        row_idx = start_capacity_row +i
        y_idx = i
        A[row_idx][y_idx] = -capacity[i]
        for j in range(n):
            x_idx = m + i * n + j
            s_c_idx = (2* m * n + m) + (i)
            #print(f"indice y = {y_idx}\n indice x = {x_idx} \n indice s = {s_idx} ====> riga {row_idx}")            
            A[row_idx][x_idx] = demands[j]
            A[row_idx][s_c_idx] = 1
        b[row_idx] = 0
    print("dimensioni matrice A: ", len(A), len(A[0]))
    
    return c, A, b
    
   
def initializeInstanceVariables(n,m) : 
    # n = numero clienti
    # m = numero facilities
    names = []
    lower_bounds = []
    upper_bounds = []
    constraint_names = []
    constraint_senses = []

    # Variables y
    for i in range(m):        
        names.append("y"+str(i))
        lower_bounds.append(0.0)
        upper_bounds.append(1.0)
    #print("lunghezza nomi dopo y: ", len(names))
    # variables x
    for i in range(m):
        for j in range(n):
            names.append("x"+str(i)+str(j))
            lower_bounds.append(0.0)
            upper_bounds.append(1.0)
    #print("lunghezza nomi dopo x: ", len(names))
    
    # variables s
    for i in range(m):
        for j in range(n):
            names.append("s"+str(i)+str(j))
            lower_bounds.append(0.0)
    #print("lunghezza nomi dopo s: ", len(names))

    # variablile c_cap
    for i in range(m):
        names.append("s_c"+str(i))
        lower_bounds.append(0.0)
    #print("lunghezza nomi dopo s_c: ", len(names))
        
    # Vincoli
    for i in range(n + m * n + m):
        constraint_names.append("c"+str(i))
        constraint_senses.append("E")
    
    return names, lower_bounds, upper_bounds,constraint_senses,constraint_names
        

def get_tableau(prob,A,b):
    '''
    This function get the final tableau of the prob (cplex.Cplex())
    
    Arguments:
        problem -- cplex.Cplex()
     
    returns:
        n_cuts 
        b_bar 
    '''
   # print(help(prob.solution))
    #print(prob.solution.get_indicator_slack())
    b_bar = np.zeros(len(b))
    for i in range(len(b)):
        # Ottieni la riga i-esima di B^-1
        binv_row = prob.solution.advanced.binvrow(i)
        b_bar[i] = np.dot(binv_row, b)

    
    # sol_status = prob.solution.get_status_string()
    # print(sol_status)
    # if sol_status in [prob.solution.status.optimal,prob.solution.status.optimal_tolerance]:
    # # La soluzione è valida, quindi puoi chiedere lo stato delle variabili
    #     col_status, row_status = prob.solution.basis.get_status()
    
    # basic_indices = [i for i, status in enumerate(col_status) if status == prob.basis.status.basic]
    # nonbasic_indices = [i for i, status in enumerate(col_status) if status != prob.basis.status.basic]
    
    # for i, basic_var_index in enumerate(basic_indices):
    #     binv_row = prob.solution.advanced.binvrow(i) # Questo è corretto
    #     b_bar_i = np.dot(binv_row, b)
    #     if abs(b_bar_i - round(b_bar_i)) > 1e-6:
    #         print("valore frazionario")
    # print(help(prob.solution.advanced))
    # print("calcolo B inversa")
    # BinvA = np.array(prob.solution.advanced.binvarow())
    # print("calcolata B inversa")
    # nrow = BinvA.shape[0]
    # ncol = BinvA.shape[1]
    # nrow = prob.linear_constraints.get_num()
    # ncol = prob.variables.get_num()
    try:
        nrow = prob.linear_constraints.get_num()
        ncol = prob.variables.get_num()
        #print(f"numero colonne {ncol}")
    except Exception as e:
        print("Errore durante l'accesso a prob:", e)
    #print(f"numero colonne {ncol}, numero righe {nrow}")
    b_bar = np.zeros(nrow)
    
    varnames = prob.variables.get_names()
    b = prob.linear_constraints.get_rhs()
    #print("b = ", b)
    mat_Binv = prob.solution.advanced.binvrow()
    Binv = np.array(mat_Binv)
    b_bar = np.matmul(Binv, b)
    #print("b_bar = ", b_bar)
    idx = 0     # Compute the nonzeros
    n_cuts = 0  # Number of fractional variables (cuts to be generated)
    #print('\n\t LP relaxation final tableau:\n')
    # Binv_A = prob.solution.advanced.binvarow() 
    
    for i in range(nrow):
        output_t = io.StringIO()
        z = prob.solution.advanced.binvarow(i)
        # z = Binv_A[i]
        # print("popolo z numero -> ", i)
        # print("prima")
        # print(ncol)
        for j in range(ncol):
            if z[j] > 0:
                print('+', end='',file=output_t)
            # print(z[j])
            zj = fractions.Fraction(z[j]).limit_denominator(1000)
            num = zj.numerator
            den = zj.denominator
            if num != 0 and num != den:
                print(f'{num}/{den} {varnames[j]} ', end='',file=output_t)
            elif num == den:
                print(f'{varnames[j]} ', end='',file=output_t)
                # print(z[j])
            val = z[j]
            # printf'z[{j}] = {val}')  # 👈 AGGIUNGI QUESTO
            if abs(val - round(val)) > 1e-6: 
                #print(f'z[{j}] = {val}')
                zj = fractions.Fraction(z[j]).limit_denominator(1000)
                num = zj.numerator
                den = zj.denominator
                # print(num, den)
                if num != 0 and num != den:
                    print(f'{num}/{den} {varnames[j]} ', end='',file=output_t)
                else :
                    print(f'{varnames[j]} ', end='',file=output_t)
            else:
                print(f'{varnames[j]} ', end='', file=output_t)
            if np.floor(z[j]+0.5) != 0:
                idx += 1
        b_bar_i = fractions.Fraction(b_bar[i]).limit_denominator()
        
        num = b_bar_i.numerator
        den = b_bar_i.denominator
        # print(n_cuts)
        print(f'= {num}/{den}',file=output_t)
        # print("z popolato")
        # print("prima")
        contents = output_t.getvalue()
        logging.info("%s",contents)
        output_t.close()

        # Count the number of cuts to be generated
        #print(f"DEBUG: b_bar[{i}] = {b_bar[i]}, floor(b_bar[{i}]) = {np.floor(b_bar[i])}")
        
        # Count the number of cuts to be generated
        if np.floor(b_bar[i]) != b_bar[i]:
            n_cuts += 1    
    logging.info("Cuts to generate: %d", n_cuts)
    return n_cuts , b_bar


# def determineOptimal(instance, cluster_type):
#     '''
#     This function determines the optimal solution of the given instance.
    
#     Arguments:
#         instance
#     '''
#     c, A, b = getProblemData(instance) 
#     nCols, nRows = (len(c), len(b))
#     # Get the instance name
#     txtname = instance.split("/")[2]
#     name = txtname.split(".txt")[0]
#     cplexlog = name+".log"
#     #Program variables section ####################################################
#     names = []
#     all_constraints = []
#     constraint_names = []
#     constraint_senses = []
#     # Variables 
#     for i in range(nCols):
#         names.append("x"+str(i))
#     # Constraint 
#     for i in range(nRows):
#         constraint_names.append("c"+str(i))
#         constraint_senses.append("L")
#     with cplex.Cplex() as mkp:
#         mkp.set_problem_name(name)
#         mkp.objective.set_sense(mkp.objective.sense.maximize)
#         mkp.set_log_stream(None)
#         mkp.set_error_stream(None)
#         mkp.set_warning_stream(None)
#         mkp.set_results_stream(None)
#         params = mkp.parameters
#         # Disable presolve 
#         params.preprocessing.presolve.set(0) 
#         # Add variables & Slack --------------------------------------------------------------------
#         mkp.variables.add(names=names, types=[mkp.variables.type.binary] * nCols)
#         # Add contraints -------------------------------------------------------------------
#         for i in range(nRows):
#             mkp.linear_constraints.add(lin_expr= [cplex.SparsePair(ind= [j for j in range(nCols)], val= A[i])],
#              rhs= [b[i]], names = [constraint_names[i]], senses = [constraint_senses[i]])
#             all_constraints.append(A[i])
#         # Add objective function -----------------------------------------------------------
#         for i in range(nCols): 
#             mkp.objective.set_linear([(i, c[i])])
#         # Resolve the problem instance
#         mkp.solve()
#         # Report the results 
#         logging.info("\n\t\t\t\t\t\t*** OPTIMAL PLI SOLUTION ***")
#         print_solution(mkp)
#         mkp.write("lp/"+cluster_type+"/"+name+"/optimal.lp")
#         mkp.solution.write("solutions/"+cluster_type+"/"+name+"/optimal.log")
#         optimal_sol= mkp.solution.get_objective_value()
#     return optimal_sol


def initialize_fract_gc(n_cuts,ncol , prob, varnames, b_bar) : 
    '''
    
    Arguments:
        n_cuts
        ncol
        nrow
        prob
        varnames
        b_bar
    returns:
        gc_lhs
        gc_rhs 
    '''
    
    cuts = np.zeros([n_cuts,ncol])
    cut_limits= []
    gc_sense = [''] * n_cuts
    gc_rhs   = np.zeros(n_cuts)
    gc_lhs   = np.zeros([n_cuts, ncol])
    rmatbeg  = np.zeros(n_cuts)
    rmatind  = np.zeros(ncol)
    rmatval  = np.zeros(ncol)
    logging.info('Generating Gomory cuts...\n')
    cut = 0  #  Index of cut to be added
    for i in range(len(b_bar)):
        idx = 0
        output = io.StringIO()
        
        if np.floor(b_bar[i]) != b_bar[i]:
            print(f'Row {i+1} gives cut -> ', end = '', file=output)
            #print("sono entrato qui dentro")
            z = np.copy(prob.solution.advanced.binvarow(i)) # Use np.copy to avoid changing the
                                                        # optimal tableau in the problem instance
            rmatbeg[cut] = idx
            for j in range(ncol):
                z[j] = z[j] - np.floor(z[j]) #calcolo l parte frazionaria
             
                if z[j] != 0:
                    rmatind[idx] = j
                    rmatval[idx] = z[j]
                    idx +=1
                # Print the cut
                if z[j] > 0:
                    print('+', end = '',file=output)
                if (z[j] != 0):
                        fj = fractions.Fraction(z[j])
                        fj = fj.limit_denominator()
                        num, den = (fj.numerator, fj.denominator)
                        print(f'{num}/{den} {varnames[j]} ', end='',file=output)
           
            gc_lhs[cut,:] = z
            cuts[cut,:]= z
            gc_rhs[cut] = b_bar[i] - np.copy(np.floor(b_bar[i])) # np.copy as above
            #print(gc_rhs[cut])
            gc_sense[cut] = 'L'
            gc_rhs_i = fractions.Fraction(gc_rhs[cut]).limit_denominator()
            num = gc_rhs_i.numerator
            den = gc_rhs_i.denominator
            print(f'>= {num}/{den}',file=output)
            cut_limits.append(gc_rhs[cut])
            cut += 1
            contents = output.getvalue()
            output.close()
            logging.info(contents)
            #index +=1
            #print("\n\tgc_rhs = ", gc_rhs)
    return gc_lhs, gc_rhs   # lhs è la parte sinistra del taglio
                            # rhs è la parte destra del taglio

def generate_gc(mkp, A, gc_lhs, gc_rhs, names) : 
    
    logging.info('*** GOMORY CUTS ***\n')
    cuts = []
    cuts_limits = []
    cut_senses = []
    #print(len(gc_lhs))
    #print("\tQUI")
    for i in range(len(gc_lhs)):
        output = io.StringIO()
  
        current_gc_lhs = gc_lhs[i] # Coefficienti del taglio corrente
        #print(current_gc_lhs)
        current_gc_rhs = gc_rhs[i] # Termine noto del taglio corrente
        cuts.append([])
        
        #lhs, rhs = get_lhs_rhs(mkp, gc_lhs[i], gc_rhs[i], A)
        # Print the cut
        cut_string_parts = []
        for j in range(len(current_gc_lhs)): # Itera su tutte le variabili nel taglio
            coefficient = current_gc_lhs[j]
            fj = fractions.Fraction(coefficient).limit_denominator()
            if len(cut_string_parts) > 0 and fj > 0:
                cut_string_parts.append("+")
            if fj == 1:
                cut_string_parts.append(f"{names[j]}")
            elif fj == -1:
                cut_string_parts.append(f"- {names[j]}")
            else:
                cut_string_parts.append(f"{fj} {names[j]}")
            cuts[i].append(float(coefficient))


        cut_senses.append('G')
        cuts_limits.append(float(current_gc_rhs)) # Aggiunge il RHS corretto
        # Completa la stringa per la stampa
        final_cut_string = " ".join(cut_string_parts)
        print(f"{final_cut_string} >= {fractions.Fraction(current_gc_rhs).limit_denominator()}", file=output)

        
    return cuts, cuts_limits, cut_senses

def get_lhs_rhs(prob, cut_row, cut_rhs, A):

    ncol = len(A[0])
    cut_row = np.append(cut_row, cut_rhs)
    b = np.array(prob.linear_constraints.get_rhs())
    A = np.append(A, b.reshape(-1, 1), axis=1)
    plotted_vars = np.nonzero(prob.objective.get_linear())[0]
    # Assumption: plotted variables are at the beginning of the initial tableau
    for i, sk in enumerate(range(len(plotted_vars), ncol)):
        cut_coef = cut_row[sk]
        cut_row -= A[i,:] * cut_coef
    lhs = cut_row[:len(plotted_vars)]
    rhs = cut_row[ncol:]
    return lhs, rhs
    
def print_solution(prob):# : cplex.Cplex()):
    '''
    This function print solution of problem (cplex.Cplex())
    
    Arguments:
        problem -- cplex.Cplex()
    
    '''
    ncol = len(prob.variables.get_cols())
    nrow = len(prob.linear_constraints.get_rows())
    varnames = prob.variables.get_names()
    slack = np.round(prob.solution.get_linear_slacks(), 3)
    x = np.round(prob.solution.get_values(), 3)

    # Log everything about the solutions found
    logging.info("\t-> Solution status = %s", prob.solution.status[prob.solution.get_status()])
    logging.info("\t-> Solution value  = %f\n", prob.solution.get_objective_value())
    logging.info("SLACKS SITUATION:")
    for i in range(nrow):
        logging.info(f'-> Row {i}:  Slack = {slack[i]}')
    logging.info("\n\t\t\t\t\t PROBLEM VARIABLES:")
    for j in range(ncol):
        logging.info(f'-> Column {j} (variable {varnames[j]}):  Value = {x[j]}')
    
    sol=  prob.solution.get_objective_value()
    sol_type= sol.is_integer()
    status = prob.solution.status[prob.solution.get_status()]
    return sol, sol_type, status