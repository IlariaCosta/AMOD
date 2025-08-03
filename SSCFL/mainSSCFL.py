import time
import csv
import os
import cplex
from amplpy import AMPL, Environment
import numpy as np
from utils import (
    compute_gap,
    is_integral,
    parse_dat_file,
    build_A_b,
    invert_matrix,
    mat_vec_mul,
    mat_mul
)
from gomory import (
    solve_with_gomory
)
from cut import (
    getProblemData,
    get_tableau,
    initializeInstanceVariables,
    initialize_fract_gc,
    generate_gc,
    print_solution
) 
 


#*************************#
#   ---     UFL     ---   #
#*************************#
 
def run_sscfl_experiment(mod_path_int, mod_path_relax, data_path):
    # 1. CARICO MODELLO INTERO
    print("=================================================================")
    ampl = AMPL()
    ampl.set_option('solver', 'cplexamp')
    ampl.set_option('solver_msg', 0)
    ampl.read(mod_path_int)
    ampl.read_data(data_path)
    ampl.set_option('cplex_options', 'mipgap=0')
   
    t0 = time.time()
    ampl.solve() 
    time_int = time.time() - t0
    
    obj_int = ampl.obj['TotalCost'].value()
    print(f"\n\tSoluzione intera  = {obj_int}");

    # 2. CARICO MODELLO RILASSATO
    print("=================================================================")
    ampl_relax = AMPL()
    ampl_relax.set_option('solver', 'cplexamp')
    ampl_relax.set_option('solver_msg', 0)
    ampl_relax.read(mod_path_relax)
    ampl_relax.read_data(data_path)

    # Rilassa da Python tutte le variabili intere
    variables = ampl_relax.get_variables()
    for var in variables:
         try:
             if var.is_integer():
               var.set_integer(False)
         except:
            continue
    
    t0 = time.time()
    ampl_relax.solve()
    y_vals = ampl_relax.get_variable('y').get_values().to_dict()
    
    time_relax = time.time() - t0
    
    obj_relax = ampl_relax.obj['TotalCost'].value()
    print(f"\n\tSoluzione rilassata = {obj_relax:.4f}");
    
    gap_relax = compute_gap(obj_relax, obj_int)
    
    x_vals = list(ampl_relax.get_variable('y').get_values().to_dict().values())

    already_integer = is_integral(x_vals)

    if already_integer:
        print ("SOLUZIONE GIA' INTERA\n");
        return {
            "istanza": os.path.basename(data_path),
            "obj_int": obj_int,
            "time_int": time_int,
            "obj_relax": obj_relax,
            "time_relax": time_relax,
            "gap_relax": 0,
            "obj_gomory_all": obj_relax,
            "time_gomory_all": 0,
            "iter_gomory_all": 0,
            "gap_gomory_all": 0,
            "obj_gomory_step": obj_relax,
            "time_gomory_step": 0,
            "iter_gomory_step": 0,
            "gap_gomory_step": 0,
            "nota": "Relax già intero"
        }
    #3. Gomory un taglio alla volta
    print("=================================================================")
    print(f"|\t\t\t\t\t\t\t\t|\n|\t\t\tGOMORY CUTS\t\t\t\t|\n|\t\t\t\t\t\t\t\t|")
    print("=================================================================")
    print("|\t\tTAGLI DI GOMORY UNO ALLA VOLTA\t\t\t|")
    print("=================================================================")
    ## COSTRUISCO INSTANZA PROBLEMA CON CPLEX
    facilities, clients, f_vector, c_param, demands, capacity = parse_dat_file(data_path)
    assert np.all(np.array(demands) != 0)
    assert np.all(np.array(c_param) >= 0)
    m = len(facilities)
    n = len(clients)
    c,A,b = getProblemData(f_vector, c_param,demands, capacity)   # c -> m + n*m variabili
    nCols, nRows =(len(c)), (len(b))                    # b -> n + m*n + m vincoli (clienti + clienti*faciliy + facility)
    
    # inizializzo instanza delle variabili con i telaviti UB e LB
    names, lower_bounds, upper_bounds,constraint_senses,constraint_names = initializeInstanceVariables(n,m) 
    
    
    nCols= nCols + n*m + m    # numero variabili totali 'y' + 'x' + 's' (clienti*facility + facility)

    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_results_stream(None)
    prob.set_problem_name("istanza1") #nominiamo istanza
    prob.objective.set_sense(prob.objective.sense.minimize)
    params = prob.parameters
    params.preprocessing.presolve.set(0) 
    params.preprocessing.linear.set(0)
    params.preprocessing.reduce.set(0)
    prob.variables.add(names=names)
    
    
    # Add variables 
    for i in range(nCols-(n*m + m)):
        prob.variables.set_lower_bounds(i, lower_bounds[i])
        prob.variables.set_upper_bounds(i, upper_bounds[i])
    # Add slack
    for i in range(nCols-(n*m + m),nCols):
        prob.variables.set_lower_bounds(i, lower_bounds[i])

    #Add slack to constraints
    A = A.tolist()
                
    for i in range(nRows):
        prob.linear_constraints.add(
            lin_expr= [cplex.SparsePair(ind= [j for j in range(nCols)], 
            val= A[i])], 
            rhs= [b[i]], 
            names = [constraint_names[i]], 
            senses = [constraint_senses[i]]
        )
        
    # aggiungo funzione obiettivo -----------------------------------------------------------
    # con le variabili reali (quindi non considero quelle di slack)
    for i in range(nCols-(n*m + m)): 
        prob.objective.set_linear([(i, c[i])])

    t0 = time.time()
    prob.solve()
    #print(f"valore soluzione cplex = {prob.solution.get_objective_value()}")
    #print("-----------------calcolo tableau ---------------------")
    n_cuts, b_bar = get_tableau(prob,A,b)
    print("\tPossibili tagli:", n_cuts)
    #b_bar vettore dei termini noti
    
    #Print(f"Mi preparo per {n_cuts} tagli")

    varnames = prob.variables.get_names()
    gc_lhs, gc_rhs = initialize_fract_gc(n_cuts,nCols , prob, varnames, b_bar)
    #gc_lhs: matrice dei coefficienti delle disuguaglianze Gomory (left-hand side)
    #gc_rhs: termini noti (right-hand side) delle disuguaglianze

    #print("Genero tagli")
    cuts, cuts_limits, cut_senses = generate_gc(prob, A, gc_lhs, gc_rhs, names) 
    # cuts: lista dei coefficienti dei tagli (uno per ciascun vincolo)
    # cuts_limits: lista dei termini noti (right-hand side) dei tagli
    # cut_senses: lista dei sensi dei vincoli (es. <= → 'L')

    #print("Print nel log")
    #print_solution(prob)

    #print("Print a schermo")
    #print(f"Numero di tagli generati: {len(cuts)}")
    counter = 0 
    obj_prev = 0
    print("=================================================================")
    for i in range(len(cuts)):
        
        # 1. Estrai il taglio corrente
        indici = [j for j, val in enumerate(cuts[i]) if val != 0]
        valori = [cuts[i][j] for j in indici]
        # 2. Aggiungilo al modello CPLEX
        prob.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=indici, val=valori)],
            senses=[cut_senses[i]],
            rhs=[cuts_limits[i]]
        )

        # 3. Risolvi il problema aggiornato
        prob.solve()

        # 4. Stampa la nuova soluzione
       # print(prob)  # Assicurati di aver definito questa funzione
        #print("\n========== SOLUZIONE CORRENTE ==========")
        try:
            obj_value_cplex = prob.solution.get_objective_value()
            var_names = prob.variables.get_names()
            var_values = prob.solution.get_values()
            if obj_value_cplex == obj_prev : 
                counter +=1
            else : 
                counter = 0 
            obj_prev = obj_value_cplex
            print("-----------------------------------------------------------------")
            if i<10:
                print(f"| {i}° iterazione\t\tValore funzione obiettivo: {obj_value_cplex:.4f}\t|")
            else:
                print(f"| {i}° iterazione\tValore funzione obiettivo: {obj_value_cplex:.4f}\t|")
                
            # print("Valori delle variabili:")
            # for name, val in zip(var_names, var_values):
            #     if abs(val) > 1e-6:  # evita di stampare zeri numerici
            #         print(f"  {name} = {val:.4f}")
        except Exception as e:
            print("Errore nel recupero della soluzione:", e)
        if counter ==10: 
            print("-----------------------------------------------------------------")
            print(f"\tDopo {i} iterazioni non ci sono miglioramenti")
            break
    time_step = time.time() - t0
    obj_step = obj_value_cplex
    gap_step = compute_gap(obj_step, obj_int)
    
    # 4. Gomory: tutti i tagli
    print("=================================================================")
    print("|\t\tTAGLI DI GOMORY TUTTI INSIEME\t\t\t|")
    print("=================================================================")
    ampl_all = AMPL()
    ampl_relax.set_option('show_stats', 0)
    ampl_relax.set_option('solver_msg', 0)
    ampl_all.read(mod_path_relax)
    ampl_all.read_data(data_path)
    obj_all, time_all, iter_all = solve_with_gomory(ampl_all, all_cuts=True)
    gap_all = compute_gap(obj_all, obj_int)
    print(f"\n\tValore soluzione Rilassata: \t\t {obj_relax:.4f}")
    print(f"\tValore soluzione rilassata dopo i tagli: {obj_all:.4f}")
    print(f"\tGap Gomory all: {gap_all:.4f}%")

    return {
        "istanza": os.path.basename(data_path),
        "obj_int": obj_int,
        "time_int": time_int,
        "obj_relax": obj_relax,
        "obj_cplex" : obj_value_cplex,
        "time_relax": time_relax,
        "gap_relax": gap_relax,
        "obj_gomory_all": obj_all,
        "time_gomory_all": time_all,
        # "iter_gomory_all": iter_all,
        "gap_gomory_all": gap_all,
        "obj_gomory_step": obj_step,
        "time_gomory_step": time_step,
        # "iter_gomory_step": iter_step,
        "gap_gomory_step": gap_step,
        "nota": ""
    }



#*************************#
#   ---     MAIN    ---   #
#*************************#

def main():
    modello_intero = "sscfl.mod"
    modello_relax = "sscfl_relax.mod"
    istanze = sorted([f for f in os.listdir() if f.startswith("cap") and f.endswith("41.dat")])
    risultati = []
    for ist in istanze:
        print("=================================================================")
        print(f"|\t\t\t\t\t\t\t\t|\n|\t\tElaborazione {ist}...\t\t\t|\n|\t\t\t\t\t\t\t\t|")
        try:
            res = run_sscfl_experiment(modello_intero, modello_relax, ist)
            risultati.append(res)
        except Exception as e:
            print(f"Errore su {ist}: {e}")
            risultati.append({"istanza": ist, "nota": f"Errore: {e}"})
            return

    # Risultati finali
    print("=================================================================")
    print("Risultati finali:")
    max_len = 15
    for res in risultati:
        for key, value in res.items():
            if type(value) is int or type(value) is float :
                print(f"{key.ljust(max_len)} :{value:.4f}")
            else:
                print(f"{key.ljust(max_len)} :{value}")

        print("-" * 40)

if __name__ == "__main__":
    main()