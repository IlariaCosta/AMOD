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
    solve_with_gomory,
    tagli
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
    ampl.read(mod_path_int)
    ampl.read_data(data_path)

    #ampl.set_option('cplex_options', 'mipgap=0')
    t0 = time.time()
    ampl.solve()
    
    time_int = time.time() - t0
    
    obj_int = ampl.obj['TotalCost'].value()
    print(f"\n\tSoluzione intera  = {obj_int}");

    # 2. CARICO MODELLO RILASSATO
    print("=================================================================")
    ampl_relax = AMPL()
    ampl_relax.set_option('solver', 'cplexamp')
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
    # print("Valori di y:")
    # for idx, val in y_vals.items():
    #     print(f"y[{idx}] = {val}")

    time_relax = time.time() - t0
    
    obj_relax = ampl_relax.obj['TotalCost'].value()
    print(f"\n\tSoluzione rilassata = {obj_relax:.4f}");
    gap_relax = compute_gap(obj_relax, obj_int)
    
    x_vals = list(ampl_relax.get_variable('y').get_values().to_dict().values())

    # print("Valori di y:")
    # for idx, val in y_vals.items():
    #     print(f"y[{idx}] = {val}")
    
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
    
    print("=================================================================")
    print(f"|\t\t\t\t\t\t\t\t|\n|\t\t\tGOMORY CUTS\t\t\t\t|\n|\t\t\t\t\t\t\t\t|")
    print("=================================================================")
    print("|\t\tTAGLI DI GOMORY UNO ALLA VOLTA\t\t\t|")
    print("=================================================================")
    ## COSTRUISCO INSTANZA PROBLEMA CON CPLEX
    facilities, clients, f_vector, c_param, demands = parse_dat_file(data_path)
    assert np.all(np.array(demands) != 0)
    assert np.all(np.array(c_param) >= 0)
    m = len(facilities)
    n = len(clients)
    c,A,b = getProblemData(f_vector, c_param,demands)   # c -> m + n*m variabili
    nCols, nRows =(len(c)), (len(b))                    # b -> n + m*n vincoli

    # inizializzo instanza delle variabili con i telaviti UB e LB
    (names, 
     lower_bounds, 
     upper_bounds,
     constraint_senses,
     constraint_names) = initializeInstanceVariables(n,m) 
    
    nCols= nCols + n*m          # numero variabili totali 'y' + 'x' + 's'
    
    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_results_stream(None)
    prob.set_problem_name("istanza1") #nominiamo istanza
    prob.objective.set_sense(prob.objective.sense.minimize)
    params = prob.parameters
    params.preprocessing.presolve.set(0) 
    params.preprocessing.linear.set(0)
    params.preprocessing.reduce.set(0)
    # print("NUMERO VARIABILI DEL PROBLEMA ",prob.variables.get_num())
    prob.variables.add(names=names)
        
    # agginiungo variabili 'y' e 'x'
    for i in range(nCols-2*n*m):
        prob.variables.set_lower_bounds(i, lower_bounds[i])
        prob.variables.set_upper_bounds(i, upper_bounds[i])
    # aggiungo variabili di slack 's'
    for i in range(nCols-2*m*n,nCols):
        prob.variables.set_lower_bounds(i, lower_bounds[i])


    A = A.tolist()
    # aggiungo i vincoli         
    for i in range(nRows):
        # print(f"Row {i}: len(A[i]) = {len(A[i])}, expected {nCols}")
        # print(f"rhs: {b[i]}, sense: {constraint_senses[i]}, name: {constraint_names[i]}")
        prob.linear_constraints.add(
            lin_expr= [cplex.SparsePair(ind= [j for j in range(nCols)], 
            val= A[i])], 
            rhs= [b[i]], 
            names = [constraint_names[i]], 
            senses = [constraint_senses[i]]
        )
        
    # aggiungo la funzione obiettivo -----------------------------------------------------------
    for i in range(nCols-n*m): 
        prob.objective.set_linear([(i, c[i])])
 
    #print("risolvo istanza cplex")
    prob.solve()
    #print(f"valore soluzione cplex = {prob.solution.get_objective_value()}")
    #print("-----------------calcolo tableau ---------------------")
    n_cuts, b_bar = get_tableau(prob,A,b)
    #print(b_bar)
    print("Possibili tagli", n_cuts)
    #b_bar vettore dei termini noti
    #rint("Mi preparo per tagli")
    varnames = prob.variables.get_names()
    gc_lhs, gc_rhs = initialize_fract_gc(n_cuts,nCols , prob, varnames, b_bar)
    
    #print(f"{gc_lhs}>={gc_rhs}")
    #for i in range(len(gc_rhs)):
    #    print(f"{gc_lhs[i]}>= {gc_rhs[i]}")
    #gc_lhs: matrice dei coefficienti delle disuguaglianze Gomory (left-hand side)
    #gc_rhs: termini noti (right-hand side) delle disuguaglianze

    #print("Genero tagli")
    cuts, cuts_limits, cut_senses = generate_gc(prob, A, gc_lhs, gc_rhs, names) 
    # cuts: lista dei coefficienti dei tagli (uno per ciascun vincolo)
    # cuts_limits: lista dei termini noti (right-hand side) dei tagli
    # cut_senses: lista dei sensi dei vincoli (es. <= → 'L')

    #print(cut_senses)
    print_solution(prob)

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
        try:
            obj_value_cplex = prob.solution.get_objective_value()
            var_names = prob.variables.get_names()
            var_values = prob.solution.get_values()
            if abs(obj_value_cplex - obj_prev)< 1e-3 :
                counter +=1
            else : 
                counter = 0 
            obj_prev = obj_value_cplex
            print("-----------------------------------------------------------------")
            if i<10:
                print(f"| {i}° iterazione   Valore funzione obiettivo: {obj_value_cplex:.4f}\t|")
            else:
                print(f"| {i}° iterazione  Valore funzione obiettivo: {obj_value_cplex:.4f}\t|")
        except Exception as e:
            print("Errore nel recupero della soluzione:", e)
        if counter ==5: 
            print("-----------------------------------------------------------------")
            print(f"La soluzione non migliora ulteriormente dopo iterazione {i}")
            break
    time_step = time.time() - t0
    obj_step = obj_value_cplex
    gap_step = compute_gap(obj_step, obj_int)
    

    
    #print("\n***** GOMORY CUTS *****");
    # 3. Gomory: tutti i tagli
    print("=================================================================")
    print("|\t\tTAGLI DI GOMORY TUTTI INSIEME\t\t\t|")
    print("=================================================================")
    ampl_all = AMPL()
    ampl_relax.set_option('solver_msg', 0)
    ampl_all.read(mod_path_relax)
    ampl_all.read_data(data_path)
    obj_all, time_all, iter_all = solve_with_gomory(ampl_all, all_cuts=True)
    gap_all = compute_gap(obj_all, obj_int)
    print(f"\n\tValore soluzione Rilassata: \t\t {obj_relax:.4f}")
    print(f"\tValore soluzione rilassata dopo i tagli: {obj_all:.4f}")
    print(f"\tGap Gomory all: {gap_all:.4f}%")


    # 4. Gomory: uno alla volta
    # ampl_step = AMPL()
    # ampl_step.read(mod_path_relax)
    # ampl_step.read_data(data_path)
    # obj_step, time_step, iter_step = solve_with_gomory(ampl_step, all_cuts=False)
    # gap_step = compute_gap(obj_step, obj_int)
    #-----------------------------------------------------------------------------------------------------------------------
    
    # tagli(ampl_relax,mod_path_int, mod_path_relax, data_path)
   


    return {
        "istanza": os.path.basename(data_path),
        "facility" : m,
        "clienti" : n,
        "obj_int": obj_int,
        "time_int": time_int,
        "obj_relax": obj_relax,
        "time_relax": time_relax,
        "gap_relax": gap_relax,
        "time_relax": time_relax,
        "gap_relax": gap_relax,
        "obj_gomory_all": obj_all,
        "time_gomory_all": time_all,
        # "iter_gomory_all": iter_all,
        "gap_gomory_all": gap_all,
        "obj_gomory_step": obj_step,
        "time_gomory_step": time_step,
        "numero tagli possibili" : n_cuts,
        "iter_gomory_step": i,
        "gap_gomory_step": gap_step,
        "nota": ""
    }



#*************************#
#   ---     MAIN    ---   #
#*************************#

def main():
    modello_intero = "ufl.mod"
    modello_relax = "ufl_relax.mod"
    istanze = sorted([f for f in os.listdir() if f.startswith("cap") and f.endswith(".dat")])
    print("File .dat trovati:", istanze)
    risultati = []
    #istanze = istanze[1:2]
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

    # Scrittura CSV
    if risultati:
        with open("risultati_ufl.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=risultati[0].keys())
            writer.writeheader()
            writer.writerows(risultati)
    print("\nTutto completato. Risultati salvati in risultati_ufl.csv.")

    # Risultati finali
    print("Risultati finali:")
    for res in risultati:
        max_len = max(len(k) for k in res.keys())  # per allineare i due punti
        for key, value in res.items():
            if type(value) is float :
                if key.startswith("gap"):
                    print(f"{key.ljust(max_len)} :{value:.8f}")
                else:
                    print(f"{key.ljust(max_len)} :{value:.4f}")
            else:
                print(f"{key.ljust(max_len)} :{value}")
        print("-" * 40)




if __name__ == "__main__":
    main()