import time
import csv
import os
import cplex
import numpy as np
from amplpy import AMPL, Environment
from utils import (
    compute_gap,
    is_integral,
    parse_dat_file,
    build_A_b,
    invert_matrix,
    mat_vec_mul,
    mat_mul
)


#*************************#
#   ---    GOMORY   ---   # 
#*************************#
 
def solve_with_gomory(ampl, all_cuts, max_iter=100, min_improvement=1e-3):
    variables = ampl.get_variables()
    for var in variables:   
        try:
          if var.is_integer():
             var.set_integer(False)
        except:
         continue


    ampl.set_option('presolve', 0)
    ampl.set_option('cut_generation', 'gomory')
    ampl.set_option('cplex_options', 'display=2 simplex display=2')
    ampl.set_option('solver', 'cplex')
    ampl.set_option('display', 1)

    already_cut = set()

    if all_cuts:
        print("\nGOMORY CUTS TUTTI INSIEME\n")
        ampl.set_option('gomory_cuts', -1)  # all available
        
        t0 = time.time()
        ampl.solve()
        elapsed = time.time() - t0
        var_values = list(ampl.get_variable('y').get_values().to_dict().values())
        print(var_values)
     
        obj = ampl.obj['TotalCost'].value()
        return obj, elapsed, 1
    else:
        # un taglio per volta
        print("\nGOMORY CUTS UNO ALLA VOLTA\n")
        iter_count = 0
        no_progress_count = 0  # Conta quante iterazioni con miglioramento trascurabile
        prev_obj = float('inf')
        t0 = time.time()

        while True:
         ampl.solve()
         iter_count += 1
         if iter_count > max_iter:
            print("Raggiunto numero massimo iterazioni. Termino.")
            break

         y_vals = ampl.get_variable('y').get_values().to_dict()
         x_vals = ampl.get_variable('x').get_values().to_dict()

         obj_curr = ampl.obj['TotalCost'].value()
         improvement = abs(prev_obj - obj_curr)

         print(f"\nIterazione {iter_count}: valore soluzione = {obj_curr}, miglioramento = {improvement:.6f}")
       
         prev_obj = obj_curr

         cut_added = False
        
         # ➤ Tagli per y (uno per volta)
         for i, val in y_vals.items():
                val_rounded = round(val)
                if abs(val - round(val)) > 1e-5:
                    bound = ">= 1" if val > 0.5 else "<= 0.5"
                    ampl.eval(f"subject to cut_y_{iter_count}_{i}: y[{i}] {bound};")
                    print(f"Taglio: y[{i}] {bound} (valore attuale = {val})")
                    already_cut.add(i)
                    cut_added = True
                    break  # SOLO UNO PER ITERAZIONE

         if not cut_added:
         # ➤ Tagli per x (uno per volta)
            for i, val in x_vals.items():
                 val_rounded = round(val)
                 if abs(val - round(val)) > 1e-5:
                            i_str = "_".join(str(k) for k in i)
                            bound = ">= 1" if val > 0.5 else "<= 0.5"
                            ampl.eval(f"subject to cut_x_{iter_count}_{i_str}: x[{i[0]}, {i[1]}] {bound};")
                            print(f"Taglio: x[{i[0]}, {i[1]}] {bound} (valore attuale = {val})")
                            already_cut.add(i)
                            cut_added = True
                            break  # SOLO UNO PER ITERAZIONE




         if not cut_added:
            print("Nessuna variabile frazionaria trovata, soluzione intera raggiunta.")
            # print("Valori di y:")
            # for idx, val in y_vals.items():
            #     print(f"y[{idx}] = {val}")
            # print("\n")

            # print("Valori di x:")
            # for idx, val in x_vals.items():
            #     print(f"x[{idx}] = {val}")
            # print("\n")

            print("Variabili frazionarie:")
            for i, val in y_vals.items():
                if abs(val - round(val)) > 1e-5:
                    print(f"y[{i}] = {val}")
            for i, val in x_vals.items():
                if abs(val - round(val)) > 1e-5:
                    print(f"x[{i}] = {val}")
            break
         
         if improvement < min_improvement:
                no_progress_count += 1
                print(f"Nessun miglioramento. ({no_progress_count}/5)") 
         else:
            no_progress_count = 0  # Reset se c'è miglioramento


         if no_progress_count >= 5:
            print("Terminato: 5 iterazioni senza miglioramento.")
            break

        

        
        elapsed = time.time() - t0
        obj = ampl.obj['TotalCost'].value()
        return obj, elapsed, iter_count
    

def tagli(ampl,mod_path_int, mod_path_relax, data_path):
     # Costruzione A, b
    print("ORA CALCOLO A E B \n");
    #facilities, clients, f_vector, _ = parse_dat_file(data_path)
    # Leggi direttamente dal file .dat
    facilities, clients, f_vector, c_param, demands = parse_dat_file(data_path)
    m = len(facilities)
    n = len(clients)


    A, b_vec = build_A_b(facilities, clients, demands)

    # Ricava il nome base del file (senza estensione .dat)
    basename = os.path.splitext(os.path.basename(data_path))[0]

    # Salvataggio di A e b
    with open(f"A_{basename}.txt", "w") as fa:
        for row in A:
            fa.write(" ".join(map(str, row)) + "\n")

    with open(f"b_{basename}.txt", "w") as fb:
        fb.write("\n".join(map(str, b_vec)))

    
    # --- Ora procediamo a costruire B, N e fare i calcoli ---
    # Ricava gli indici delle variabili basiche e non basiche
    basic_indices = []
    nonbasic_indices = []

    # Assumiamo che ampl.getVariables() ti dia variabili in ordine di colonne di A
    # Dovrai adattare se l'ordine è diverso

    variables = ampl.get_variables()
    basic_indices = []
    nonbasic_indices = []
    col_idx = 0

    variables = ampl.get_variables()

    # print("Tipo variables:", type(variables))
    # print("variables:", variables)

    # Proviamo a fare list() per vedere cosa contiene
    # variables_list = list(variables)
    # print("Lista variables:", variables_list)

    # for i, elem in enumerate(variables_list):
        # print(f"Elemento {i}: tipo {type(elem)} - valore {elem}")

    # Se elem è una tupla (nome, Variable), possiamo usarlo:
    # for var_name, var in variables_list:
        # print(f"Var name: {var_name}, Var type: {type(var)}")


    #from amplpy import Solution

    #solution = Solution(ampl)

    basic_indices = []
    nonbasic_indices = []
    col_idx = 0

    variables_list = list(ampl.get_variables())

    status = ampl.get_parameter('solve_result').value()
    # print("Solve result:", status)
    
    for var_name, var in variables_list:
        # print("\n")
        #print(f"var_name: {var_name}, tipo var: {type(var)}")
        #print(list(var.get_values().to_dict()))  ## get_values da un dataframe
        
        try:
            keys = list(var.get_values().to_dict())
            #print(f"  keys: {keys}")
        except Exception as e:
            # print(f"  Nessun keys() per {var_name}: {e}")
            keys = None

        if keys:
            val = var.get_values().to_dict()
            # print("qui")
            
            
            
            for idx in keys:
                status = var[idx].sstatus()
                # print(idx)
                # status = var[idx].sstatus()
                # print(var[idx].sstatus())
                
                #if status is 'none':
                # if val[idx] !=0:
                    # status = 'Basic' 
                    # print(val[idx])
                    # print('\n')
                # else :
                    # status = 'Nonbasic'
            
                #status = solution.get_var_status(var[idx])
                # print(f"  {var_name}[{idx}] status: {status}")
                
                #print(var[idx].get_values())
                if status == 'bas':
                    basic_indices.append(col_idx)
                else:
                    nonbasic_indices.append(col_idx)
                col_idx += 1
        else:

            try:
                status = var.sstatus()
                val = list(var.get_values())
                #status = solution.get_var_status(var)
                print(f"  qui si {var_name} status (scalare): {status}")
            except Exception as e:
                # print(f"  Errore get_var_status per {var_name}: {e}")
                status = 'Nonbasic'
                
            # if val !=0:
                # status = 'Basic' 
                
            # else :
                # status = 'Nonbasic'
                
            if status == 'bas':
                basic_indices.append(col_idx)
            else:
                nonbasic_indices.append(col_idx)
            col_idx += 1

    # print("Basic indices:", basic_indices)
    # print("Nonbasic indices:", nonbasic_indices)
    # print(len(basic_indices))



    


    B = [[row[i] for i in basic_indices] for row in A]
    # print(len(A), len(A[0]))
    
    N = [[row[i] for i in nonbasic_indices] for row in A]
    # print(len(N), len(N[0]))


    # Inverti B
    try:

        # print(np.linalg.matrix_rank(B))
        B_inv = np.linalg.inv(B)
        # B_inv = invert_matrix(B)

    except ValueError:
        print("⚠️ Matrice B non invertibile")
        B_inv = None  # O gestisci l'errore come preferisci

    if B_inv is not None:
        print(len(N), len(N[0]))
        rhs = np.dot(B_inv , b_vec)
        # rhs = mat_vec_mul(B_inv, b_vec)

        # Calcola -B_inv * N
        # BN = mat_mul(B_inv, N)
        BN = np.dot(B_inv,N)
        reduced_matrix = [[-x for x in row] for row in BN]

        # Stampa per debug (opzionale)
        # print("rhs (B_inv * b):", rhs)
        # print("reduced_matrix (-B_inv * N):")
        # for r in reduced_matrix:
            # print(r)

        # A questo punto puoi procedere a cercare la riga con la massima parte frazionaria e
        # costruire il taglio di Gomory come da tua logica
        
        

    return