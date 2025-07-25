import time
import csv
import os
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

import gomory
 


#*************************#
#   ---     UFL     ---   #
#*************************#
 
def run_sscfl_experiment(mod_path_int, mod_path_relax, data_path):
    # 1. Soluzione ottima intera
    print("Soluzione intera ->");
    ampl = AMPL()
    ampl.set_option('solver', 'cplex')
    ampl.read(mod_path_int)
    ampl.read_data(data_path)

    ampl.set_option('cplex_options', 'mipgap=0')
    t0 = time.time()
    ampl.solve()
    
    time_int = time.time() - t0
    
    obj_int = ampl.obj['TotalCost'].value()


    # 2. Rilassamento ottenuto da modello intero
    print("Soluzione rilassata ->");
    ampl_relax = AMPL()
    ampl_relax.set_option('solver', 'cplex')
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

    print("\n***** GOMORY CUTS *****");
    # 3. Gomory: tutti i tagli
    # ampl_all = AMPL()
    # ampl_all.read(mod_path_relax)
    # ampl_all.read_data(data_path)
    # obj_all, time_all, iter_all = solve_with_gomory(ampl_all, all_cuts=True)
    # gap_all = compute_gap(obj_all, obj_int)
    # print(f"Obj rilassato: {obj_relax}")
    # print(f"Obj dopo tagli Gomory: {obj_all}")
    # print(f"Gap Gomory all: {gap_all:.4f}%")


    # 4. Gomory: uno alla volta
    # ampl_step = AMPL()
    # ampl_step.read(mod_path_relax)
    # ampl_step.read_data(data_path)
    # obj_step, time_step, iter_step = solve_with_gomory(ampl_step, all_cuts=False)
    # gap_step = compute_gap(obj_step, obj_int)
    #-----------------------------------------------------------------------------------------------------------------------

    # Costruzione A, b
    print("ORA CALCOLO A E B \n");
    #facilities, clients, f_vector, _ = parse_dat_file(data_path)
    # Leggi direttamente dal file .dat
    facilities, clients, f_vector, _ = parse_dat_file(data_path)
    m = len(facilities)
    n = len(clients)


    A, b_vec = build_A_b(facilities, clients)

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

    print("Tipo variables:", type(variables))
    print("variables:", variables)

    # Proviamo a fare list() per vedere cosa contiene
    variables_list = list(variables)
    print("Lista variables:", variables_list)

    for i, elem in enumerate(variables_list):
        print(f"Elemento {i}: tipo {type(elem)} - valore {elem}")

    # Se elem è una tupla (nome, Variable), possiamo usarlo:
    for var_name, var in variables_list:
        print(f"Var name: {var_name}, Var type: {type(var)}")


    from amplpy import Solution

    solution = Solution(ampl)

    basic_indices = []
    nonbasic_indices = []
    col_idx = 0

    variables_list = list(ampl.get_variables())

    for var_name, var in variables_list:
        print(f"var_name: {var_name}, tipo var: {type(var)}")
        try:
            keys = list(var.keys())
            print(f"  keys: {keys}")
        except Exception as e:
            print(f"  Nessun keys() per {var_name}: {e}")
            keys = None

        if keys:
            for idx in keys:
                status = solution.get_var_status(var[idx])
                print(f"  {var_name}[{idx}] status: {status}")
                if status == 'Basic':
                    basic_indices.append(col_idx)
                else:
                    nonbasic_indices.append(col_idx)
                col_idx += 1
        else:
            try:
                status = solution.get_var_status(var)
                print(f"  {var_name} status (scalare): {status}")
            except Exception as e:
                print(f"  Errore get_var_status per {var_name}: {e}")
                status = 'Nonbasic'
            if status == 'Basic':
                basic_indices.append(col_idx)
            else:
                nonbasic_indices.append(col_idx)
            col_idx += 1

    print("Basic indices:", basic_indices)
    print("Nonbasic indices:", nonbasic_indices)







    B = [[row[i] for i in basic_indices] for row in A]
    N = [[row[i] for i in nonbasic_indices] for row in A]


    # Inverti B
    try:
        B_inv = invert_matrix(B)
    except ValueError:
        print("⚠️ Matrice B non invertibile")
        B_inv = None  # O gestisci l'errore come preferisci

    if B_inv is not None:
        # Calcola rhs = B_inv * b_vec
        rhs = mat_vec_mul(B_inv, b_vec)

        # Calcola -B_inv * N
        BN = mat_mul(B_inv, N)
        reduced_matrix = [[-x for x in row] for row in BN]

        # Stampa per debug (opzionale)
        print("rhs (B_inv * b):", rhs)
        print("reduced_matrix (-B_inv * N):")
        for r in reduced_matrix:
            print(r)

        # A questo punto puoi procedere a cercare la riga con la massima parte frazionaria e
        # costruire il taglio di Gomory come da tua logica






    return {
        "istanza": os.path.basename(data_path),
        "obj_int": obj_int,
        "time_int": time_int,
        "obj_relax": obj_relax,
        "time_relax": time_relax,
        "gap_relax": gap_relax,
        # "obj_gomory_all": obj_all,
        # "time_gomory_all": time_all,
        # "iter_gomory_all": iter_all,
        # "gap_gomory_all": gap_all,
        #  "obj_gomory_step": obj_step,
        # "time_gomory_step": time_step,
        # "iter_gomory_step": iter_step,
        #  "gap_gomory_step": gap_step,
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
        print(f"\nElaborazione {ist}...")
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
            print(f"{key.ljust(max_len)} : {value}")
        print("-" * 40)




if __name__ == "__main__":
    main()