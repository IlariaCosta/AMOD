import time
import csv
import os
from amplpy import AMPL, Environment
 
def is_integral(var_values):
    return all(abs(val - round(val)) < 1e-8 for val in var_values)
 
def compute_gap(relaxed_val, optimal_val):
    return abs((relaxed_val - optimal_val) / optimal_val) * 100

def parse_dat_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    facilities = []
    clients = []
    f_param = {}
    c_param = []

    reading = None

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("set FACILITIES"):
            facilities = list(map(int, line.split(":=")[1].replace(";", "").split()))
            continue
        elif line.startswith("set CLIENTS"):
            clients = list(map(int, line.split(":=")[1].replace(";", "").split()))
            continue
        elif line.startswith("param f"):
            reading = "f"
            continue
        elif line.startswith("param c"):
            reading = "c"
            c_param = []
            continue
        elif line.startswith("param d"):
            reading = None  # ignoriamo
            continue

        # Lettura dei costi f
        if reading == "f":
            if ";" in line:
                line = line.replace(";", "")
                reading = None
            parts = line.split()
            if len(parts) == 2:
                f_param[int(parts[0])] = float(parts[1])  

        # Lettura della matrice dei costi c
        elif reading == "c":
            if ":=" in line:
                continue
            if ";" in line:
                line = line.replace(";", "")
                if line.strip():
                    c_param.append(list(map(float, line.split())))  
                reading = None
            else:
                if line.strip():
                    c_param.append(list(map(float, line.split())))  

    # Ordina f_param secondo l'ordine dei facilities
    f_vector = [f_param[i] for i in sorted(facilities)]
    return facilities, clients, f_vector, c_param

def build_A_b(facilities, clients):
    m = len(facilities)
    n = len(clients)

    var_count = m * n + m
    row_count = n + m * n

    A = [[0 for _ in range(var_count)] for _ in range(row_count)]
    b = [0 for _ in range(row_count)]

    def x_index(i, j):
        return i * n + j

    def y_index(i):
        return m * n + i

    # Cliente j deve essere assegnato a una sola facility
    for j in range(n):
        for i in range(m):
            A[j][x_index(i, j)] = 1
        b[j] = 1

    # x_ij ≤ y_i
    row = n
    for i in range(m):
        for j in range(n):
            A[row][x_index(i, j)] = 1
            A[row][y_index(i)] = -1
            b[row] = 0
            row += 1

    print("A e b pronte\n")

    # # Stampa matrice A
    # print("Matrice A:")
    # for row in A:
    #     print(row)

    # # Stampa vettore b
    # print("\nVettore b:")
    # print(b)

    return A, b




#*************************#
#   ---    GOMORY   ---   # -------------------------------------------------------------------------------------------------
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
    

#*************************#
#   ---     SSCFL   ---   #
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
    facilities, clients, f_vector, _ = parse_dat_file(data_path)
    A, b_vec = build_A_b(facilities, clients)

    # Ricava il nome base del file (senza estensione .dat)
    basename = os.path.splitext(os.path.basename(data_path))[0]

    # Salvataggio di A e b
    with open(f"A_{basename}.txt", "w") as fa:
        for row in A:
            fa.write(" ".join(map(str, row)) + "\n")

    with open(f"b_{basename}.txt", "w") as fb:
        fb.write("\n".join(map(str, b_vec)))



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