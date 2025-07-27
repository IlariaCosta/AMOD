import time
import csv
import os
import cplex
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
from gomory import (
    solve_with_gomory,
    tagli
)
import cut
 


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
    
    tagli(ampl_relax,mod_path_int, mod_path_relax, data_path)
   


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