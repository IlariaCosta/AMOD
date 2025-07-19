import time
import csv
import os
from amplpy import AMPL, Environment
 
def is_integral(var_values):
    return all(abs(val - round(val)) < 1e-8 for val in var_values)
 
def compute_gap(relaxed_val, optimal_val):
    return abs((relaxed_val - optimal_val) / optimal_val) * 100
 
def solve_with_gomory(ampl, all_cuts, max_iter=100):
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

    if all_cuts:
        print("GOMORY CUTS TUTTI INSIEME\n")
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
        print("GOMORY CUTS UNO ALLA VOLTA\n")
        iter_count = 0
        t0 = time.time()
        while True:
         ampl.solve()
         iter_count += 1

         y_vals = ampl.get_variable('y').get_values().to_dict()
        #  print("Valori di y:")
        #  for idx, val in y_vals.items():
        #     print(f"y[{idx}] = {val}")
        #  print("\n")

         x_vals = ampl.get_variable('x').get_values().to_dict()
        #  print("Valori di x:")
        #  for idx, val in x_vals.items():
        #     print(f"x[{idx}] = {val}")
        #  print("\n")

         print(f"\nIterazione {iter_count}: valore soluzione = {ampl.obj['TotalCost'].value()}")

         cut_added = False
            #  for i, val in y_vals.items():
            #     if val > 0 and val < 1 :
            #         floor_val = int(val)
            #         cut_name = f"gomory_cut_{iter_count}_{i}"
            #         ampl.eval(f"subject to {cut_name}: y[{i}] <= {floor_val};")
            #         print(f"➕ Aggiunto taglio {cut_name}: y[{i}] <= {floor_val}")
            #         cut_added = True
            #         break

            #  for i, val in x_vals.items():
            #     if val > 0 and val < 1 :
            #         floor_val = int(val)
            #         i_str = "_".join(str(k) for k in i)
            #         cut_name = f"gomory_cut_{iter_count}_{i_str}"
            #         ampl.eval(f"subject to {cut_name}: x[{i[0]}, {i[1]}] <= {floor_val};")
            #         print(f"➕ Aggiunto taglio {cut_name}: x[{i[0]}, {i[1]}] <= {floor_val}")
            #         cut_added = True
            #         break

         EPSILON = 1e-5  # tolleranza per verificare se una variabile è frazionaria

         # ➤ Cerca frazionarie in y
         for i, val in y_vals.items():
            if 0 < val < 1 :
                floor_val = int(val)
                cut_name = f"gomory_cut_{iter_count}_y_{i}"
                ampl.eval(f"subject to {cut_name}: y[{i}] <= {floor_val};")
                print(f"➕ Aggiunto taglio {cut_name}: y[{i}] <= {floor_val}")
                cut_added = True
                break  # esci dopo il primo taglio

         # ➤ Cerca frazionarie in x
         for i, val in x_vals.items():
            if 0 < val < 1 :
                floor_val = int(val)
                i_str = "_".join(str(k) for k in i)
                cut_name = f"gomory_cut_{iter_count}_x_{i_str}"
                ampl.eval(f"subject to {cut_name}: x[{i[0]}, {i[1]}] <= {floor_val};")
                print(f"➕ Aggiunto taglio {cut_name}: x[{i[0]}, {i[1]}] <= {floor_val}")
                cut_added = True
                break  # esci dopo il primo taglio


         if not cut_added:
            print("✅ Nessuna variabile frazionaria trovata, soluzione intera raggiunta.")
            break

         if iter_count >= max_iter:
            print("⚠️ Raggiunto numero massimo iterazioni. Termino.")
            break

        
        elapsed = time.time() - t0
        obj = ampl.obj['TotalCost'].value()
        return obj, elapsed, iter_count
    
 
def run_ufl_experiment(mod_path_int, mod_path_relax, data_path):
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
    print ("%d\n", obj_int);
    

   # 2. Rilassamento ottenuto da modello intero
    print("Soluzione rilassata ->");
    ampl_relax = AMPL()
    ampl_relax.set_option('solver', 'cplex')
    ampl_relax.read(mod_path_relax)
    ampl_relax.read_data(data_path)

    # 🔁 Rilassa da Python tutte le variabili intere
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
    print("%d\n", obj_relax)
    
    gap_relax = compute_gap(obj_relax, obj_int)
    
    x_vals = list(ampl_relax.get_variable('y').get_values().to_dict().values())

    # print("Valori di y:")
    # for idx, val in y_vals.items():
    #     print(f"y[{idx}] = {val}")
    
    already_integer = is_integral(x_vals)


    if already_integer:
        print ("SOLUZIONE INTERA\n");
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

    print("***** GOMORY CUTS *****");
    # 3. Gomory: tutti i tagli
    ampl_all = AMPL()
    ampl_all.read(mod_path_relax)
    ampl_all.read_data(data_path)
    obj_all, time_all, iter_all = solve_with_gomory(ampl_all, all_cuts=True)
    gap_all = compute_gap(obj_all, obj_int)
    print(f"Obj rilassato: {obj_relax}")
    print(f"Obj dopo tagli Gomory: {obj_all}")
    print(f"Gap Gomory all: {gap_all:.4f}%")


    # 4. Gomory: uno alla volta
    ampl_step = AMPL()
    ampl_step.read(mod_path_relax)
    ampl_step.read_data(data_path)
    obj_step, time_step, iter_step = solve_with_gomory(ampl_step, all_cuts=False)
    gap_step = compute_gap(obj_step, obj_int)

    return {
        "istanza": os.path.basename(data_path),
        "obj_int": obj_int,
        "time_int": time_int,
        "obj_relax": obj_relax,
        "time_relax": time_relax,
        "gap_relax": gap_relax,
        "obj_gomory_all": obj_all,
        "time_gomory_all": time_all,
        "iter_gomory_all": iter_all,
        "gap_gomory_all": gap_all,
        #  "obj_gomory_step": obj_step,
        # "time_gomory_step": time_step,
        # "iter_gomory_step": iter_step,
        #  "gap_gomory_step": gap_step,
        "nota": ""
    }

 
def main():
    modello_intero = "sscfl.mod"
    modello_relax = "sscfl_relax.mod"
    istanze = sorted([f for f in os.listdir() if f.startswith("cap") and f.endswith(".dat")])
    print("File .dat trovati:", istanze)
    risultati = []
    istanze = istanze[1:2]
    for ist in istanze:
        print(f"\n🔄 Elaborazione {ist}...")
        try:
            res = run_ufl_experiment(modello_intero, modello_relax, ist)
            risultati.append(res)
        except Exception as e:
            print(f"❌ Errore su {ist}: {e}")
            risultati.append({"istanza": ist, "nota": f"Errore: {e}"})
            return

    # Scrittura CSV
    if risultati:
        with open("risultati_ufl.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=risultati[0].keys())
            writer.writeheader()
            writer.writerows(risultati)
    print("\n✅ Tutto completato. Risultati salvati in risultati_ufl.csv.")
    print("Risultati finali:")
    for res in risultati:
        print(res)
        print("\n")



if __name__ == "__main__":
    main()