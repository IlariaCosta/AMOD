import time
import csv
import os
from amplpy import AMPL, Environment
 
def is_integral(var_values):
    return all(abs(val - round(val)) < 1e-5 for val in var_values)
 
def compute_gap(relaxed_val, optimal_val):
    return abs((relaxed_val - optimal_val) / optimal_val) * 100
 
def solve_with_gomory(ampl, all_cuts=True, max_iter=100):
    for var in ampl.get_variables().values():
        for var in ampl.get_variables().values():
            if var.num_instances() > 1:
                for i in var.get_values().index():
                    var[i].set_integer(False)
            else:
                var.set_integer(False)

    ampl.set_option('presolve', 0)
    ampl.set_option('cut_generation', 'gomory')
    ampl.set_option('solver', 'cplex')
 
    if all_cuts:
        ampl.set_option('gomory_cuts', -1)  # all available
        t0 = time.time()
        ampl.solve()
        elapsed = time.time() - t0
        obj = ampl.obj['TotalCost'].value()
        return obj, elapsed, 1
    else:
        # un taglio per volta
        iter_count = 0
        t0 = time.time()
        while True:
            ampl.set_option('gomory_cuts', 1)
            ampl.solve()
            iter_count += 1
            var_values = list(ampl.get_variable('x').get_values().to_dict().values())
            if is_integral(var_values) or iter_count >= max_iter:
                break
        elapsed = time.time() - t0
        obj = ampl.obj['TotalCost'].value()
        return obj, elapsed, iter_count
 
def run_ufl_experiment(mod_path_int, mod_path_relax, data_path):
    # 1. Soluzione ottima intera
    ampl = AMPL()
    ampl.set_option('solver', 'cplex')
    ampl.read(mod_path_int)
    ampl.read_data(data_path)

    ampl.set_option('cplex_options', 'mipgap=0')
    t0 = time.time()
    ampl.solve()
    time_int = time.time() - t0
    obj_int = ampl.obj['TotalCost'].value()

    # 2. Rilassamento lineare con modello già rilassato
    ampl_relax = AMPL()
    ampl_relax.set_option('solver', 'cplex')
    ampl_relax.read(mod_path_relax)
    ampl_relax.read_data(data_path)

    t0 = time.time()
    ampl_relax.solve()
    time_relax = time.time() - t0
    obj_relax = ampl_relax.obj['TotalCost'].value()
    gap_relax = compute_gap(obj_relax, obj_int)

    x_vals = list(ampl_relax.get_variable('x').get_values().to_dict().values())
    already_integer = is_integral(x_vals)

    if already_integer:
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

    # 3. Gomory: tutti i tagli
    ampl_all = AMPL()
    ampl_all.read(mod_path_relax)
    ampl_all.read_data(data_path)
    obj_all, time_all, iter_all = solve_with_gomory(ampl_all, all_cuts=True)
    gap_all = compute_gap(obj_all, obj_int)

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
        "obj_gomory_step": obj_step,
        "time_gomory_step": time_step,
        "iter_gomory_step": iter_step,
        "gap_gomory_step": gap_step,
        "nota": ""
    }

 
def main():
    modello_intero = "ufl.mod"
    modello_relax = "ufl_relax.mod"
    istanze = sorted([f for f in os.listdir() if f.startswith("cap") and f.endswith(".dat")])
    print("File .dat trovati:", istanze)
    risultati = []
    for ist in istanze:
        print(f"\n🔄 Elaborazione {ist}...")
        try:
            res = run_ufl_experiment(modello_intero, modello_relax, ist)
            risultati.append(res)
        except Exception as e:
            print(f"❌ Errore su {ist}: {e}")
            risultati.append({"istanza": ist, "nota": f"Errore: {e}"})

    # Scrittura CSV
    if risultati:
        with open("risultati_ufl.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=risultati[0].keys())
            writer.writeheader()
            writer.writerows(risultati)
    print("\n✅ Tutto completato. Risultati salvati in risultati_ufl.csv.")


if __name__ == "__main__":
    main()