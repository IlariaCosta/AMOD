import time


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
    
