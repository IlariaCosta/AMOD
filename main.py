import os
import subprocess

MAX_ITER = 5
TOL = 1e-5

for iter in range(1, MAX_ITER + 1):
    print(f"=== Iterazione {iter} ===")

    # 1. Lancia AMPL per risolvere il rilassamento LP
    with open("solve.run", "w") as f:
        f.write(f"""
reset;
model ufl_relax.mod;
data cap101.dat;
param obj_lp;
option cplex_options 'cuts=0';
option solver cplex;
solve;
printf {{i in FACILITIES}} "%s %f\\n", i, y[i] > "y_sol.txt";
let obj_lp := TotalCost;
display obj_lp;
display y;
        """)

    subprocess.run(["ampl", "solve.run"])

    # 2. Lancia script Python per generare il taglio
    subprocess.run(["python", "build_cut.py", "y_sol.txt", str(iter)], stdout=open("cut_iter.mod", "w"))

    # 3. Applica il taglio e risolve di nuovo
    with open("addcut.run", "w") as f:
        f.write(f"""
model ufl_relax.mod;
data cap101.dat;
include cut_iter.mod;
option cplex_options 'cuts=0';
option solver cplex;
solve;
display TotalCost;
        """)
    subprocess.run(["ampl", "addcut.run"])
