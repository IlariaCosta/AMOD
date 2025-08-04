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
    #ampl.set_option('cplex_options', 'display=2 simplex display=2')
    ampl.set_option('solver', 'cplexamp')
    ampl.set_option('display', 1)

    already_cut = set()

    #if all_cuts:
    ampl.set_option('gomory_cuts', -1)  # all available
    ampl.set_option('solver_msg', 0)
    t0 = time.time()
    ampl.solve()
    elapsed = time.time() - t0
    var_values = list(ampl.get_variable('y').get_values().to_dict().values())
    #print(var_values)
    
    obj = ampl.obj['TotalCost'].value()
    return obj, elapsed, 1
    
