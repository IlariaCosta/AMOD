from amplpy import AMPL

def relax_all_integer_vars(ampl):
    vars_map = ampl.get_variables()
    for var_name in vars_map:
        var = vars_map[var_name]
        if var.is_integer():
            var.set_integer(False)

# Uso esempio
ampl = AMPL()
ampl.read("ufl.mod")
ampl.read_data("cap101.dat")

relax_all_integer_vars(ampl)  # Rilassa tutte le variabili intere

ampl.solve()

# Stampa i valori della variabile y
y_vals = ampl.get_variable('y').get_values()
print(y_vals)
