import sys

def load_y_values(file_path):
    """
    Legge il file con i valori delle variabili y[i] (apertura facility)
    Restituisce un dizionario: { 'F1': 0.35, 'F2': 0.92, ... }
    """
    y_vals = {}
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue  # ignora righe malformate
            facility, val = parts
            try:
                y_vals[facility] = float(val)
            except ValueError:
                continue  # ignora valori non numerici
    return y_vals

def is_number(s):
    """
    Funzione helper per controllare se una stringa è un numero (utile per il formatting)
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def find_most_fractional(y_vals, num=3, tol=1e-5):
    """
    Trova le variabili y[i] più frazionarie, cioè più lontane da 0 o 1.
    Restituisce una lista con i nomi delle facility più critiche (es. ['F1', 'F3']).
    """
    fractional = []
    for i, y in y_vals.items():
        frac_part = abs(y - round(y))
        if frac_part > tol:  # se la variabile non è intera
            fractional.append((i, frac_part))
    # Ordina in ordine decrescente di "frazionarietà"
    fractional.sort(key=lambda x: -x[1])
    return [i for i, _ in fractional[:num]]  # restituisce i top `num` facility

def generate_cut(facilities, iter_num=1):
    """
    Genera il vincolo AMPL per forzare l'apertura di almeno una delle facility date.
    Es: y['F1'] + y['F2'] + y['F3'] >= 1;
    """
    if not facilities:
        return f"# Nessuna variabile frazionaria trovata, nessun taglio generato nella iterazione {iter_num}.\n"
    
    cut = f"s.t. GomoryCut{iter_num}: "
    cut += " + ".join([f"y[{f}]" if is_number(f) else f"y['{f}']" for f in facilities])
    cut += " >= 1;\n"
    return cut

if __name__ == "__main__":
    # Controllo che venga passato almeno il file y_sol.txt
    if len(sys.argv) < 2:
        print("Uso: python build_cut.py y_sol.txt [iter_num]")
        sys.exit(1)

    y_file = sys.argv[1]  # Nome del file y_sol.txt passato da AMPL
    iter_num = int(sys.argv[2]) if len(sys.argv) >= 3 else 1  # Numero iterazione opzionale

    # Carica i valori di y[i]
    y_vals = load_y_values(y_file)

    # Trova le 3 facility più "frazionarie"
    most_frac = find_most_fractional(y_vals, num=3)

    # Genera il taglio e stampalo su stdout (che AMPL catturerà)
    cut_str = generate_cut(most_frac, iter_num)
    print(cut_str)
