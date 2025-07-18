import re
import os
 
def correggi_riga_set(riga):
    # Trova righe che iniziano con "set NOME := ..."
    match = re.match(r"^\s*set\s+(\w+)\s*:=\s*(.*);", riga)
    if not match:
        return riga
    nome_set = match.group(1)
    valori = match.group(2).strip().split()
    # Aggiungi apostrofi se il valore è numerico
    valori_str = [f"{v}" if v.isdigit() else v for v in valori]
    nuova_riga = f"set {nome_set} := {' '.join(valori_str)};\n"
    return nuova_riga
 
def correggi_file_dat(percorso_file):
    with open(percorso_file, "r") as f:
        righe = f.readlines()
 
    righe_corrette = [correggi_riga_set(r) for r in righe]
 
    # Salva backup
    os.rename(percorso_file, percorso_file + ".bak")
 
    with open(percorso_file, "w") as f:
        f.writelines(righe_corrette)
 
    print(f"✅ Corretto: {percorso_file} (backup in {percorso_file}.bak)")
 
def correggi_tutti_i_dat(nomi_file):
    for nome in nomi_file:
        if nome.endswith(".dat") and os.path.isfile(nome):
            correggi_file_dat(nome)
 
if __name__ == "__main__":
    # Inserisci qui i nomi dei file da correggere
    file_dat = [
        "cap101.dat",
        "cap102.dat",
        "cap103.dat",
        "cap104.dat",
        "cap71.dat",
        "cap72.dat",
        "cap73.dat",
        "cap74.dat"
    ]
    correggi_tutti_i_dat(file_dat)