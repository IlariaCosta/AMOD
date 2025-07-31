import re

def process_data_file(input_filename, output_filename, max_facilities=16, max_clients=25):
    """
    Legge un file di dati in formato AMPL, restringe il numero di facilities e clients,
    e scrive il risultato in un nuovo file.

    Args:
        input_filename (str): Il percorso del file di dati di input.
        output_filename (str): Il percorso del file di dati di output.
        max_facilities (int): Il numero massimo di facilities da mantenere.
        max_clients (int): Il numero massimo di clients da mantenere.
    """
    with open(input_filename, 'r') as infile:
        content = infile.read()

    output_lines = []
    current_section = None
    
    # Pattern per le sezioni
    facilities_pattern = re.compile(r'set FACILITIES := (.+?);', re.DOTALL)
    clients_pattern = re.compile(r'set CLIENTS := (.+?);', re.DOTALL)
    param_f_pattern = re.compile(r'param f :=(.+?);', re.DOTALL)
    param_d_pattern = re.compile(r'param d :=(.+?);', re.DOTALL)
    param_c_pattern = re.compile(r'param c : (.+?) :=(.+?);', re.DOTALL)

    # Elabora la sezione FACILITIES
    match_facilities = facilities_pattern.search(content)
    if match_facilities:
        facilities_str = match_facilities.group(1).strip()
        facilities_list = facilities_str.split()
        truncated_facilities = facilities_list[:max_facilities]
        output_lines.append(f'set FACILITIES := {" ".join(truncated_facilities)};')
        content = facilities_pattern.sub('', content, 1) # Rimuove la sezione elaborata

    # Elabora la sezione CLIENTS
    match_clients = clients_pattern.search(content)
    if match_clients:
        clients_str = match_clients.group(1).strip()
        clients_list = clients_str.split()
        truncated_clients = clients_list[:max_clients]
        output_lines.append(f'set CLIENTS := {" ".join(truncated_clients)};')
        content = clients_pattern.sub('', content, 1) # Rimuove la sezione elaborata

    # Elabora la sezione param f
    match_param_f = param_f_pattern.search(content)
    if match_param_f:
        param_f_data = match_param_f.group(1).strip().split('\n')
        output_lines.append('\nparam f :=')
        for line in param_f_data:
            parts = line.strip().split()
            if parts and int(parts[0]) <= max_facilities:
                output_lines.append(f'  {line.strip()}')
        output_lines.append(';')
        content = param_f_pattern.sub('', content, 1)

    # Elabora la sezione param d
    match_param_d = param_d_pattern.search(content)
    if match_param_d:
        param_d_data = match_param_d.group(1).strip().split('\n')
        output_lines.append('\nparam d :=')
        for line in param_d_data:
            parts = line.strip().split()
            if parts and int(parts[0]) <= max_clients:
                output_lines.append(f'  {line.strip()}')
        output_lines.append(';')
        content = param_d_pattern.sub('', content, 1)

    # Elabora la sezione param c
    match_param_c = param_c_pattern.search(content)
    if match_param_c:
        c_header = match_param_c.group(1).strip()
        c_data = match_param_c.group(2).strip().split('\n')
        
        # Filtra l'header delle facilities se necessario, anche se dovrebbe essere già limitato
        c_header_parts = c_header.split()
        truncated_c_header_parts = [p for p in c_header_parts if int(p) <= max_facilities]
        output_lines.append(f'\nparam c : {" ".join(truncated_c_header_parts)} :=')
        
        for line in c_data:
            parts = line.strip().split()
            if parts and int(parts[0]) <= max_clients:
                client_id = parts[0]
                costs = parts[1:]
                truncated_costs = costs[:max_facilities]
                output_lines.append(f'{client_id} {" ".join(truncated_costs)}')
        output_lines.append(';')
        content = param_c_pattern.sub('', content, 1)

    # Scrivi il contenuto nel nuovo file
    with open(output_filename, 'w') as outfile:
        outfile.write('\n'.join(output_lines))
    
    print(f"File '{input_filename}' processato con successo. Output salvato in '{output_filename}'.")

if __name__ == "__main__":
    input_file = "cap41.dat"  # Nome del tuo file di input
    output_file = "cap41_restricted.dat" # Nome del file di output
    
    process_data_file(input_file, output_file, max_facilities=16, max_clients=25)