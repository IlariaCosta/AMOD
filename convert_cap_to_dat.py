def convert_cap_file_to_ampl_dat(input_file, output_file="output.dat"):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    m, n = map(int, lines[0].split())
    facilities = []
    idx = 1

    # Leggi capacità e costo fisso per ogni magazzino
    for _ in range(m):
        cap, cost = map(float, lines[idx].split())
        facilities.append((int(cap), int(cost)))
        idx += 1

    demands = []
    costs = []

    while len(demands) < n and idx < len(lines):
        # Leggi la domanda
        demand = float(lines[idx])
        demands.append(int(demand))
        idx += 1

        # Leggi i costi: m numeri in totale (anche su più righe)
        cost_values = []
        while len(cost_values) < m:
            cost_values += list(map(float, lines[idx].split()))
            idx += 1
        costs.append(cost_values)

    with open(output_file, 'w') as out:
        out.write("set I := " + " ".join(str(i + 1) for i in range(m)) + ";\n")
        out.write("set J := " + " ".join(str(j + 1) for j in range(n)) + ";\n\n")

        out.write("param f :=\n")
        for i, (_, fixed_cost) in enumerate(facilities):
            out.write(f"  {i+1} {fixed_cost}\n")
        out.write(";\n\n")
        
        out.write("param capacity :=\n")
        for i, (capacity, _) in enumerate(facilities):
            val = int(capacity) if (isinstance(capacity, float) and capacity.is_integer()) else capacity
            out.write(f"  {i+1} {val}\n")

        out.write(";\n\n")
        
        out.write("param d :=\n")
        for j, demand in enumerate(demands):
            out.write(f"  {j+1} {demand}\n")
        out.write(";\n\n")

        out.write("param c : " + " ".join(str(i + 1) for i in range(m)) + " :=\n")
        for j, row in enumerate(costs):
            out.write(f"{j+1} " + " ".join(str(int(c) if c.is_integer() else str(c)) for c in row) + "\n")
        out.write(";\n")

    print(f"✔ Conversione completata: {output_file}")

# Esegui direttamente
convert_cap_file_to_ampl_dat("cap41.txt", "cap41.dat")
convert_cap_file_to_ampl_dat("cap42.txt", "cap42.dat")
convert_cap_file_to_ampl_dat("cap43.txt", "cap43.dat")
# convert_cap_file_to_ampl_dat("cap74.txt", "cap74.dat")
# convert_cap_file_to_ampl_dat("cap101.txt", "cap101.dat")
# convert_cap_file_to_ampl_dat("cap102.txt", "cap102.dat")
# convert_cap_file_to_ampl_dat("cap103.txt", "cap103.dat")
# convert_cap_file_to_ampl_dat("cap104.txt", "cap104.dat")
