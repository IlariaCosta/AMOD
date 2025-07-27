

def is_integral(var_values):
    return all(abs(val - round(val)) < 1e-8 for val in var_values)
 
def compute_gap(relaxed_val, optimal_val):
    return abs((relaxed_val - optimal_val) / optimal_val) * 100

def invert_matrix(m):
    n = len(m)
    # crea matrice identità
    identity = [[float(i == j) for j in range(n)] for i in range(n)]
    
    # copia della matrice per non modificare l'originale
    mat = [row[:] for row in m]

    for i in range(n):
        # trova pivot
        pivot = mat[i][i]
        if abs(pivot) < 1e-15:
            # cerca righe da scambiare
            for r in range(i+1, n):
                if abs(mat[r][i]) > 1e-15:
                    mat[i], mat[r] = mat[r], mat[i]
                    identity[i], identity[r] = identity[r], identity[i]
                    pivot = mat[i][i]
                    break
            else:
                raise ValueError("Matrice singolare, non invertibile")

        # normalizza riga pivot
        for j in range(n):
            mat[i][j] /= pivot
            identity[i][j] /= pivot

        # elimina altre righe
        for r in range(n):
            if r != i:
                factor = mat[r][i]
                for c in range(n):
                    mat[r][c] -= factor * mat[i][c]
                    identity[r][c] -= factor * identity[i][c]

    return identity

def mat_vec_mul(mat, vec):
    result = []
    for row in mat:
        s = 0.0
        for a, b in zip(row, vec):
            s += a * b
        result.append(s)
    return result

def mat_mul(A, B):
    rows = len(A)
    cols = len(B[0])
    inner = len(B)
    result = [[0]*cols for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            s = 0
            for k in range(inner):
                s += A[i][k] * B[k][j]
            result[i][j] = s
    return result


def parse_dat_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    facilities = []
    clients = []
    f_param = {}
    c_param = []
    demands_dict = {}
    reading = None

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("set FACILITIES"):
            facilities = list(map(int, line.split(":=")[1].replace(";", "").split()))
            continue
        elif line.startswith("set CLIENTS"):
            clients = list(map(int, line.split(":=")[1].replace(";", "").split()))
            continue
        elif line.startswith("param f"):
            reading = "f"
            continue
        elif line.startswith("param c"):
            reading = "c"
            c_param = []
            continue
        elif line.startswith("param d"):
            reading = "d"  # ignoriamo
            continue

        # Lettura dei costi f
        if reading == "f":
            if ";" in line:
                line = line.replace(";", "")
                reading = None
            parts = line.split()
            if len(parts) == 2:
                f_param[int(parts[0])] = float(parts[1])  

        # Lettura della matrice dei costi c
        elif reading == "c":
            if ":=" in line:
                continue
            if ";" in line:
                line = line.replace(";", "")
                if line.strip():
                    c_param.append(list(map(float, line.split())))  
                reading = None
            else:
                if line.strip():
                    c_param.append(list(map(float, line.split()))) 

        # Lettura delle domande d
        elif reading == "d":
            if ";" in line:
                line = line.replace(";", "")
                reading = None
            parts = line.split()
            if len(parts) == 2:
                demands_dict[int(parts[0])] = float(parts[1])

    # Ordina f_param secondo l'ordine dei facilities
    f_vector = [f_param[i] for i in sorted(facilities)]
    demands = [demands_dict[i] for i in sorted(demands_dict.keys())]
    return facilities, clients, f_vector, c_param, demands

def build_A_b(facilities, clients, demands):
    m = len(facilities)
    n = len(clients)
    slack = n*m;
    var_count = m * n + m + slack
    row_count = n + m * n
    print("lunghezza domanda")
    
    A = [[0 for _ in range(var_count)] for _ in range(row_count)]
    # print("dimensioni A")
    # print(len(A), len(A[0]))
    b = [0 for _ in range(row_count)]

    def x_index(i, j):
        return i * n + j

    def y_index(i):
        return m * n + i
    
    def slack_index(i,j):
        return m * n + m + (i*n +j)  # slack iniziano dopo x e y
        
    # Cliente j deve essere assegnato a una sola facility
    for j in range(n):
        for i in range(m):
            A[j][x_index(i, j)] = 1
        b[j] = demands[j]

    # x_ij ≤ y_i
    row = n
    
    
    for i in range(m):      #ciente
        for j in range(n):  #facil
            A[row][x_index(i, j)] = 1
            A[row][y_index(i)] = -1
            A[row][slack_index(i,j)] = 1  # coefficiente slack +1
            b[row] = 0
            row += 1

    print("A e b pronte\n")
    # print(b)
    # # Stampa matrice A
    # print("Matrice A:")
    # for row in A:
    #     print(row)

    # # Stampa vettore b
    # print("\nVettore b:")
    # print(b)

    return A, b



