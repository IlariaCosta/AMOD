set FACILITIES;
set CLIENTS;
set CUTS;
param f {FACILITIES} >= 0;                # costo apertura facility
param d {CLIENTS} >= 0;                   # domanda cliente
param c {CLIENTS, FACILITIES} >= 0;      # costo trasporto cliente-facility
param capacity {FACILITIES} >= 0;         # capacità facility

var x {CLIENTS, FACILITIES} binary;       # assegnazione clienti (single source)
var y {FACILITIES} binary;                 # apertura facility

# Ogni cliente è assegnato a una sola facility
s.t. Assign {j in CLIENTS}:
    sum {i in FACILITIES} x[j,i] = 1;

# Non assegnare clienti a facility chiuse
s.t. OpenLink {i in FACILITIES, j in CLIENTS}:
    x[j,i] <= y[i];

# Capacità: domanda totale assegnata a facility i non supera capacità i
s.t. Capacity {i in FACILITIES}:
    sum {j in CLIENTS} d[j] * x[j,i] <= capacity[i] * y[i];

# Funzione obiettivo: costi apertura + costi trasporto
minimize TotalCost:
    sum {i in FACILITIES} f[i]*y[i] + sum {j in CLIENTS, i in FACILITIES} c[j,i]*x[j,i];
