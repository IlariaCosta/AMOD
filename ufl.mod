# modello UFL - Uncapacitated Facility Location

set FACILITIES;
set CLIENTS;

param f {FACILITIES} >= 0;                # costo apertura facility
param d {CLIENTS} >= 0;                   # domanda cliente (non usata nel modello ma puoi aggiungere)
param c {CLIENTS, FACILITIES} >= 0;      # costo trasporto cliente-facility

var x {CLIENTS, FACILITIES} >= 0, <= 1;  # assegnazione clienti (frazi o intera)
var y {FACILITIES} binary;                # apertura facility (variabile intera)

# Ogni cliente è assegnato a una e una sola facility
s.t. Assign {j in CLIENTS}:
    sum{i in FACILITIES} x[j,i] = 1;

# Non assegnare clienti a facility chiuse
s.t. OpenLink {i in FACILITIES, j in CLIENTS}:
    x[j,i] <= y[i];

# Funzione obiettivo: costi apertura + costi trasporto
minimize TotalCost:
    sum{i in FACILITIES} f[i]*y[i] + sum{j in CLIENTS, i in FACILITIES} c[j,i]*x[j,i];
