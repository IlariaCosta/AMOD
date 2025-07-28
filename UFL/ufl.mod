# modello UFL - Uncapacitated Facility Location

set FACILITIES;
set CLIENTS;

param f {FACILITIES} >= 0;                # costo apertura facility
param d {CLIENTS} >= 0;                   # domanda cliente (non usata nel modello ma puoi aggiungere)
param c {CLIENTS, FACILITIES}  >= 0;      # costo trasporto cliente-facility

param M := max {j in CLIENTS} d[j];						# definisco big M
#var x {CLIENTS, FACILITIES} binary;  # assegnazione clienti (frazi o intera)
var x {FACILITIES,CLIENTS}  >=0;  # assegnazione clienti (frazi o intera)
var y {FACILITIES} binary;                # apertura facility (variabile intera)


/*# Ogni cliente è assegnato a una e una sola facility
s.t. Assign {j in CLIENTS}:
    sum{i in FACILITIES} x[j,i] = 1;
*/

# La domanda del cliente deve essere soddisfatta
s.t. demand {j in CLIENTS}:
	sum{i in FACILITIES}  x[i,j] = d[j];
	
# Non assegnare clienti a facility chiuse
s.t. OpenLink {i in FACILITIES, j in CLIENTS}:
     x[i,j] <= M*y[i];


	
# Funzione obiettivo: costi apertura + costi trasporto
minimize TotalCost:
    sum{i in FACILITIES} f[i]*y[i] + sum{i in FACILITIES, j in CLIENTS} c[j,i]* x[i,j];
