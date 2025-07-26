from pyomo.environ import *
import time

def estrai_matrici(instance):
    # Estrae matrice A, vettore b, vettore c dal modello Pyomo

    # Lista variabili (x e y)
    var_list = list(instance.component_data_objects(Var, descend_into=True))
    n_vars = len(var_list)

    # Lista constraints lineari (salvo solo quelli con corpo lineare)
    constr_list = []
    for c in instance.component_data_objects(Constraint, descend_into=True):
        if c.body.polynomial_degree() == 1:
            constr_list.append(c)
    n_cons = len(constr_list)

    # Matrice A (n_cons x n_vars), vettore b (rhs), vettore c (obiettivo)
    A = []
    b = []
    for c in constr_list:
        coefs = []
        # per ogni variabile calcoliamo coefficiente nel vincolo c
        for v in var_list:
            coef = linear_coef(c.body, v)
            coefs.append(coef)
        A.append(coefs)
        # rhs
        b.append(c.upper())  # supponendo vincolo <= b
    # c: coeff obiettivo
    c_obj = []
    for v in var_list:
        c_obj.append(instance.obj.expr.coeff(v))

    return A, b, c_obj, var_list

def linear_coef(expr, var):
    # Se expr è semplicemente la variabile stessa
    if expr == var:
        return 1.0
    
    # Se expr è una somma di termini
    if hasattr(expr, 'args'):
        coef = 0.0
        for term in expr.args:
            # Se il termine è tipo coef*var, cioè MulExpression
            # Pyomo rappresenta un termine lineare come MulExpression (coef * var)
            # Qui controlliamo se term è multiplo di var
            if hasattr(term, 'args'):
                # Cerca di trovare il var tra args
                if var in term.args:
                    # l'altro argomento è il coefficiente
                    for arg in term.args:
                        if arg != var:
                            # arg è coefficiente numerico o Pyomo Param
                            try:
                                coef += float(arg)
                            except:
                                pass
                elif term == var:
                    coef += 1.0
            elif term == var:
                coef += 1.0
        return coef

    # Se non trova la variabile
    return 0.0


def main():
    # Carica modello e dati
    model = AbstractModel()
    model_path = 'ufl.mod'
    data_path = 'cap71.dat'

    # Carico modello da file (pyomo script)
    # Se non vuoi riscrivere modello in python, puoi usare
    # pyomo.environ's 'Model' per importare file mod/dat (ma è complesso)
    # Qui per semplicità riscrivo modello in pyomo

    # Definizione modello Pyomo corrispondente (semplificata, basata su tuo modello)
    model.FACILITIES = Set()
    model.CLIENTS = Set()
    model.f = Param(model.FACILITIES, within=NonNegativeReals)
    model.d = Param(model.CLIENTS, within=NonNegativeReals)
    model.c = Param(model.CLIENTS, model.FACILITIES, within=NonNegativeReals)

    model.M = Param(initialize=1e5, mutable=True)  # Big M, potrai sovrascrivere

    model.x = Var(model.CLIENTS, model.FACILITIES, domain=NonNegativeReals)
    model.y = Var(model.FACILITIES, domain=Binary)

    def openlink_rule(model, i, j):
        return model.x[j, i] <= model.M * model.y[i]
    model.OpenLink = Constraint(model.FACILITIES, model.CLIENTS, rule=openlink_rule)

    def demand_rule(model, j):
        return sum(model.x[j, i] for i in model.FACILITIES) == model.d[j]
    model.demand = Constraint(model.CLIENTS, rule=demand_rule)

    def obj_rule(model):
        return sum(model.f[i] * model.y[i] for i in model.FACILITIES) + \
               sum(model.c[j, i] * model.x[j, i] for j in model.CLIENTS for i in model.FACILITIES)
    model.obj = Objective(rule=obj_rule, sense=minimize)

    instance = model.create_instance(data_path)

    # Calcolo big M (massimo domanda cliente)
    max_d = max(value(instance.d[j]) for j in instance.CLIENTS)
    instance.M = max_d

    # Risolvo
    solver = SolverFactory('glpk')  # Assicurati che GLPK sia installato
    start = time.time()
    results = solver.solve(instance)
    elapsed = time.time() - start

    # Estrazione e stampa
    print("Tempo risoluzione: {:.3f} sec".format(elapsed))

    # Stampo valori variabili y e x
    print("\nSoluzione variabili:")
    for i in instance.FACILITIES:
        print(f"y[{i}] = {value(instance.y[i])}")
    for j in instance.CLIENTS:
        for i in instance.FACILITIES:
            print(f"x[{j},{i}] = {value(instance.x[j,i])}")

    # Funzione obiettivo
    print("\nValore funzione obiettivo:", value(instance.obj))

if __name__ == '__main__':
    main()
