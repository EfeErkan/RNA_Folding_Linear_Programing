'''
This Python file contains the Gurobi and Dynamic Programming Models

Group No: 2
Authors: Efe Erkan, Recep Uysal, Uygar Onat Erol
'''
from gurobipy import Model, GRB, quicksum, tupledict
import gurobipy as gp
import numpy as np

env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()

STACKING_PAIRS_ENERGY = {("A-U", "A-U"):-1.1, ("A-U", "C-G"):-2.1, ("A-U", "G-C"):-2.2, ("A-U", "U-A"):-0.6,
                         ("C-G", "A-U"):-2.1, ("C-G", "C-G"):-2.4, ("C-G", "G-C"):-3.3, ("C-G", "U-A"):-1.4,
                         ("G-C", "A-U"):-2.2, ("G-C", "C-G"):-3.3, ("G-C", "G-C"):-3.4, ("G-C", "U-A"):-1.5,
                         ("U-A", "A-U"):-0.6, ("U-A", "C-G"):-1.4, ("U-A", "G-C"):-1.5, ("U-A", "U-A"):-0.3}

def calculate_energy_pair(RNA:str, i:int, j:int):
    first_pair = f"{RNA[i - 1]}-{RNA[j - 1]}"
    second_pair = f"{RNA[i]}-{RNA[j - 2]}"
    if (first_pair, second_pair) in STACKING_PAIRS_ENERGY.keys():
        return STACKING_PAIRS_ENERGY[(first_pair, second_pair)]
    else:
        return 0

def isComplementary(ch1, ch2):
    if ch1 == 'A' and ch2 == 'U':
        return 1
    elif ch1 == 'U' and ch2 == 'A':
        return 1
    elif ch1 == 'G' and ch2 == 'C':
        return 2
    elif ch1 == 'C' and ch2 == 'G':
        return 2
    return 0

def isValidPairing(RNA:str, i:int, j:int, distance_limit:int = 4):
    if abs(i - j) < distance_limit:
        return False
    if isComplementary(RNA[i - 1], RNA[j - 1]) == 0:
        return False
    return True

def Energy(base1, base2):
    if (base1, base2) in [('A', 'U'), ('U', 'A')]:
        return -1.33
    elif (base1, base2) in [('G', 'C'), ('C', 'G')]:
        return -1.45
    else:
        return 0
    
def Optimal_Solution(optimized_model, var):
    solution = []
    for v in optimized_model.getVars():
        if v.VarName[0] == var and v.X == 1:
            index = v.VarName.index("_")
            solution.append((int(v.VarName[1:index]), int(v.VarName[index + 1:])))
    return solution

def RNA_Folding_MAX_PAIRS(RNA:str, distance_limit:int = 4): # Part A
    
    model = Model("RNA", env=env)
    X = tupledict()

    # Decision variables
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            X[i, j] = model.addVar(vtype=GRB.BINARY, name=f"x{i}_{j}")
    
    # Objective function   
    model.setObjective(X.sum(), GRB.MAXIMIZE)

    # At most 1 pairing        
    for i in range(1, len(RNA) + 1):
        model.addConstr(quicksum(X[j, i] for j in range(1, i)) + quicksum(X[i, j] for j in range(i + 1, len(RNA) + 1)) <= 1)

    # Compplementary characters => A-U, G-C and distance limitation
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            if isValidPairing(RNA, i, j, distance_limit) == False:
                model.addConstr(X[i, j] == 0)

    # No cross pairing            
    for i in range(1, len(RNA) + 1):
        for k in range(i + 1, len(RNA) + 1):
            for j in range(k + 1, len(RNA) + 1):
                for l in range(j + 1, len(RNA) + 1):
                    if isValidPairing(RNA, i, j, distance_limit) and isValidPairing(RNA, k, l, distance_limit):
                        model.addConstr(X[i, j] + X[k, l] <= 1)

    model.optimize()
    return {"Optimal_Result":model.objVal, "Pairings":Optimal_Solution(model, 'x')}

def RNA_Folding_MIN_Energy(RNA:str, A_U_energy:float = -1.33, G_C_ENERGY:float = -1.45, distance_limit:int = 4): # Part B and C
    model = Model("RNA", env=env)
    X = tupledict()

    # Decision variables
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            X[i, j] = model.addVar(vtype=GRB.BINARY, name=f"x{i}_{j}")
    
    # Objective function   
    model.setObjective(quicksum(X[i, j] * (A_U_energy if isComplementary(RNA[i - 1], RNA[j - 1]) == 1 else (G_C_ENERGY if isComplementary(RNA[i - 1], RNA[j - 1]) == 2 else 0)) for i in range(1, len(RNA) + 1) for j in range(i + 1, len(RNA) + 1)), GRB.MINIMIZE)

    # At most 1 pairing        
    for i in range(1, len(RNA) + 1):
        model.addConstr(quicksum(X[j, i] for j in range(1, i)) + quicksum(X[i, j] for j in range(i + 1, len(RNA) + 1)) <= 1)

    # Complementary characters => A-U, G-C and distance limitation
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            if isValidPairing(RNA, i, j, distance_limit) == False:
                model.addConstr(X[i, j] == 0)

    # No cross pairing            
    for i in range(1, len(RNA) + 1):
        for k in range(i + 1, len(RNA) + 1):
            for j in range(k + 1, len(RNA) + 1):
                for l in range(j + 1, len(RNA) + 1):
                    if isValidPairing(RNA, i, j, distance_limit) and isValidPairing(RNA, k, l, distance_limit):
                        model.addConstr(X[i, j] + X[k, l] <= 1)

    model.optimize()
    return {"Optimal_Result":model.objVal, "Pairings":Optimal_Solution(model, 'x')}

def RNA_Folding_MIN_Stack_Energy(RNA:str, distance_limit:int = 4): # Part D
    model = Model("RNA", env=env)
    X = tupledict()
    S = tupledict()

    # Decision variables
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            X[i, j] = model.addVar(vtype=GRB.BINARY, name=f"x{i}_{j}")
            S[i, j] = model.addVar(vtype=GRB.BINARY, name=f"s{i}_{j}")
    
    # Objective function   
    model.setObjective(quicksum(S[i, j] * calculate_energy_pair(RNA, i, j) for i in range(1, len(RNA) + 1) for j in range(i + 1, len(RNA) + 1)), GRB.MINIMIZE)

    # At most 1 pairing        
    for i in range(1, len(RNA) + 1):
        model.addConstr(quicksum(X[j, i] for j in range(1, i)) + quicksum(X[i, j] for j in range(i + 1, len(RNA) + 1)) <= 1)

    # Complementary characters => A-U, G-C and distance limitation
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            if isValidPairing(RNA, i, j, distance_limit) == False:
                model.addConstr(X[i, j] == 0)  

    # No cross pairing            
    for i in range(1, len(RNA) + 1):
        for k in range(i + 1, len(RNA) + 1):
            for j in range(k + 1, len(RNA) + 1):
                for l in range(j + 1, len(RNA) + 1):
                    if isValidPairing(RNA, i, j, distance_limit) and isValidPairing(RNA, k, l, distance_limit):
                        model.addConstr(X[i, j] + X[k, l] <= 1)
    
    # Stacked pairs check 1
    for i in range(1, len(RNA) + 1):
        for j in range(i + 3, len(RNA) + 1):
            model.addConstr(X[i, j] + X[i + 1, j - 1] - S[i, j] <= 1)

    # Stacked pairs check 2
    for i in range(1, len(RNA) + 1):
        for j in range(i + 3, len(RNA) + 1):
            model.addConstr(2 * S[i, j] - X[i, j] - X[i + 1, j - 1] <= 0)
               
    model.optimize()
    return {"Optimal_Result":model.objVal, "Pairings":Optimal_Solution(model, 'x')}

def RNA_Folding_MIN_Stack_Energy_Pseudoknots(RNA:str, distance_limit:int = 4): # Part E
    model = Model("RNA", env=env)
    X = tupledict()
    S = tupledict()
    Y = tupledict()
    Q = tupledict()

    # Decision variables
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            X[i, j] = model.addVar(vtype=GRB.BINARY, name=f"x{i}_{j}")
            S[i, j] = model.addVar(vtype=GRB.BINARY, name=f"s{i}_{j}")
            Y[i, j] = model.addVar(vtype=GRB.BINARY, name=f"y{i}_{j}")
            Q[i, j] = model.addVar(vtype=GRB.BINARY, name=f"q{i}_{j}")
    
    # Objective function   
    model.setObjective(quicksum(S[i, j] * calculate_energy_pair(RNA, i, j) + Q[i, j] * calculate_energy_pair(RNA, i, j) for i in range(1, len(RNA) + 1) for j in range(i + 1, len(RNA) + 1)), GRB.MINIMIZE)
    
    # At most 1 pairing        
    for i in range(1, len(RNA) + 1):
        model.addConstr(quicksum(X[j, i] for j in range(1, i)) + quicksum(X[i, j] for j in range(i + 1, len(RNA) + 1)) + quicksum(Y[j, i] for j in range(1, i)) + quicksum(Y[i, j] for j in range(i + 1, len(RNA) + 1)) <= 1)
    
    # Complementary characters => A-U, G-C and distance limitation
    for i in range(1, len(RNA) + 1):
        for j in range(i + 1, len(RNA) + 1):
            if isValidPairing(RNA, i, j, distance_limit) == False:
                model.addConstr(X[i, j] == 0)  
                model.addConstr(Y[i, j] == 0)  
    
    # No cross pairing            
    for i in range(1, len(RNA) + 1):
        for k in range(i + 1, len(RNA) + 1):
            for j in range(k + 1, len(RNA) + 1):
                for l in range(j + 1, len(RNA) + 1):
                    if isValidPairing(RNA, i, j, distance_limit) and isValidPairing(RNA, k, l, distance_limit):
                        model.addConstr(X[i, j] + X[k, l] <= 1)
                        model.addConstr(Y[i, j] + Y[k, l] <= 1)

    # Stacked pairs check for upward pairs
    for i in range(1, len(RNA) + 1):
        for j in range(i + 3, len(RNA) + 1):
            model.addConstr(X[i, j] + X[i + 1, j - 1] - S[i, j] <= 1)
            
    for i in range(1, len(RNA) + 1):
        for j in range(i + 3, len(RNA) + 1):
            model.addConstr(2 * S[i, j] - X[i, j] - X[i + 1, j - 1] <= 0)

    # Stacked pairs check for downward pairs
    for i in range(1, len(RNA) + 1):
        for j in range(i + 3, len(RNA) + 1):
            model.addConstr(Y[i, j] + Y[i + 1, j - 1] - Q[i, j] <= 1)
            
    for i in range(1, len(RNA) + 1):
        for j in range(i + 3, len(RNA) + 1):
            model.addConstr(2 * Q[i, j] - Y[i, j] - Y[i + 1, j - 1] <= 0)
            
    model.optimize()
    return {"Optimal_Result":model.objVal, "Upward_Pairings":Optimal_Solution(model, 'x'), "Downward_Pairings":Optimal_Solution(model, 'y')}

def RNA_Folding_MIN_Energy_DP(RNA:str, distance_limit:int = 4): # Part F
    n = len(RNA)
    W = np.zeros((n,n))
    
    P = np.empty((n,n), dtype=object)
    for i in range(n):
        for j in range(n):
            P[i,j] = []
    
    for k in range(distance_limit, n):
        for i in range(n - k):
            j = i + k
            min = W[i, j - 1]
            P[i, j] = P[i, j - 1]
            
            for t in range(i, j - distance_limit + 1):
                case = W[i, t - 1] + Energy(RNA[t], RNA[j]) + W[t + 1, j - 1]
                if (t + 2, j) in P[t + 1, j - 1]:
                    case += calculate_energy_pair(RNA, t + 1, j + 1)
                if case < min:
                    min = case
                    P[i, j] = P[i, t - 1] + [(t + 1, j + 1)] + P[t + 1, j - 1]
            W[i, j] = min
            
    return {"Optimal_Result":W[0, n - 1], "Pairings":P[0, n - 1]}