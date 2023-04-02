from gurobipy import Model, GRB, quicksum, tupledict
import gurobipy as gp

env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()

STACKING_PAIRS_ENERGY = {("A-U", "A-U"):-1.1, ("A-U", "C-G"):-2.1, ("A-U", "G-C"):-2.2, ("A-U", "U-A"):-0.6,
                         ("C-G", "A-U"):-2.1, ("C-G", "C-G"):-2.4, ("C-G", "G-C"):-3.3, ("C-G", "U-A"):-1.4,
                         ("G-C", "A-U"):-2.2, ("G-C", "C-G"):-3.3, ("G-C", "G-C"):-3.4, ("G-C", "U-A"):-1.5,
                         ("U-A", "A-U"):-0.6, ("U-A", "C-G"):-1.4, ("U-A", "G-C"):-1.5, ("U-A", "U-A"):-0.3}

def calculate_energy_pair(RNA:str, i:int, j:int):
    return STACKING_PAIRS_ENERGY[(f"{RNA[i]}-{RNA[j]}", f"{RNA[i + 1]}-{RNA[j - 1]}")]

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
    return model.objVal

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
    num_of_pairs = quicksum(X[i, j] for i in range(1, len(RNA) + 1) for j in range(i + 1, len(RNA) + 1)).getValue()
    return model.objVal, num_of_pairs

def RNA_Folding_MIN_Stack_Energy(RNA:str, distance_limit:int = 4): # Part D
    pass