from gurobipy import Model, GRB, quicksum, tupledict
import gurobipy as gp

env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()

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

def RNA_Folding_MIN_Energy(RNA:str, A_U_energy:float = -1.33, G_C_ENERGY:float = -1.45, distance_limit:int = 4): # Part B
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
    return model.objVal
