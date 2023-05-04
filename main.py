'''
This Python file is for running each model

Group No: 2
Authors: Efe Erkan, Recep Uysal, Uygar Onat Erol
'''
import time
from rna_folding import RNA_Folding_MAX_PAIRS, RNA_Folding_MIN_Energy, RNA_Folding_MIN_Stack_Energy, RNA_Folding_MIN_Stack_Energy_Pseudoknots, RNA_Folding_MIN_Energy_DP

def main():
    
    RNA = "AAGCUCCUUGGUCCGGGCCAUUAAAGCCCGGAGCAGGAACAGAUGCUAUCCUAUUAAAGUAGGGUUACCC"

    # Part A
    timm1 = time.time()
    result1 = RNA_Folding_MAX_PAIRS(RNA, distance_limit=4)
    print(f"Max Pairs Group distance limit 4: {result1['Optimal_Result']}, Pairings:{result1['Pairings']}, time: {(time.time() - timm1) * 1000} \n")

    # Part B
    timm2 = time.time()
    result2 = RNA_Folding_MIN_Energy(RNA, distance_limit=4)
    print(f"Min Energy Group distance limit 4: {result2['Optimal_Result']}, Pairings:{result2['Pairings']}, time: {(time.time() - timm2) * 1000} \n")

    # Part C
    timm3 = time.time()
    result3 = RNA_Folding_MIN_Energy(RNA, distance_limit=7)
    print(f"Min Energy Group distance limit 7: {result3['Optimal_Result']}, Number of pairs: {result3['Pairings']}, time: {(time.time() - timm3) * 1000} \n")

    # Part D
    timm4 = time.time()
    result4 = RNA_Folding_MIN_Stack_Energy(RNA)
    print(f"Min Stack Energy Group: {result4['Optimal_Result']}, Pairings: {result4['Pairings']}, time: {(time.time() - timm4) * 1000} \n")

    # Part E
    timm5 = time.time()
    result5 = RNA_Folding_MIN_Stack_Energy_Pseudoknots(RNA)
    print(f"Min Energy Group with Pseudo-Knots: {result5['Optimal_Result']}, UpWard Pairings: {result5['Upward_Pairings']}, Downward Pairings: {result5['Downward_Pairings']}, time: {(time.time() - timm5) * 1000} \n")

    # Part F
    timm6 = time.time()
    result6 = RNA_Folding_MIN_Energy_DP(RNA, distance_limit=4)
    print(f"Dynamic Programming: {result6['Optimal_Result']}, Pairings:{result6['Pairings']}, time: {(time.time() - timm6) * 1000} \n")

if __name__ == '__main__':
    main()