import pandas as pd
from rna_folding import RNA_Folding_MAX_PAIRS, RNA_Folding_MIN_Energy, RNA_Folding_MIN_Stack_Energy, RNA_Folding_MIN_Stack_Energy_Pseudoknots, RNA_Folding_MIN_Energy_DP

def main():
    
    RNA = "AAGCUCCUUGGUCCGGGCCAUUAAAGCCCGGAGCAGGAACAGAUGCUAUCCUAUUAAAGUAGGGUUACCC"

    # result1 = RNA_Folding_MIN_Energy(RNA, distance_limit=4)
    # result2 = RNA_Folding_MIN_Energy(RNA, distance_limit=7)
    # result3 = RNA_Folding_MIN_Stack_Energy(RNA)
    result4 = RNA_Folding_MIN_Stack_Energy_Pseudoknots(RNA)
    # result5 = RNA_Folding_MIN_Energy_DP(RNA, distance_limit=4)

    # print(f"Min Energy Group {index + 1} distance limit 4: {result1['Optimal_Result']}, Pairings:{result1['Pairings']}")
    # print(f"Min Energy Group {index + 1} distance limit 7: {result2[0]}, Number of pairs: {result2[1]}")
    # print(f"Min Stack Energy Group {index + 1}: {result3['Optimal_Result']}, Pairings: {result3['Pairings']}")
    print(f"Min Energy Group with Pseudo-Knots {index + 1}: {result4['Optimal_Result']}, UpWard Pairings: {result4['Upward_Pairings']}, Downward Pairings: {result4['Downward_Pairings']}")
    # print(f"Dynamic Programming: {result5['Optimal_Result']}, Pairings:{result5['Pairings']}")

if __name__ == '__main__':
    main()