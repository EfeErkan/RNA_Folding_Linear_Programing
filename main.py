import pandas as pd
from rna_folding import RNA_Folding_MAX_PAIRS, RNA_Folding_MIN_Energy

def main():
    RNA_df = pd.read_excel('data.xlsx')
    
    for index, row in RNA_df.iterrows():
        RNA = row['RNA']
        result1 = RNA_Folding_MIN_Energy(RNA)
        result2 = RNA_Folding_MIN_Energy(RNA, distance_limit=7)
        print(f"Min energy Group {index + 1} distance limit 4: {result1[0]}, Number of pairs: {result1[1]}")
        print(f"Min energy Group {index + 1} distance limit 7: {result2[0]}, Number of pairs: {result2[1]}")

if __name__ == '__main__':
    main()