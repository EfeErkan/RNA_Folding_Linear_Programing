import pandas as pd
from rna_folding import RNA_Folding_MAX_PAIRS, RNA_Folding_MIN_Energy

def main():
    RNA_df = pd.read_excel('data.xlsx')
    
    for index, row in RNA_df.iterrows():
        RNA = row['RNA']
        print(f"Min energy Group {index + 1}: {RNA_Folding_MIN_Energy(RNA)}")

if __name__ == '__main__':
    main()