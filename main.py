import pandas as pd
from rna_folding import RNA_Folding_MAX_PAIRS

def main():
    RNA_df = pd.read_excel('data.xlsx')
    
    for index, row in RNA_df.iterrows():
        RNA = row['RNA']
        print(f"Max pairs Group {index + 1}: {RNA_Folding_MAX_PAIRS(RNA)}")

if __name__ == '__main__':
    main()