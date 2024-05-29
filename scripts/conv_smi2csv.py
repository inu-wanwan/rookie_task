import pandas as pd
from rdkit import Chem
import argparse

def convert_smi_to_csv(smi_file, csv_file):
    # SMILESファイルを読み込む
    smiles_list = []
    with open(smi_file, 'r') as file:
        for line in file:
            smiles = line.strip()
            smiles_list.append(smiles)
    
    # DataFrameに変換
    df = pd.DataFrame(smiles_list, columns=['SMILES'])
    
    # CSVファイルとして保存
    df.to_csv(csv_file, index=False, encoding='utf-8')

def main():
    parser = argparse.ArgumentParser(description='Convert .smi file to .csv file.')
    parser.add_argument('input_file', type=str, help='Path to the input .smi file')
    parser.add_argument('output_file', type=str, help='Path to the output .csv file')
    args = parser.parse_args()
    
    convert_smi_to_csv(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
