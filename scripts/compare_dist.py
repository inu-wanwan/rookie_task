import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors
import argparse

def calculate_properties(smiles_list):
    molecular_weights = []
    qed_values = []
    logp_values = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecular_weights.append(Descriptors.MolWt(mol))
            qed_values.append(QED.qed(mol))
            logp_values.append(Crippen.MolLogP(mol))
    
    return molecular_weights, qed_values, logp_values

def make_dataframe(file):
    with open(file, 'r') as file:
        gen_smiles_list = [line.strip() for line in file.readlines()]

    molecular_weights, qed_values, logp_values = calculate_properties(gen_smiles_list)

    data = {
        'SMILES': gen_smiles_list,
        'MolecularWeight': molecular_weights,
        'QED': qed_values,
        'MolLogP': logp_values
    }
    return pd.DataFrame(data)

def main(gen_file, train_file, output_image):

    gen_df = make_dataframe(gen_file)
    train_df = make_dataframe(train_file)

    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.hist([gen_df['MolecularWeight'], train_df['MolecularWeight']], bins=10, density=True, color=['b', 'r'], edgecolor='b', alpha=0.7, label=['Generated Data', 'Train Data'])
    plt.title('Molecular Weight Distribution')
    plt.xlabel('Molecular Weight')
    plt.ylabel('Frequency')

    plt.subplot(1, 3, 2)
    plt.hist([gen_df['QED'], train_df['QED']], bins=10, density=True, color=['b', 'r'], edgecolor='b', alpha=0.7, label=['Generated Data', 'Train Data'])
    plt.title('QED Distribution')
    plt.xlabel('QED')
    plt.ylabel('Frequency')

    plt.subplot(1, 3, 3)
    plt.hist([gen_df['MolLogP'], train_df['MolLogP']], bins=10, density=True, color=['b', 'r'], edgecolor='b', alpha=0.7, label=['Generated Data', 'Train Data'])
    plt.title('MolLogP Distribution')
    plt.xlabel('MolLogP')
    plt.ylabel('Frequency')

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_image)
    plt.show()

    print(f"Plot saved as {output_image}")

    return gen_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate molecular properties from SMILES and plot distributions.")
    parser.add_argument("input_file_1", help="Input file containing SMILES")
    parser.add_argument("input_file_2", help="Input file containing SMILES")
    parser.add_argument("output_image", help="Output image file to save the plot")
    args = parser.parse_args()

    df = main(args.input_file_1, args.input_file_2, args.output_image)
    print(df)
