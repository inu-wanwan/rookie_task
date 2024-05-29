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

def main(input_file, output_image):
    with open(input_file, 'r') as file:
        smiles_list = [line.strip() for line in file.readlines()]

    molecular_weights, qed_values, logp_values = calculate_properties(smiles_list)

    data = {
        'SMILES': smiles_list,
        'MolecularWeight': molecular_weights,
        'QED': qed_values,
        'MolLogP': logp_values
    }
    df = pd.DataFrame(data)

    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.hist(df['MolecularWeight'], bins=10, color='b', edgecolor='k', alpha=0.7)
    plt.title('Molecular Weight Distribution')
    plt.xlabel('Molecular Weight')
    plt.ylabel('Frequency')

    plt.subplot(1, 3, 2)
    plt.hist(df['QED'], bins=10, color='g', edgecolor='k', alpha=0.7)
    plt.title('QED Distribution')
    plt.xlabel('QED')
    plt.ylabel('Frequency')

    plt.subplot(1, 3, 3)
    plt.hist(df['MolLogP'], bins=10, color='r', edgecolor='k', alpha=0.7)
    plt.title('MolLogP Distribution')
    plt.xlabel('MolLogP')
    plt.ylabel('Frequency')

    plt.tight_layout()
    plt.savefig(output_image)
    plt.show()

    print(f"Plot saved as {output_image}")

    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate molecular properties from SMILES and plot distributions.")
    parser.add_argument("input_file", help="Input file containing SMILES")
    parser.add_argument("output_image", help="Output image file to save the plot")
    args = parser.parse_args()

    df = main(args.input_file, args.output_image)
    print(df)
