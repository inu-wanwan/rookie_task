import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_sdf(smiles_file, sdf_file):
    # SMILESファイルの読み込み
    suppl = Chem.SmilesMolSupplier(smiles_file, titleLine=False)
    
    # SDFファイルの書き込み準備
    writer = Chem.SDWriter(sdf_file)
    
    for mol in suppl:
        if mol is None:
            continue
        # 3D座標の生成
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        # SDFファイルに分子を書き込む
        writer.write(mol)
    
    writer.close()
    print(f"SMILESファイル {smiles_file} をSDFファイル {sdf_file} に変換しました。")

def main():
    parser = argparse.ArgumentParser(description="Convert SMILES file to SDF file")
    parser.add_argument("smiles_file", type=str, help="Path to the input SMILES file")
    parser.add_argument("sdf_file", type=str, help="Path to the output SDF file")
    
    args = parser.parse_args()
    
    smiles_to_sdf(args.smiles_file, args.sdf_file)

if __name__ == "__main__":
    main()
