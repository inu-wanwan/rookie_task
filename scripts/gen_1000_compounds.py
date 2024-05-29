import os
from moses.vae import VAE
from rdkit import Chem

# モデルとボキャブラリのパス
model_path = 'workspace/model.pt'
vocab_path = 'workspace/vocab.txt'
seeds_path = 'workspace/data/seeds.smi'
output_path = 'workspace/gen/generated_compounds.smi'

# VAEモデルのロード
vae = VAE()
vae.load(model_path, vocab_path)

# シード化合物の読み込み
with open(seeds_path, 'r') as f:
    seeds = [line.strip() for line in f]

# 各シードから100個の化合物を生成
generated_compounds = []
for seed in seeds:
    seed_mol = Chem.MolFromSmiles(seed)
    if seed_mol is None:
        continue
    generated = vae.sample(100, seed=seed)
    generated_compounds.extend(generated)

# 生成された化合物をファイルに書き出し
with open(output_path, 'w') as f:
    for smi in generated_compounds:
        f.write(smi + '\n')

print(f'Generated {len(generated_compounds)} compounds.')
