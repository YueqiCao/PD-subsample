'''
Train Poincare embeddings of the WordNet Mammal relations.

Install the additional training dependency in the project environment:
    python -m pip install gensim==4.4.0

Run from the repository directory:
    python scripts/train_poincare.py --dimension 2 --epochs 20 50 150 200

To train the other dimensions used in the original notebooks:
    python scripts/train_poincare.py --dimension 5 --epochs 20 50 150 200
'''

import argparse
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/poincare"
SEED = 20211121


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dimension", type=int, default=2)
    parser.add_argument("--epochs", type=int, nargs="+", default=[20, 50, 150, 200])
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/poincare_training")
    args = parser.parse_args()
    from gensim.models.poincare import PoincareModel, PoincareRelations

    args.output.mkdir(parents=True, exist_ok=True)
    words = np.load(DATA / "words.npy")
    word_path = args.output / "words.npy"
    np.save(word_path, words)
    relations = list(PoincareRelations(str(DATA / "wordnet_mammal_hypernyms.tsv"), delimiter="\t"))
    for epochs in args.epochs:
        destination = args.output / f"embedding_dim{args.dimension}_epoch{epochs}.npy"
        model = PoincareModel(train_data=relations, size=args.dimension, burn_in=10, seed=SEED, workers=1)
        model.train(epochs=epochs, print_every=500)
        vectors = np.array([model.kv[word] for word in words])
        np.save(destination, vectors)
        print(f"dimension={args.dimension}, epochs={epochs}: vectors saved to {destination}", flush=True)


if __name__ == "__main__":
    main()
