import os
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Dict

# Constants for fingerprinting
FP_LEN = 2048
FP_RAD = 2
SCORE_SCALE = 5.0

class SCScorer:
    def __init__(self, score_scale=SCORE_SCALE, model_path=None, verbose=False):
        self.vars = []
        self.score_scale = score_scale
        self._restored = False
        self.verbose = verbose
        self._initialize(model_path)

    def _initialize(self, model_path):
        if model_path is None:
            model_path = os.path.join(os.path.dirname(__file__), "models", "model.ckpt-10654.as_numpy.pickle")
        self.restore(model_path)

    def restore(self, model_path):
        with open(model_path, 'rb') as f:
            self.vars = pickle.load(f, encoding='latin1')
        self._restored = True
        if self.verbose:
            print(f"Restored model weights from {model_path}")
            for i, var in enumerate(self.vars):
                print(f"Var {i}: shape = {var.shape}")

        def mol_to_fp(mol):
            if mol is None:
                return np.zeros((FP_LEN,), dtype=bool)
            return np.array(AllChem.GetMorganFingerprintAsBitVect(
                mol, FP_RAD, nBits=FP_LEN, useChirality=True), dtype=bool)

        self.mol_to_fp = mol_to_fp

    def smi_to_fp(self, smi):
        if not smi:
            return np.zeros((FP_LEN,), dtype=bool)
        mol = Chem.MolFromSmiles(smi)
        return self.mol_to_fp(mol)

    def apply(self, x):
        if not self._restored:
            raise ValueError("Model weights must be restored before applying.")
        for i in range(0, len(self.vars), 2):
            W = self.vars[i]
            b = self.vars[i + 1]
            if x.shape[0] != W.shape[0]:
                raise ValueError(f"Input dim {x.shape[0]} ≠ weight dim {W.shape[0]}")
            x = np.matmul(x, W) + b
            if i < len(self.vars) - 2:
                x = np.maximum(x, 0)  # ReLU
        return 1 + (self.score_scale - 1) * (1 / (1 + np.exp(-x)))

    def get_score_from_smi(self, smi):
        fp = self.smi_to_fp(smi)
        score = self.apply(fp)
        return float(score.item())


class SynthesizabilityScorer:
    def __init__(self, model_path=None, verbose=False):
        self.model = SCScorer(model_path=model_path, verbose=verbose)

    def score(self, smiles_list: List[str]) -> Dict[str, float]:
        results = {}
        for smi in smiles_list:
            try:
                score = self.model.get_score_from_smi(smi)
                results[smi] = round(score, 3)
            except Exception:
                results[smi] = None
        return results


if __name__ == "__main__":
    scorer = SynthesizabilityScorer()
    smiles = ["CCO", "CCN", "C1=CC=CC=C1", "C1CCCC1", "C1=CCCCC1"]
    scores = scorer.score(smiles)
    for smi, score in scores.items():
        print(f"SMILES: {smi}, Score: {score}")