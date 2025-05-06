from typing import List, Dict

class FormulaCombiner:
    def __init__(self, linker_count: int = 3):
        self.linker_count = linker_count 

    def combine(self, metal_sbu: str, smiles_score_map: Dict[str, float]) -> List[Dict]:
        results = []
        for smiles, score in smiles_score_map.items():
            if score is None:
                continue  # 跳过不可评分结构
            formula = f"{metal_sbu} + {self.linker_count}x {smiles}"
            results.append({
                "metal_sbu": metal_sbu,
                "linker": smiles,
                "scscore": score,
                "formula": formula
            })
        return results

if __name__ == "__main__":
    metal = "Zn4O"
    scored_smiles = {
        "OC(=O)c1ccc(cc1)C(=O)O": 2.13,
        "OC(=O)CC(C(=O)O)C1=CC=CC=C1": 2.45
    }

    combiner = FormulaCombiner(linker_count=3)
    results = combiner.combine(metal, scored_smiles)

    for r in results:
        print(r)