from rdkit import Chem
from typing import List

def parse_smiles_from_text(text: str) -> List[str]:
    lines = text.strip().splitlines()
    smiles_candidates = [line.strip() for line in lines if line.strip()]
    return smiles_candidates

def validate_smiles(smiles_list: List[str]) -> List[str]:
    valid_smiles = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid_smiles.append(smi)
    return valid_smiles

if __name__ == "__main__":
    from prompt_builder import PromptBuilder
    from llm_interface import LLMClient

    # 构建 prompt
    pb = PromptBuilder()
    examples = pb.load_examples_from_yaml(r"linkergen/config/default_linkers.yaml")
    prompt = pb.build_prompt(metal_sbu="Zn4O", n=5, examples=examples)

    # 调用 LLM
    llm = LLMClient()
    raw_output = llm.generate(prompt)

    # 提取 + 校验 SMILES
    smiles_candidates = parse_smiles_from_text(raw_output)
    valid_smiles = validate_smiles(smiles_candidates)
    print("Valid SMILES:")
    print(valid_smiles)
    print("Invalid SMILES:")
    print(set(smiles_candidates) - set(valid_smiles))