from pathlib import Path
from typing import List, Dict

from linkergen.core.prompt_builder import PromptBuilder
from linkergen.core.llm_interface import LLMClient
from linkergen.core.smiles_validator import parse_smiles_from_text, validate_smiles
from linkergen.core.synthesizability import SynthesizabilityScorer
from linkergen.core.formula_combiner import FormulaCombiner

class LinkerGenPipeline:
    def __init__(
        self,
        template_dir: str = None,
        config_path: str = None,
        model_path: str = None,
        linker_count: int = 3
    ):
        self.prompt_builder = PromptBuilder(template_dir=template_dir)
        self.llm_client = LLMClient(config_path=config_path)
        self.validator = SynthesizabilityScorer(model_path=model_path)
        self.combiner = FormulaCombiner(linker_count=linker_count)

    def run(self, metal_sbu: str, n: int = 5, example_yaml: str = None) -> List[Dict]:
        # 1. Load example linkers
        if example_yaml is None:
            root = Path(__file__).resolve().parent.parent
            example_yaml = root / "config" / "default_linkers.yaml"
        examples = self.prompt_builder.load_examples_from_yaml(str(example_yaml))

        # 2. Build prompt
        prompt = self.prompt_builder.build_prompt(metal_sbu=metal_sbu, n=n, examples=examples)

        # 3. Call LLM
        raw_output = self.llm_client.generate(prompt)

        # 4. Parse & validate SMILES
        candidates = parse_smiles_from_text(raw_output)
        valid_smiles = validate_smiles(candidates)

        # 5. Score with SCScore
        scored = self.validator.score(valid_smiles)

        # 6. Combine with metal SBU
        return self.combiner.combine(metal_sbu, scored)

