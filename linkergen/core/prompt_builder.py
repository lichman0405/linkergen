from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import yaml
from typing import List

class PromptBuilder:
    def __init__(self, template_dir: str = None):
        if template_dir is None:
            # 自动定位当前脚本所在的父级目录下的 templates 文件夹
            current_dir = Path(__file__).resolve().parent.parent
            template_dir = current_dir / "templates"
        self.env = Environment(loader=FileSystemLoader(str(template_dir)))

    def load_examples_from_yaml(self, path: str) -> List[str]:
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
        return data.get("examples", [])

    def build_prompt(
        self,
        metal_sbu: str,
        n: int = 5,
        examples: List[str] = None,
        template_name: str = "basic_prompt.j2"
    ) -> str:
        if not examples:
            examples = [
                "OC(=O)c1ccc(cc1)C(O)=O",
                "OC(=O)c1cc(O)c(cc1O)C(O)=O",
                "OC(=O)CCC(=O)O"
            ]
        template = self.env.get_template(template_name)
        return template.render(metal_sbu=metal_sbu, num=n, examples=examples)


if __name__ == "__main__":
    #current_dir = Path(__file__).resolve().parent.parent
    #config_file = current_dir / "config" / "default_linkers.yaml"

    pb = PromptBuilder()  # 自动定位模板
    examples = pb.load_examples_from_yaml(r"linkergen\config\default_linkers.yaml")
    prompt = pb.build_prompt(metal_sbu="Zn4O", n=5, examples=examples)
    print(prompt)
