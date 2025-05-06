from openai import OpenAI
from pathlib import Path
import yaml

class LLMClient:
    def __init__(self, config_path: str = None):
        if config_path is None:
            # 自动定位 secrets.yaml
            current_dir = Path(__file__).resolve().parent.parent
            config_path = current_dir / "config" / "secrets.yaml"
        with open(config_path, "r") as f:
            secrets = yaml.safe_load(f)

        self.api_key = secrets["deepseek"]["api_key"]
        self.base_url = secrets["deepseek"].get("base_url", "https://api.deepseek.com/v1")
        self.model = secrets["deepseek"].get("model", "deepseek-chat")

        self.client = OpenAI(api_key=self.api_key, base_url=self.base_url)

    def generate(self, prompt: str, system_message: str = "You are a MOF linker generator assistant.") -> str:
        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user", "content": prompt}
            ],
            stream=False
        )
        return response.choices[0].message.content.strip()

if __name__ == "__main__":
    from prompt_builder import PromptBuilder

    pb = PromptBuilder()
    examples = pb.load_examples_from_yaml(r"linkergen/config/default_linkers.yaml")
    prompt = pb.build_prompt(metal_sbu="Zn4O", n=5, examples=examples)

    llm = LLMClient()
    output = llm.generate(prompt)
    print(output)
