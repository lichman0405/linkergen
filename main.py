import argparse
from linkergen.core.generator import LinkerGenPipeline

def main():
    parser = argparse.ArgumentParser(description="LinkerGen: MOF linker generator pipeline")
    parser.add_argument("--metal", type=str, required=True, help="Target metal SBU (e.g. Zn4O, Cu2Paddlewheel)")
    parser.add_argument("--n", type=int, default=5, help="Number of linkers to generate")
    parser.add_argument("--linker_count", type=int, default=3, help="Number of linker units per MOF formula (e.g. Zn4O + 3x linker)")
    parser.add_argument("--template_dir", type=str, default=None, help="Path to Jinja2 templates")
    parser.add_argument("--config_path", type=str, default=None, help="Path to secrets.yaml")
    parser.add_argument("--model_path", type=str, default=None, help="Path to SCScore model file")
    parser.add_argument("--examples", type=str, default=None, help="Path to default_linkers.yaml")

    args = parser.parse_args()

    pipeline = LinkerGenPipeline(
        template_dir=args.template_dir,
        config_path=args.config_path,
        model_path=args.model_path,
        linker_count=args.linker_count
    )

    results = pipeline.run(
        metal_sbu=args.metal,
        n=args.n,
        example_yaml=args.examples
    )

    print("\nGenerated Linkers:\n")
    for item in results:
        print(f"- [Formula] {item['formula']}")
        print(f"  SMILES: {item['linker']}")
        print(f"  SCScore: {item['scscore']}")
        print("")

if __name__ == "__main__":
    main()
