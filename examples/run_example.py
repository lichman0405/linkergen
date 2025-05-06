import sys
from pathlib import Path

# 将项目根目录添加到 Python 模块路径中
sys.path.append(str(Path(__file__).resolve().parent.parent))

from linkergen.core.generator import LinkerGenPipeline

# 以下是原有调用逻辑
pipeline = LinkerGenPipeline()
results = pipeline.run(metal_sbu="Zn4O", n=5)

for item in results:
    print(item)
