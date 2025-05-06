# LinkerGen：用于 MOF 连接体生成的模块化流水线

**LinkerGen** 是一个用于金属有机框架（MOFs）中有机连接体（linkers）设计与筛选的工程化工具。它结合了大型语言模型（LLM）、RDKit 和 SCScore，可用于自动生成 SMILES 表达式、验证结构有效性、评估合成复杂性，并输出结构化的连接体信息。

---

## 功能特点

- 基于 few-shot 提示构建，支持 LLM（如 OpenAI GPT-4）生成连接体
- 使用 RDKit 检查 SMILES 表达式的结构合法性
- 使用 SCScore（NumPy 重写版）评估分子合成复杂性, 有兴趣可以直接看[scscore-numpy](https://github.com/lichman0405/scscore-numpy.git)
- 可指定连接体数量（如 Zn4O + 3x linker）
- 支持 YAML 格式的 few-shot 示例管理
- 提供命令行界面（CLI）快速调用

---

## 项目结构

```
linkergen/
├── main.py                  # 命令行主入口
├── examples/run_example.py  # 示例脚本
├── linkergen/core/          # 各模块核心功能代码
├── linkergen/config/        # 示例连接体、API 密钥配置
├── linkergen/templates/     # Prompt 模板
├── requirements.txt         # 依赖清单
└── README.md
```

---

## 快速使用

### 安装依赖

```bash
pip install -r requirements.txt
```

### 运行示例

```bash
python main.py --metal Zn4O --n 6 --linker_count 2 
```

参数说明：

- `--metal`：目标金属簇，如 Zn4O、Cu2Paddlewheel
- `--n`：希望生成的连接体数量
- `--linker_count`：每个配方中连接体的数目（默认 3）
- `--examples`：few-shot 示例 YAML 路径
- `--config_path`：OpenAI API 密钥配置文件路径
- `--model_path`：SCScore 模型权重路径

---

## 配置说明

### `config/secrets.yaml`

```yaml
openai:
  api_key: "sk-..."
  base_url: "https://api.openai.com/v1"
  model: "gpt-4"
```

### `config/default_linkers.yaml`

包含用于 prompt 中引导模型生成的 SMILES 示例。

---

## 输出示例

```text
Generated Linkers:

- [Formula] Zn4O + 3x OC(=O)c1ccc(cc1)C(=O)O
  SMILES: OC(=O)c1ccc(cc1)C(=O)O
  SCScore: 2.15
```

---

## 后续计划

- 支持连接体结构自动分类（芳香族、杂环等）
- 集成 Allchemy 或 ASKCOS 评估合成路线
- 自动提取 linker 语料库（从 MOFDB / CoRE MOF）

---

## 协议

本项目基于 MIT 协议发布，欢迎自由使用和修改。