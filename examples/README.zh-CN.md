# MoleCode 示例

[English](README.md) | **中文**

MoleCode 的可运行示例。在仓库根目录下运行其中任意一个:

```bash
pip install -e .
python examples/01_molecule_roundtrip.py
```

| 脚本 | 演示内容 | 需要 LLM? |
| --- | --- | --- |
| `01_molecule_roundtrip.py` | SMILES → MoleCode 图 → SMILES(无损) | 否 |
| `02_polymer_roundtrip.py` | PSMILES → MoleCode 图(`×n`)→ PSMILES;嵌段共聚物 | 否 |
| `03_markush_roundtrip.py` | 马库什 `{}` 缩写节点;缩写感知的图同构 | 否 |
| `04_understanding.py` | **理解** — 数原子 / 分子式 / 环 | 可选 |
| `05_generation.py` | **生成** — 约束下的从头设计 | 可选 |
| `06_editing.py` | **编辑** — 局部图编辑(增/删/替换) | 可选 |
| `07_reasoning.py` | **推理** — 反应产物预测 | 可选 |
| `08_image_to_molecode.py` | **OCSR** — 分子图片 → MoleCode 图 | 视觉模型 |

`08_image_to_molecode.py` 演示 OCSR(光学化学结构识别):它把一个示例分子渲染成 PNG
(或使用你传入的图片路径),将图片连同 Markush 版 MoleCode 系统提示一起发给**支持视觉的**
模型,再把返回的 MoleCode 图解析回结构。把 `MOLECODE_MODEL` 设为视觉模型(如 `gpt-4o-mini`、
`gpt-4o`);未设 API key 时为空跑(dry-run)。

## 四类 LLM 任务

示例 04–07 各自构建你会发给任意 LLM 的确切提示词,并通过 `examples/_llm.py`(对库内
`molecode.llm.LLMClient` 的薄封装)发送。它们**默认离线运行**——未设置 API key 时,
会打印拼装好的 system + user 提示词(一次"空跑"),让你确切看到 MoleCode 发送了什么。
要真正调用模型,设置环境变量即可(任何 OpenAI 兼容端点都行):

```bash
export MOLECODE_API_KEY="sk-..."
export MOLECODE_BASE_URL="https://api.openai.com/v1"   # 或你的服务商
export MOLECODE_MODEL="gpt-4o-mini"                    # 或任意对话模型
python examples/04_understanding.py
```

在你自己的代码里,可直接使用该客户端(api key 与 url 由你提供):

```python
from molecode import LLMClient
client = LLMClient(api_key="sk-...", base_url="https://api.openai.com/v1", model="gpt-4o-mini")
reply = client.chat(user_prompt, system=system_prompt)
```

MoleCode 自身从不调用 LLM——这个客户端只是一个可选的、零依赖的便利工具。可复用的核心组件:

* `molecode.prompts.MOLECULE_SYSTEM_PROMPT` / `MARKUSH_SYSTEM_PROMPT` —— 作为系统提示
  交给模型的文法规范;
* `molecode.molecule.mol_to_mermaid(...)` —— 把你的分子转成模型读取的图;
* `molecode.molecule.mermaid_to_mol(...)` —— 用 RDKit 确定性地解析并**校验**模型的输出。
