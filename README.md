# CIDD: Collaborative Intelligence for Structure-Based Drug Design

[![NeurIPS 2025](https://img.shields.io/badge/NeurIPS-2025-blue.svg)](https://openreview.net/pdf?id=7k7cubl1iL)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

[**Paper**](https://openreview.net/pdf?id=7k7cubl1iL)

> **⚠️ UNDER CONSTRUCTION:** This repository is currently still under construction. Updates will be provided soon.

> **💡 Important - Flexible LLM Integration:** The core of CIDD is the `send_request()` function in `src/generation/cidd_generation.py`. **You can customize this function to use ANY LLM API of your choice** (OpenAI, Azure, Claude, local models, etc.) as long as it can call a language model and return responses. The provided implementation supports Azure OpenAI and DeepSeek (via Volcengine Ark) by default, but feel free to adapt it to your preferred LLM service.

---

## Installation

### Prerequisites
- Python 3.8+
- Conda (recommended for RDKit and OpenBabel)
- LLM API key (Azure OpenAI or DeepSeek via Volcengine Ark)

### Setup

```bash
# Clone repository
git clone https://github.com/bowen-gao/CIDD.git
cd CIDD

# Create environment
conda create -n cidd python=3.8
conda activate cidd

# Install dependencies
conda install -c conda-forge rdkit openbabel
pip install -r requirements.txt
pip install -e .

# Configure LLM API (see note above for customization)
# Option 1: Azure OpenAI (Default)
export AZURE_OPENAI_API_KEY="your_azure_api_key"
export AZURE_OPENAI_URL="https://your-resource.openai.azure.com/openai/deployments/your-deployment/chat/completions?api-version=2024-08-01-preview"

# Option 2: DeepSeek via Volcengine Ark
export VOLCENGINE_ARK_API_KEY="your_volcengine_ark_api_key"
export VOLCENGINE_MODEL="deepseek-v3-241226"  # or deepseek-r1-250120, or endpoint ID like ep-xxxxx
```

**LLM Configuration Options:**
- **Azure OpenAI** (default): `AZURE_OPENAI_API_KEY`, `AZURE_OPENAI_URL` (full endpoint URL)
- **DeepSeek via Volcengine Ark**: `VOLCENGINE_ARK_API_KEY`, `VOLCENGINE_MODEL` (model name or endpoint ID)
  - Model names: `deepseek-v3-241226`, `deepseek-r1-250120`
  - Or use endpoint ID format: `ep-xxxxx` (from Volcengine console)
- **Custom**: Modify `send_request()` in `src/generation/cidd_generation.py` to use your preferred LLM

---

## Running CrossDocked Benchmark

### Data Preparation

1. **CrossDocked2020 Test Set**: Download and set `CIDD_DATA_ROOT="/path/to/crossdocked_data/test_set"`

2. **Initial Molecules**: Download from [Google Drive](https://drive.google.com/drive/folders/1A3Mthm9ksbfUnMCe5T2noGsiEV1RfChH) (provided by Keyue Qiu)

### Run Benchmark

```bash
cd examples

# Set environment variables
export CIDD_DATA_ROOT="/path/to/crossdocked_data/test_set"
export AZURE_OPENAI_API_KEY="your_api_key"
export AZURE_OPENAI_URL="your_azure_endpoint_url"

# Run benchmark
python bench_crossdocked.py
```

### Configuration

Edit `bench_crossdocked.py` to adjust:
- Path to SBIG initial molecules (`.pt` file)
- Number of molecules to process
- Number of parallel workers

Key environment variables:
- `CIDD_DATA_ROOT`: Path to CrossDocked2020 test set (required)
- `CIDD_OUTPUT_DIR`: Output directory (default: `./results`)
- LLM API keys: See Installation section above

---

## Citation

```bibtex
@inproceedings{gao2025cidd,
  title={CIDD: Collaborative Intelligence for Structure-Based Drug Design Empowered by LLMs},
  author={Gao, Bowen and others},
  booktitle={Advances in Neural Information Processing Systems},
  year={2025},
  url={https://openreview.net/pdf?id=7k7cubl1iL}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file for details.
