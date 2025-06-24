# Cancer Research Multi-Agent System

A sophisticated multi-agent architecture for cancer research, specifically designed to discover novel small molecule inhibitors for the EGFR L858R mutation in non-small cell lung cancer.

## 🏗️ Architecture Overview

The system consists of four specialized agents working in collaboration:

### 🔬 Sub-Agent 1: Research Aggregator
- **Role**: Comprehensive research and information gathering
- **Tools**: 
  - ArXiv Research Agent (academic papers)
  - Google Search Agent (latest discoveries and solutions)
- **Output**: Curated list of current discoveries, solutions, and important findings

### 👨‍🔬 Sub-Agent 2: Medicinal Chemist
- **Role**: Molecular design and chemical innovation
- **Tools**: Research Agent, Web Search Agent
- **Capabilities**: Proposes novel molecular structures and chemical modifications

### 👩‍🔬 Sub-Agent 3: Computational Biologist  
- **Role**: Validation and critical analysis
- **Tools**: Research Agent, Web Search Agent
- **Capabilities**: Critiques proposals, validates feasibility, suggests improvements

### 🧪 Sub-Agent 4: Solution Writer & Simulator
- **Role**: AutoDock Vina molecular docking simulation
- **Tools**: AutoDock Vina for molecular docking
- **Target**: EGFR L858R mutated protein
- **Metric**: Binding energy optimization

## 🔄 Workflow Process

1. **Research Phase**: Sub-Agent 1 aggregates current knowledge
2. **Collaborative Discussion**: Sub-Agents 2 & 3 discuss and design solutions
3. **Molecular Docking**: Sub-Agent 4 tests against EGFR L858R using AutoDock Vina
4. **Iterative Improvement**: Loop continues until binding energy target is met

## 🎯 Success Criteria

- **Target**: Binding energy ≤ -10.0 kcal/mol
- **Protein**: EGFR with L858R mutation
- **Application**: Non-small cell lung cancer treatment

## 🚀 Quick Start

### Prerequisites

1. **Google API Key**: Required for Gemini 2.5 Flash model
2. **AutoDock Vina Dependencies**: RDKit, Vina Python package, Open Babel

### Installation

```bash
cd research-agent/cure-for-cancer-agent
pip install -r requirements.txt
```

### Configuration

1. Copy the example environment file:
```bash
cp .env.example .env
```

2. Edit `.env` with your API keys:
```bash
GOOGLE_API_KEY=your_google_api_key_here
INDUCTIVA_API_KEY=your_inductiva_api_key_here  # Optional
```

### Running the System

```python
from main import main

# Run the complete pipeline
results = main()
```

Or run directly:
```bash
python -m main
```

## 📊 Expected Output

The system will output:
- Research findings from ArXiv and web search
- Collaborative discussion between agents
- Molecular proposals (SMILES format)
- GROMACS simulation results
- Binding energy measurements
- Iterative improvements until goal is achieved

## 🛠️ Technical Implementation

### API Integration

- **Gemini 2.5 Flash**: Function calling for tool integration
- **Google Search**: Built-in grounding for web search
- **ArXiv API**: Academic paper retrieval and analysis
- **Inductiva API**: GROMACS molecular dynamics simulations

### Key Features

- **Separation of Concerns**: Different models for different tool types to avoid API conflicts
- **Iterative Feedback Loop**: Continuous improvement based on simulation results
- **Collaborative Intelligence**: Multi-agent discussion for creative solutions
- **Real Simulation**: Actual GROMACS testing via Inductiva cloud platform

### Error Handling

- Graceful fallback to mock simulations if Inductiva API unavailable
- Robust error handling for API failures
- Validation of molecular structures (SMILES format)

## 📁 Project Structure

```
research-agent/cure-for-cancer-agent/
├── __init__.py
├── config.py          # Configuration and constants
├── tools.py           # Tool implementations (ArXiv, GROMACS, etc.)
├── agents.py          # Agent implementations
├── main.py            # Main orchestrator
├── requirements.txt   # Python dependencies
├── .env.example       # Environment variables template
└── README.md          # This file
```

## 🔧 Customization

### Modifying Target Protein
Edit `config.py`:
```python
TARGET_PROTEIN = "EGFR"
TARGET_MUTATION = "L858R"
```

### Adjusting Success Criteria
Edit `config.py`:
```python
GOAL_BINDING_ENERGY = -10.0  # kcal/mol
```

### Changing Discussion Rounds
Edit `config.py`:
```python
DISCUSSION_ROUNDS = 3
```

## 📚 References

- [ArXiv API Documentation](https://info.arxiv.org/help/api/user-manual.html)
- [Google Search Grounding](https://ai.google.dev/gemini-api/docs/google-search)
- [Gemini Function Calling](https://ai.google.dev/gemini-api/docs/function-calling)
- [Inductiva GROMACS Guide](https://inductiva.ai/guides/gromacs/setup-test)
- [Inductiva Platform Overview](https://inductiva.ai/guides/documentation/intro/how_it_works)

## 🤝 Contributing

This system is designed to be modular and extensible. You can:
- Add new research tools
- Implement additional agent types
- Modify simulation parameters
- Extend the feedback loop logic

## ⚠️ Important Notes

- The system uses Gemini 2.5 Flash with function calling
- Google Search and custom tools are separated to avoid API conflicts
- Mock simulations are used when Inductiva API is unavailable
- The system is optimized for EGFR L858R but can be adapted for other targets