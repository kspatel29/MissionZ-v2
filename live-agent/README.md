# News Agent System

A sophisticated multi-agent architecture for news research, bias analysis, and verification. The system uses 4 specialized agents working in collaboration to process news topics, verify claims, discuss bias, and format structured output.

## 🏗️ Architecture Overview

The system consists of four specialized agents working in collaboration:

### 🔬 Sub-Agent 1: News Researcher
- **Role**: Research new news updates from the past 48 hours
- **Tools**: Google Search, URL Context, Current Date/Time
- **Output**: Single news update that doesn't exist in current timeline

### 🔍 Sub-Agent 2: Bias Analyzer  
- **Role**: Analyze news reports for bias, credibility, and source verification
- **Tools**: Google Search, URL Context, Current Date/Time
- **Output**: Bias analysis with credibility assessment

### ⚖️ Sub-Agent 3: Counter-Claims Agent
- **Role**: Provide counter-claims, verify sources, manage discussion confidence
- **Tools**: Google Search, URL Context, Current Date/Time  
- **Output**: Counter-claims with confidence scoring

### 📝 Sub-Agent 4: Report Formatter
- **Role**: Create structured JSON output from discussion results
- **Tools**: Current Date/Time
- **Output**: Formatted JSON with event_time, news_report, bias, credibility_score, sources

## 🔄 Workflow

1. **Research Phase**: Sub-Agent 1 finds new news updates not in timeline
2. **Discussion Phase**: Sub-Agents 2 & 3 discuss bias and verify claims until confidence threshold reached
3. **Formatting Phase**: Sub-Agent 4 creates structured JSON output
4. **Storage Phase**: Final report saved to news_timeline table

## 📊 Database Schema

The system uses three Supabase tables:

- **news_space**: Stores news space information
- **news_timeline**: Stores final news reports with structured data
- **news_conversation_turns**: Logs all agent interactions

## 🚀 Installation

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Configure environment variables in `.env`:
```bash
GOOGLE_API_KEY=your_google_api_key
SUPABASE_URL=your_supabase_url
SUPABASE_KEY=your_supabase_key
MODEL_NAME=gemini-2.5-flash
DISCUSSION_ROUNDS=3
CONFIDENCE_THRESHOLD=0.8
NEWS_LOOKBACK_HOURS=48
```

## 📖 Usage

### Single Processing Cycle
```bash
python main.py --topic "israel vs iran war" --single
```

### Continuous Monitoring
```bash
python main.py --topic "ukraine russia conflict" --continuous --interval 30
```

### Space Summary
```bash
python main.py --topic "any topic" --space-id "your-space-id" --summary
```

### Test Workflow
```bash
python main.py test
```

## 🔧 Configuration

Key configuration options in `config.py`:

- `DISCUSSION_ROUNDS`: Number of discussion rounds between agents 2 & 3
- `CONFIDENCE_THRESHOLD`: Confidence score needed to exit discussion (0.0-1.0)
- `NEWS_LOOKBACK_HOURS`: Hours to look back for existing timeline entries

## 📄 Output Format

The system outputs structured JSON with:

```json
{
  "event_time": "2025-07-07T10:30:00Z",
  "news_report": "Verified news summary",
  "bias": "Description of identified bias",
  "bias_percentage": 25,
  "credibility_score": 85,
  "sources": ["url1", "url2"],
  "verification_status": "verified",
  "confidence_level": "high",
  "claims_verified": ["claim1", "claim2"],
  "claims_disputed": ["disputed_claim"],
  "additional_context": "Important context"
}
```

## 🛠️ Development

### Project Structure
```
live-agent/
├── __init__.py
├── config.py          # Configuration settings
├── database.py        # Supabase database integration
├── tools.py           # Shared tools for all agents
├── agents.py          # All agent implementations
├── orchestrator.py    # Main workflow orchestrator
├── main.py           # Entry point
├── requirements.txt   # Dependencies
└── README.md         # This file
```

### Adding New Tools

1. Add tool function to `tools.py`
2. Add to `GEMINI_TOOLS` list
3. Update agent prompts to use new tool

### Customizing Agents

Each agent can be customized by modifying their prompts and logic in `agents.py`.

## 🔍 Monitoring

All agent interactions are logged to the `news_conversation_turns` table for monitoring and debugging.

## ⚠️ Notes

- Requires valid Google API key for Gemini and search
- Requires Supabase database setup
- Confidence scoring uses pattern matching from agent responses
- System designed for continuous operation
