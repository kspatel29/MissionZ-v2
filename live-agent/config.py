"""
Configuration settings for the news agent system.
"""
import os
from dotenv import load_dotenv

load_dotenv()

# API Keys
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

if not GOOGLE_API_KEY:
    raise ValueError("FATAL: GOOGLE_API_KEY not found in .env file")

if not SUPABASE_URL or not SUPABASE_KEY:
    print("WARNING: SUPABASE_URL or SUPABASE_KEY not found in .env file. Database logging disabled.")

# Model Configuration
MODEL_NAME = os.getenv("MODEL_NAME", "gemini-2.5-flash")

# News Agent Configuration
DISCUSSION_ROUNDS = int(os.getenv("DISCUSSION_ROUNDS", "3"))
CONFIDENCE_THRESHOLD = float(os.getenv("CONFIDENCE_THRESHOLD", "0.8"))
NEWS_LOOKBACK_HOURS = int(os.getenv("NEWS_LOOKBACK_HOURS", "48"))

# Application Configuration
DEBUG = os.getenv("DEBUG", "True").lower() == "true"
HOST = os.getenv("HOST", "0.0.0.0")
PORT = int(os.getenv("PORT", "8000"))

# Agent Names
AGENT_NAMES = {
    "RESEARCHER": "sub_agent_1_researcher",
    "BIAS_ANALYZER": "sub_agent_2_bias_analyzer", 
    "COUNTER_CLAIMS": "sub_agent_3_counter_claims",
    "FORMATTER": "sub_agent_4_formatter"
}

print(f"✅ Configuration loaded - Model: {MODEL_NAME}")
print(f"✅ Discussion rounds: {DISCUSSION_ROUNDS}, Confidence threshold: {CONFIDENCE_THRESHOLD}")
print(f"✅ News lookback: {NEWS_LOOKBACK_HOURS} hours")
