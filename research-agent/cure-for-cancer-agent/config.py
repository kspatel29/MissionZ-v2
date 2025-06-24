"""
Configuration settings for the cancer research agent system.
"""
import os
from dotenv import load_dotenv

load_dotenv()

# API Keys
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_ANON_KEY = os.getenv("SUPABASE_ANON_KEY")

if not GOOGLE_API_KEY:
    raise ValueError("FATAL: GOOGLE_API_KEY not found in .env file")

if not SUPABASE_URL or not SUPABASE_ANON_KEY:
    print("WARNING: SUPABASE_URL or SUPABASE_ANON_KEY not found in .env file. Database logging disabled.")

# Model Configuration
MODEL_NAME = "gemini-2.5-flash"
MAX_ITERATIONS = 100
GOAL_BINDING_ENERGY = -16.0  # kcal/mol (more negative is better)

# Research Configuration
ARXIV_MAX_RESULTS = 5
DISCUSSION_ROUNDS = 7

# EGFR Target Information
TARGET_PROTEIN = "EGFR"
TARGET_MUTATION = "L858R"
PROBLEM_STATEMENT = f"Design a small molecule inhibitor for the {TARGET_MUTATION} mutated {TARGET_PROTEIN} protein, a key factor in non-small cell lung cancer, to achieve better binding affinity than existing drugs."