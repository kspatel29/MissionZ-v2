"""
Configuration settings for the cancer research agent system.
"""
import os
from dotenv import load_dotenv

load_dotenv()

# API Keys
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")

if not GOOGLE_API_KEY:
    raise ValueError("FATAL: GOOGLE_API_KEY not found in .env file")

# Model Configuration
MODEL_NAME = "gemini-2.5-flash"
MAX_ITERATIONS = 10
GOAL_BINDING_ENERGY = -10.0  # kcal/mol (more negative is better)

# Research Configuration
ARXIV_MAX_RESULTS = 5
DISCUSSION_ROUNDS = 7

# EGFR Target Information
TARGET_PROTEIN = "EGFR"
TARGET_MUTATION = "L858R"
PROBLEM_STATEMENT = f"Design a small molecule inhibitor for the {TARGET_MUTATION} mutated {TARGET_PROTEIN} protein, a key factor in non-small cell lung cancer, to achieve better binding affinity than existing drugs."