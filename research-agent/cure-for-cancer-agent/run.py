#!/usr/bin/env python3
"""
Simple runner script for the Cancer Research Multi-Agent System.
"""

import sys
import os

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from main import main

if __name__ == "__main__":
    print("🚀 Starting Cancer Research Multi-Agent System...")
    try:
        results = main()
        print("\n✅ Pipeline completed successfully!")
        print(f"🎯 Goal achieved: {'Yes' if results['success'] else 'No'}")
        print(f"🧬 Best molecule: {results['final_result']['molecule_smiles']}")
        print(f"📊 Best binding energy: {results['best_binding_energy']:.2f} kcal/mol")
    except KeyboardInterrupt:
        print("\n⏹️ Pipeline interrupted by user")
    except Exception as e:
        print(f"\n❌ Pipeline failed: {e}")
        sys.exit(1)