#!/usr/bin/env python3
"""
Test imports for the cancer research system.
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_imports():
    """Test all imports."""
    try:
        print("🔍 Testing config import...")
        from config import MODEL_NAME, GOAL_BINDING_ENERGY, INDUCTIVA_API_KEY
        print(f"✅ Config: Model={MODEL_NAME}, Goal={GOAL_BINDING_ENERGY}")
        
        print("\n🔍 Testing tools import...")
        from tools import ArxivSearchTool, GromacsSimulationTool
        print("✅ Tools imported successfully")
        
        print("\n🔍 Testing mock simulation...")
        tool = GromacsSimulationTool()
        result = tool.run_simulation('CCO')
        print(f"✅ Mock simulation: {result['binding_energy']} kcal/mol")
        
        print("\n🔍 Testing ArXiv tool...")
        arxiv_tool = ArxivSearchTool()
        print("✅ ArXiv tool initialized")
        
        print("\n✅ All imports successful!")
        return True
        
    except Exception as e:
        print(f"❌ Import failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_imports()
    sys.exit(0 if success else 1)