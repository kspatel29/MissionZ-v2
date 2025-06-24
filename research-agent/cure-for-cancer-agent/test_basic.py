#!/usr/bin/env python3
"""
Basic test script that doesn't require google-generativeai library.
"""

import sys
import os

def test_basic_imports():
    """Test basic imports without google-generativeai."""
    print("🧪 Testing basic imports...")
    
    try:
        import requests
        print("✅ requests imported successfully")
    except ImportError as e:
        print(f"❌ requests import failed: {e}")
        return False
    
    try:
        import xml.etree.ElementTree as ET
        print("✅ xml.etree.ElementTree imported successfully")
    except ImportError as e:
        print(f"❌ xml.etree.ElementTree import failed: {e}")
        return False
    
    try:
        from dotenv import load_dotenv
        print("✅ python-dotenv imported successfully")
    except ImportError as e:
        print(f"❌ python-dotenv import failed: {e}")
        return False
    
    return True

def test_config():
    """Test config module."""
    print("\n🔧 Testing config module...")
    
    try:
        from config import MODEL_NAME, GOAL_BINDING_ENERGY, TARGET_PROTEIN, TARGET_MUTATION
        print(f"✅ Config loaded successfully")
        print(f"  - Model: {MODEL_NAME}")
        print(f"  - Goal: {GOAL_BINDING_ENERGY} kcal/mol")
        print(f"  - Target: {TARGET_PROTEIN} {TARGET_MUTATION}")
        return True
    except Exception as e:
        print(f"❌ Config test failed: {e}")
        return False

def test_arxiv_tool():
    """Test ArXiv tool without google-generativeai."""
    print("\n📚 Testing ArXiv tool...")
    
    try:
        from tools import ArxivSearchTool
        
        tool = ArxivSearchTool()
        print("✅ ArxivSearchTool initialized")
        
        # Test with a simple query (this will make a real API call)
        print("🔍 Testing ArXiv search...")
        result = tool.search("EGFR inhibitor")
        print(f"✅ ArXiv search completed: {len(result)} characters returned")
        
        return True
    except Exception as e:
        print(f"❌ ArXiv tool test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_gromacs_tool():
    """Test GROMACS simulation tool."""
    print("\n🧪 Testing GROMACS simulation tool...")
    
    try:
        from tools import GromacsSimulationTool
        
        tool = GromacsSimulationTool()
        print("✅ GromacsSimulationTool initialized")
        
        # Test mock simulation
        test_molecules = ["CCO", "CC(C)C", "C1=CC=CC=C1"]
        
        for molecule in test_molecules:
            result = tool.run_simulation(molecule)
            print(f"✅ Mock simulation for {molecule}: {result['binding_energy']} kcal/mol")
        
        return True
    except Exception as e:
        print(f"❌ GROMACS tool test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all basic tests."""
    print("🧬 Cancer Research System - Basic Tests")
    print("=" * 45)
    
    tests = [
        test_basic_imports,
        test_config,
        test_arxiv_tool,
        test_gromacs_tool
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("📊 Test Results:")
    print(f"✅ Passed: {passed}/{total}")
    
    if passed == total:
        print("🎉 All basic tests passed!")
        print("📝 Note: Google Gemini API tests require proper GLIBC version")
        return True
    else:
        print("❌ Some tests failed.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)