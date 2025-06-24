#!/usr/bin/env python3
"""
Test script to verify the setup and configuration.
"""

import os
import sys
from dotenv import load_dotenv

def test_environment():
    """Test environment variables and API keys."""
    print("🔍 Testing environment setup...")
    
    load_dotenv()
    
    google_key = os.getenv("GOOGLE_API_KEY")
    inductiva_key = os.getenv("INDUCTIVA_API_KEY")
    
    if google_key:
        print("✅ GOOGLE_API_KEY found")
    else:
        print("❌ GOOGLE_API_KEY not found in .env file")
        return False
    
    if inductiva_key:
        print("✅ INDUCTIVA_API_KEY found")
    else:
        print("⚠️ INDUCTIVA_API_KEY not found (will use mock simulation)")
    
    return True

def test_imports():
    """Test all required imports."""
    print("\n🔍 Testing imports...")
    
    try:
        import google.generativeai as genai
        print("✅ google.generativeai imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import google.generativeai: {e}")
        return False
    
    try:
        import requests
        print("✅ requests imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import requests: {e}")
        return False
    
    try:
        import inductiva
        print("✅ inductiva imported successfully")
    except ImportError as e:
        print(f"⚠️ Failed to import inductiva: {e} (will use mock simulation)")
    
    return True

def test_tools():
    """Test tool functionality."""
    print("\n🔍 Testing tools...")
    
    try:
        from tools import ArxivSearchTool, GromacsSimulationTool
        
        # Test ArXiv tool
        arxiv_tool = ArxivSearchTool()
        print("✅ ArxivSearchTool initialized")
        
        # Test GROMACS tool
        gromacs_tool = GromacsSimulationTool()
        print("✅ GromacsSimulationTool initialized")
        
        # Test mock simulation
        result = gromacs_tool.run_simulation("CCO")  # Simple ethanol molecule
        print(f"✅ Mock simulation test: {result['binding_energy']} kcal/mol")
        
        return True
        
    except Exception as e:
        print(f"❌ Tool test failed: {e}")
        return False

def test_agents():
    """Test agent initialization."""
    print("\n🔍 Testing agents...")
    
    try:
        from config import GOOGLE_API_KEY
        import google.generativeai as genai
        
        if not GOOGLE_API_KEY:
            print("❌ Cannot test agents without GOOGLE_API_KEY")
            return False
        
        genai.configure(api_key=GOOGLE_API_KEY)
        
        from agents import SubAgent1Researcher, CollaborativeDiscussion, SubAgent4SolutionWriter
        
        # Test basic initialization (don't run actual API calls)
        print("✅ SubAgent1Researcher class imported")
        print("✅ CollaborativeDiscussion class imported") 
        print("✅ SubAgent4SolutionWriter class imported")
        
        return True
        
    except Exception as e:
        print(f"❌ Agent test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🧪 Cancer Research Multi-Agent System - Setup Test")
    print("=" * 50)
    
    tests = [
        test_environment,
        test_imports,
        test_tools,
        test_agents
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
        print("🎉 All tests passed! System is ready to run.")
        return True
    else:
        print("❌ Some tests failed. Please check the configuration.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)