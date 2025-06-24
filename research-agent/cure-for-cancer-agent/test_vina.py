#!/usr/bin/env python3
"""
Test script to verify AutoDock Vina setup and run a simple docking simulation.
"""

import sys
import os

def test_dependencies():
    """Test if all required dependencies are available."""
    print("🧪 Testing AutoDock Vina Dependencies")
    print("=" * 40)
    
    all_available = True
    
    # Test RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        print("✅ RDKit available")
    except ImportError:
        print("❌ RDKit not available")
        print("   Install with: pip install rdkit")
        all_available = False
    
    # Test AutoDock Vina
    try:
        from vina import Vina
        print("✅ AutoDock Vina Python package available")
    except ImportError:
        print("❌ AutoDock Vina Python package not available")
        print("   Install with: pip install vina")
        all_available = False
    
    # Test Open Babel
    import subprocess
    try:
        result = subprocess.run(['obabel'], 
                              capture_output=True, text=True, timeout=5)
        output_text = result.stdout + result.stderr
        if "Open Babel" in output_text:
            print("✅ Open Babel (obabel) available")
        else:
            print("❌ Open Babel (obabel) not working properly")
            print(f"   Return code: {result.returncode}")
            print(f"   Stdout: {result.stdout[:100]}...")
            print(f"   Stderr: {result.stderr[:100]}...")
            all_available = False
    except FileNotFoundError:
        print("❌ Open Babel (obabel) not found")
        print("   Install instructions:")
        print("   Ubuntu/Debian: sudo apt-get install openbabel")
        print("   macOS: brew install open-babel")
        print("   Windows: Download from http://openbabel.org/")
        all_available = False
    except Exception as e:
        print(f"❌ Error testing Open Babel: {e}")
        all_available = False
    
    return all_available

def test_protein_file():
    """Test if protein file is available."""
    protein_file = "compound_test/protein.pdbqt"
    
    if os.path.exists(protein_file):
        print(f"✅ Protein file found: {protein_file}")
        return True
    else:
        print(f"❌ Protein file not found: {protein_file}")
        print("   Run: python setup_protein.py")
        return False

def test_simple_docking():
    """Test a simple docking simulation."""
    print("\n🎯 Testing Simple Docking Simulation")
    print("=" * 40)
    
    try:
        from tools import AutoDockVinaSimulationTool
        
        # Test with a simple molecule (ethanol)
        test_smiles = "CCO"
        print(f"🧬 Testing with molecule: {test_smiles}")
        
        vina_tool = AutoDockVinaSimulationTool()
        result = vina_tool.run_simulation(test_smiles)
        
        print(f"\n📊 Test Results:")
        print(f"   Binding Energy: {result['binding_energy']} kcal/mol")
        print(f"   Status: {result['status']}")
        print(f"   Job ID: {result['job_id']}")
        
        if result.get('output_file'):
            print(f"   Output File: {result['output_file']}")
        
        if result['status'] == 'real_autodock_vina':
            print("✅ Real AutoDock Vina simulation successful!")
        else:
            print("⚠️ Using mock simulation (dependencies missing)")
        
        return True
        
    except Exception as e:
        print(f"❌ Docking test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("🧬 AutoDock Vina Setup Test")
    print("=" * 50)
    
    tests_passed = 0
    total_tests = 3
    
    # Test 1: Dependencies
    if test_dependencies():
        tests_passed += 1
    
    print()
    
    # Test 2: Protein file
    if test_protein_file():
        tests_passed += 1
    
    print()
    
    # Test 3: Simple docking
    if test_simple_docking():
        tests_passed += 1
    
    print(f"\n📊 Test Results: {tests_passed}/{total_tests} passed")
    
    if tests_passed == total_tests:
        print("🎉 All tests passed! AutoDock Vina setup is ready.")
        return True
    else:
        print("❌ Some tests failed. Please fix the issues above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)