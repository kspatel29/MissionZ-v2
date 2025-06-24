#!/usr/bin/env python3
"""
Test script to verify AutoDock Vina docking with known EGFR inhibitors.
"""

import os
from tools import AutoDockVinaSimulationTool

def test_known_inhibitors():
    """Test docking with known EGFR inhibitors."""
    
    # Known EGFR inhibitors with their SMILES
    known_inhibitors = {
        "Gefitinib": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "Erlotinib": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",
        "Osimertinib": "COc1cc(N(C)CCN(C)C)ccc1Nc1nccc(-c2cn(C)c3ccccc23)n1",
        "Lapatinib": "CS(=O)(=O)CCNCc1oc(-c2ccc(F)cc2)cc1-c1cccc(Cl)c1Nc1ncnc2cc(OCC3CCOCC3)ccc12"
    }
    
    print("🧪 Testing AutoDock Vina with Known EGFR Inhibitors")
    print("=" * 60)
    
    vina_tool = AutoDockVinaSimulationTool()
    
    for name, smiles in known_inhibitors.items():
        print(f"\n🧬 Testing {name}...")
        print(f"   SMILES: {smiles}")
        
        try:
            result = vina_tool.run_simulation(smiles)
            
            print(f"   📊 Binding Energy: {result['binding_energy']} kcal/mol")
            print(f"   ✅ Status: {result['status']}")
            
            if result['binding_energy'] != 0.0:
                print(f"   🎉 Success! Non-zero binding energy achieved")
            else:
                print(f"   ⚠️ Warning: Zero binding energy")
                
        except Exception as e:
            print(f"   ❌ Error: {e}")
    
    print(f"\n📋 Test completed!")

if __name__ == "__main__":
    test_known_inhibitors()