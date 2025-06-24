#!/usr/bin/env python3
"""
Script to find the center of the original inhibitor (AQ4) from the 2ITV structure.
This determines the binding site coordinates for AutoDock Vina.
"""

def find_binding_center():
    """Find the center of the AQ4 inhibitor in the original 2ITV structure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolTransforms
        
        print("🔍 Finding binding site center from original inhibitor (AQ4)...")
        
        # Load the original PDB file with the inhibitor
        pdb_file = 'compound_test/2itv.pdb'
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        
        if mol is None:
            print("❌ Could not load PDB file")
            return None
        
        # Find the inhibitor (its residue name is 'AQ4')
        inhibitor_atoms = []
        for atom in mol.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            if pdb_info and pdb_info.GetResidueName() == 'AQ4':
                inhibitor_atoms.append(atom.GetIdx())
        
        if not inhibitor_atoms:
            print("⚠️ No AQ4 inhibitor found in PDB file")
            print("📍 Using default EGFR L858R binding site coordinates")
            return [22.597, -0.341, 27.054]
        
        print(f"✅ Found {len(inhibitor_atoms)} atoms in AQ4 inhibitor")
        
        # Get the center of the inhibitor
        conformer = mol.GetConformer()
        center = rdMolTransforms.ComputeCentroid(conformer, atomIndices=inhibitor_atoms)
        
        center_coords = [round(center.x, 3), round(center.y, 3), round(center.z, 3)]
        
        print("📍 Vina Search Box Center:")
        print(f"   center_x = {center_coords[0]:.3f}")
        print(f"   center_y = {center_coords[1]:.3f}")
        print(f"   center_z = {center_coords[2]:.3f}")
        
        print("\n📦 Recommended box size (in Angstroms):")
        print("   size_x = 25.0")
        print("   size_y = 25.0")
        print("   size_z = 25.0")
        
        return center_coords
        
    except ImportError:
        print("⚠️ RDKit not available for binding site detection")
        print("📍 Using default EGFR L858R binding site coordinates")
        return [22.597, -0.341, 27.054]
    except Exception as e:
        print(f"❌ Error finding binding site center: {e}")
        print("📍 Using default EGFR L858R binding site coordinates")
        return [22.597, -0.341, 27.054]

if __name__ == "__main__":
    center = find_binding_center()
    print(f"\n🎯 Final coordinates: {center}")