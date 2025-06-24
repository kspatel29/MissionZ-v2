#!/usr/bin/env python3
"""
Script to download and prepare the EGFR L858R protein structure for AutoDock Vina.
"""

import os
import requests
import subprocess
from pathlib import Path

def download_protein_structure():
    """Download the EGFR L858R protein structure from PDB."""
    
    compound_dir = Path("compound_test")
    compound_dir.mkdir(exist_ok=True)
    
    pdb_file = compound_dir / "2itv.pdb"
    protein_file = compound_dir / "protein.pdb"
    protein_pdbqt = compound_dir / "protein.pdbqt"
    
    # Download PDB file if it doesn't exist
    if not pdb_file.exists():
        print("📥 Downloading EGFR L858R structure (PDB ID: 2ITV)...")
        url = "https://files.rcsb.org/download/2ITV.pdb"
        
        response = requests.get(url)
        response.raise_for_status()
        
        with open(pdb_file, 'w') as f:
            f.write(response.text)
        print(f"✅ Downloaded: {pdb_file}")
    else:
        print(f"✅ PDB file already exists: {pdb_file}")
    
    # Clean PDB file (remove HETATM lines for ligands and water)
    if not protein_file.exists():
        print("🧹 Cleaning PDB file (removing ligands and water)...")
        
        with open(pdb_file, 'r') as infile, open(protein_file, 'w') as outfile:
            for line in infile:
                # Keep only ATOM lines (protein atoms), skip HETATM (ligands, water)
                if line.startswith('ATOM'):
                    outfile.write(line)
                elif line.startswith('END'):
                    outfile.write(line)
        
        print(f"✅ Cleaned protein file: {protein_file}")
    else:
        print(f"✅ Clean protein file already exists: {protein_file}")
    
    # Convert to PDBQT format using Open Babel
    if not protein_pdbqt.exists():
        print("🔄 Converting to PDBQT format...")
        
        try:
            # Check if obabel is available
            result = subprocess.run(['obabel', '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode != 0:
                print("❌ Open Babel (obabel) not found. Please install it:")
                print("   Ubuntu/Debian: sudo apt-get install openbabel")
                print("   macOS: brew install open-babel")
                print("   Windows: Download from http://openbabel.org/")
                return False
            
            # Convert PDB to PDBQT
            cmd = f'obabel {protein_file} -O {protein_pdbqt} -xr'
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            
            print(f"✅ PDBQT protein file created: {protein_pdbqt}")
            
        except subprocess.CalledProcessError as e:
            print(f"❌ Error converting to PDBQT: {e}")
            print(f"Command output: {e.stdout}")
            print(f"Command error: {e.stderr}")
            return False
        except FileNotFoundError:
            print("❌ Open Babel (obabel) not found in PATH")
            return False
    else:
        print(f"✅ PDBQT protein file already exists: {protein_pdbqt}")
    
    return True

def find_binding_site_center():
    """Find the center of the original inhibitor for binding site definition."""
    
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolTransforms
        
        pdb_file = Path("compound_test/2itv.pdb")
        
        if not pdb_file.exists():
            print("❌ Original PDB file not found")
            return None
        
        print("🔍 Finding binding site center from original inhibitor (AQ4)...")
        
        # Load the original PDB file with the inhibitor
        mol = Chem.MolFromPDBFile(str(pdb_file), removeHs=False)
        
        if mol is None:
            print("❌ Could not load PDB file with RDKit")
            return None
        
        # Find the inhibitor (its residue name is 'AQ4')
        inhibitor_atoms = []
        for atom in mol.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            if pdb_info and pdb_info.GetResidueName() == 'AQ4':
                inhibitor_atoms.append(atom.GetIdx())
        
        if not inhibitor_atoms:
            print("⚠️ No AQ4 inhibitor found in PDB file")
            # Use default coordinates for EGFR L858R
            center = [22.597, -0.341, 27.054]
            print(f"📍 Using default binding site center: {center}")
            return center
        
        # Get the center of the inhibitor
        conformer = mol.GetConformer()
        center = rdMolTransforms.ComputeCentroid(conformer, atomIndices=inhibitor_atoms)
        
        center_coords = [round(center.x, 3), round(center.y, 3), round(center.z, 3)]
        
        print(f"📍 Binding site center found: {center_coords}")
        print(f"📦 Recommended search box size: [25.0, 25.0, 25.0]")
        
        return center_coords
        
    except ImportError:
        print("⚠️ RDKit not available for binding site detection")
        print("📍 Using default binding site center: [22.597, -0.341, 27.054]")
        return [22.597, -0.341, 27.054]
    except Exception as e:
        print(f"❌ Error finding binding site center: {e}")
        print("📍 Using default binding site center: [22.597, -0.341, 27.054]")
        return [22.597, -0.341, 27.054]

def main():
    """Main setup function."""
    print("🧬 EGFR L858R Protein Setup for AutoDock Vina")
    print("=" * 50)
    
    # Download and prepare protein
    if download_protein_structure():
        print("\n✅ Protein preparation completed successfully!")
        
        # Find binding site center
        center = find_binding_site_center()
        
        print(f"\n📋 AutoDock Vina Configuration:")
        print(f"   Protein file: compound_test/protein.pdbqt")
        print(f"   Search center: {center}")
        print(f"   Search box size: [25.0, 25.0, 25.0]")
        
        print(f"\n🎯 Ready for molecular docking simulations!")
        
    else:
        print("\n❌ Protein preparation failed!")
        return False
    
    return True

if __name__ == "__main__":
    main()