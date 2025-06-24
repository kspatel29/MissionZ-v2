"""
Tool implementations for the cancer research agent system.
"""
import requests
import xml.etree.ElementTree as ET
from typing import Dict, List, Any
try:
    from .config import ARXIV_MAX_RESULTS
except ImportError:
    from config import ARXIV_MAX_RESULTS


class ArxivSearchTool:
    """Tool for searching ArXiv research papers."""
    
    def __init__(self):
        self.base_url = 'http://export.arxiv.org/api/query?'
    
    def search(self, query: str) -> str:
        """
        Searches arXiv for research papers and returns key findings.
        
        Args:
            query: Search query string
            
        Returns:
            Formatted string with paper titles, links, and summaries
        """
        print(f"🔍 Executing ArXiv Search Tool for: '{query}'")
        
        # Construct search query with proper encoding using urllib
        import urllib.parse
        import urllib.request
        
        # Clean and encode the query
        clean_query = query.replace('"', '').strip()
        encoded_query = urllib.parse.quote(clean_query)
        search_url = f'{self.base_url}search_query=all:{encoded_query}&start=0&max_results={ARXIV_MAX_RESULTS}'
        
        try:
            print(f"🌐 Fetching from: {search_url}")
            
            # Use urllib instead of requests
            with urllib.request.urlopen(search_url, timeout=30) as response:
                content = response.read().decode('utf-8')
            
            root = ET.fromstring(content)
            
            entries = []
            for entry in root.findall('{http://www.w3.org/2005/Atom}entry'):
                title_elem = entry.find('{http://www.w3.org/2005/Atom}title')
                summary_elem = entry.find('{http://www.w3.org/2005/Atom}summary')
                link_elem = entry.find('{http://www.w3.org/2005/Atom}id')
                
                if title_elem is not None and summary_elem is not None and link_elem is not None:
                    title = title_elem.text.strip().replace('\n', ' ')
                    summary = summary_elem.text.strip().replace('\n', ' ')
                    # Truncate summary to 200 chars for conciseness
                    summary = summary[:200] + "..." if len(summary) > 200 else summary
                    link = link_elem.text.strip()
                    
                    entries.append(f"📄 {title}\n🔗 {link}\n📝 {summary}\n{'-'*30}")
            
            result = "\n".join(entries) if entries else "No relevant results found on arXiv for that query."
            print(f"✅ Found {len(entries)} ArXiv papers")
            return result
            
        except Exception as e:
            error_msg = f"❌ Error occurred during arXiv search: {e}"
            print(error_msg)
            return error_msg


class WebSearchTool:
    """Tool for web search using Google Search API."""
    
    def search(self, query: str) -> str:
        """
        Performs web search and returns summarized findings.
        Note: This will be handled by Gemini's built-in Google Search tool.
        """
        # This is a placeholder - the actual implementation uses Gemini's built-in tool
        return f"Web search results for: {query}"


class AutoDockVinaSimulationTool:
    """Tool for running AutoDock Vina molecular docking simulations."""
    
    def __init__(self):
        self.protein_file = "research-agent/cure-for-cancer-agent/compound_test/protein.pdbqt"
        self.search_center = [22.597, -0.341, 27.054]  # EGFR L858R binding site center
        self.search_box_size = [25.0, 25.0, 25.0]  # Search box size in Angstroms
        self.working_dir = "research-agent/cure-for-cancer-agent/compound_test"
        
        # Check if required dependencies are available
        self.use_real_vina = self._check_dependencies()
        self.use_obabel = False  # Will be set in _check_dependencies
        
        if not self.use_real_vina:
            print("⚠️ Using mock AutoDock Vina simulation (dependencies not available)")
    
    def _check_dependencies(self) -> bool:
        """Check if AutoDock Vina and required dependencies are available."""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            print("✅ RDKit available")
            
            # Check if vina python package is available
            try:
                from vina import Vina
                print("✅ AutoDock Vina Python package available")
                
                # Try to check for Open Babel (optional)
                try:
                    import subprocess
                    result = subprocess.run(['obabel', '--version'], 
                                          capture_output=True, text=True, timeout=5)
                    if result.returncode == 0:
                        print("✅ Open Babel (obabel) available")
                        self.use_obabel = True
                    else:
                        print("⚠️ Open Babel not available - using RDKit for format conversion")
                        self.use_obabel = False
                except:
                    print("⚠️ Open Babel not available - using RDKit for format conversion")
                    self.use_obabel = False
                
                return True
            except ImportError:
                print("⚠️ AutoDock Vina Python package not found")
                return False
                
        except ImportError:
            print("⚠️ RDKit not available")
            return False
        except Exception as e:
            print(f"⚠️ Error checking dependencies: {e}")
            return False
    
    def _create_pdbqt_from_rdkit(self, mol, filename: str) -> bool:
        """Create a PDBQT file from RDKit molecule (simplified version)."""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Add hydrogens if not present
            mol_with_h = Chem.AddHs(mol)
            
            # Generate conformer if not present
            if mol_with_h.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol_with_h)
            
            # Write to PDB first
            pdb_filename = filename.replace('.pdbqt', '.pdb')
            Chem.MolToPDBFile(mol_with_h, pdb_filename)
            
            # Create a simplified PDBQT file (without proper charges/atom types)
            # This is a basic conversion - real PDBQT requires proper preparation
            with open(pdb_filename, 'r') as pdb_file, open(filename, 'w') as pdbqt_file:
                for line in pdb_file:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Convert PDB line to basic PDBQT format
                        # This is simplified - real PDBQT needs proper atom types and charges
                        pdbqt_line = line.rstrip() + "    0.00  0.00    A\n"
                        pdbqt_file.write(pdbqt_line)
                    elif line.startswith('END'):
                        pdbqt_file.write(line)
            
            print(f"✅ Created simplified PDBQT file: {filename}")
            return True
            
        except Exception as e:
            print(f"❌ Error creating PDBQT file: {e}")
            return False
    
    def run_simulation(self, molecule_smiles: str, target_protein: str = "EGFR_L858R") -> Dict[str, Any]:
        """
        Runs an AutoDock Vina simulation to calculate binding energy.
        
        Args:
            molecule_smiles: SMILES string of the molecule
            target_protein: Target protein identifier
            
        Returns:
            Dictionary with simulation results including binding energy
        """
        print(f"🧬 Running AutoDock Vina simulation for molecule: {molecule_smiles}")
        
        if self.use_real_vina:
            return self._run_real_vina_simulation(molecule_smiles, target_protein)
        else:
            return self._run_mock_simulation(molecule_smiles)
    
    def _run_real_vina_simulation(self, molecule_smiles: str, target_protein: str) -> Dict[str, Any]:
        """Run actual AutoDock Vina simulation."""
        try:
            import subprocess
            import os
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from vina import Vina
            
            print("🔧 Preparing ligand from SMILES...")
            
            # Change to working directory
            original_dir = os.getcwd()
            os.chdir(self.working_dir)
            
            try:
                # 1. Prepare the Ligand (from SMILES to 3D PDBQT)
                mol = Chem.MolFromSmiles(molecule_smiles)
                if mol is None:
                    raise ValueError(f"Invalid SMILES string: {molecule_smiles}")
                
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates
                embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                if embed_result != 0:
                    print("⚠️ Warning: 3D embedding may not be optimal")
                
                # Optimize geometry
                AllChem.MMFFOptimizeMolecule(mol)
                
                # Save to PDB file
                Chem.MolToPDBFile(mol, 'ligand.pdb')
                print("✅ 3D ligand structure generated")
                
                # Convert PDB to PDBQT format
                print("🔄 Converting to PDBQT format...")
                if self.use_obabel:
                    # Use Open Babel if available
                    subprocess.run(
                        'obabel ligand.pdb -O ligand.pdbqt -p 7.4',
                        shell=True, check=True, capture_output=True
                    )
                    print("✅ Ligand prepared with Open Babel: ligand.pdbqt")
                else:
                    # Use RDKit-based conversion
                    if self._create_pdbqt_from_rdkit(mol, 'ligand.pdbqt'):
                        print("✅ Ligand prepared with RDKit: ligand.pdbqt")
                    else:
                        raise Exception("Failed to create PDBQT file")
                
                # 2. Configure and Run Vina
                print("🔬 Starting AutoDock Vina docking...")
                v = Vina(sf_name='vina')
                
                # Check if protein file exists
                if not os.path.exists('protein.pdbqt'):
                    raise FileNotFoundError("protein.pdbqt not found in compound_test directory")
                
                v.set_receptor('protein.pdbqt')
                v.set_ligand_from_file('ligand.pdbqt')
                
                # Set search space
                v.compute_vina_maps(center=self.search_center, box_size=self.search_box_size)
                
                # 3. Perform Docking
                print("🎯 Performing molecular docking...")
                v.dock(exhaustiveness=8, n_poses=1)
                
                # 4. Get Results
                energies = v.energies(n_poses=1)
                best_energy = energies[0][0]  # Get the energy of the best pose
                
                # Save docked pose
                output_file = f'ligand_docked_{hash(molecule_smiles) % 10000}.pdbqt'
                v.write_poses(output_file, n_poses=1, overwrite=True)
                
                print("✅ AutoDock Vina simulation completed!")
                
                return {
                    "binding_energy": round(best_energy, 2),
                    "simulation_time": "real_vina",
                    "status": "real_autodock_vina",
                    "job_id": f"vina_{hash(molecule_smiles) % 10000}",
                    "output_file": os.path.join(self.working_dir, output_file),
                    "details": {
                        "note": "Real AutoDock Vina simulation completed",
                        "molecule_smiles": molecule_smiles,
                        "target_protein": target_protein,
                        "search_center": self.search_center,
                        "search_box_size": self.search_box_size,
                        "exhaustiveness": 8
                    }
                }
                
            finally:
                # Return to original directory
                os.chdir(original_dir)
                
        except Exception as e:
            print(f"❌ Error in real AutoDock Vina simulation: {e}")
            print("🔄 Falling back to mock simulation...")
            return self._run_mock_simulation(molecule_smiles)
    
    def _run_mock_simulation(self, molecule_smiles: str) -> Dict[str, Any]:
        """Run mock simulation for testing purposes."""
        import hashlib
        
        # Create a deterministic but varied binding energy based on SMILES
        hash_value = int(hashlib.md5(molecule_smiles.encode()).hexdigest()[:8], 16)
        
        # Base binding energy around -7 kcal/mol (typical for AutoDock Vina)
        base_energy = -7.0
        
        # Add complexity bonus (longer, more complex molecules might bind better)
        complexity_bonus = min(len(molecule_smiles) / 80.0, 3.0)
        
        # Add realistic variation
        variation = (hash_value % 400) / 100.0  # 0 to 4 kcal/mol variation
        
        # Calculate final binding energy
        binding_energy = base_energy - complexity_bonus - variation
        
        return {
            "binding_energy": round(binding_energy, 2),
            "simulation_time": "mock_vina",
            "status": "mock_autodock_vina",
            "job_id": f"mock_vina_{hash_value}",
            "output_file": None,
            "details": {
                "note": "Mock AutoDock Vina simulation for testing purposes",
                "molecule_smiles": molecule_smiles,
                "complexity_score": complexity_bonus,
                "variation": variation
            }
        }


# Tool function wrappers for Gemini function calling
def arxiv_search_tool(query: str) -> str:
    """Function wrapper for ArXiv search tool."""
    tool = ArxivSearchTool()
    return tool.search(query)


def autodock_vina_simulation_tool(molecule_smiles: str) -> str:
    """Function wrapper for AutoDock Vina simulation tool."""
    tool = AutoDockVinaSimulationTool()
    result = tool.run_simulation(molecule_smiles)
    
    return f"""
🧪 **AutoDock Vina Simulation Results**
📊 **Binding Energy:** {result['binding_energy']} kcal/mol
⏱️ **Simulation Time:** {result['simulation_time']}
✅ **Status:** {result['status']}
🆔 **Job ID:** {result['job_id']}
📁 **Output File:** {result.get('output_file', 'N/A')}
"""