#!/usr/bin/env python3
"""
3D Visualization tool for AutoDock Vina docking results.
Creates interactive HTML files for viewing protein-ligand complexes.
"""

import os
import json
from pathlib import Path
from datetime import datetime

class VisualizationTool:
    """Tool for creating 3D visualizations of docking results."""
    
    def __init__(self):
        # Check if we're already in compound_test directory
        current_dir = Path.cwd()
        if current_dir.name == "compound_test":
            # We're already in compound_test, so just create visualizations
            self.vis_dir = Path("visualizations")
        else:
            # We're in the parent directory, so create compound_test/visualizations
            self.vis_dir = Path("compound_test/visualizations")
        
        self.vis_dir.mkdir(parents=True, exist_ok=True)
    
    def create_3d_visualization(self, protein_file: str, ligand_file: str, 
                              molecule_smiles: str, binding_energy: float,
                              job_id: str) -> str:
        """
        Create an interactive 3D visualization HTML file.
        
        Args:
            protein_file: Path to protein PDBQT file
            ligand_file: Path to docked ligand PDBQT file
            molecule_smiles: SMILES string of the molecule
            binding_energy: Binding energy in kcal/mol
            job_id: Unique job identifier
            
        Returns:
            Path to the created HTML file
        """
        try:
            # Read protein and ligand files
            protein_content = self._read_pdbqt_file(protein_file)
            ligand_content = self._read_pdbqt_file(ligand_file)
            
            if not protein_content or not ligand_content:
                raise Exception("Could not read protein or ligand files")
            
            # Create HTML visualization
            html_file = self._create_html_visualization(
                protein_content, ligand_content, molecule_smiles, 
                binding_energy, job_id
            )
            
            print(f"✅ 3D visualization created: {html_file}")
            return str(html_file)
            
        except Exception as e:
            print(f"❌ Error creating 3D visualization: {e}")
            return None
    
    def _read_pdbqt_file(self, file_path: str) -> str:
        """Read and return the content of a PDBQT file."""
        try:
            with open(file_path, 'r') as f:
                return f.read()
        except Exception as e:
            print(f"❌ Error reading file {file_path}: {e}")
            return None
    
    def _create_html_visualization(self, protein_content: str, ligand_content: str,
                                 molecule_smiles: str, binding_energy: float,
                                 job_id: str) -> Path:
        """Create the HTML file with 3D visualization."""
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        html_filename = f"docking_result_{job_id}_{timestamp}.html"
        html_file = self.vis_dir / html_filename
        
        # Create HTML content with 3Dmol.js
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EGFR L858R Docking Result - {job_id}</title>
    <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            text-align: center;
        }}
        .info-panel {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            padding: 20px;
            background-color: #f8f9fa;
        }}
        .info-box {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .viewer-container {{
            height: 600px;
            position: relative;
            border-top: 2px solid #eee;
        }}
        .controls {{
            position: absolute;
            top: 10px;
            right: 10px;
            z-index: 1000;
            background: rgba(255,255,255,0.9);
            padding: 10px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.2);
        }}
        .control-btn {{
            display: block;
            width: 100%;
            margin: 5px 0;
            padding: 8px 12px;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 12px;
        }}
        .control-btn:hover {{
            background: #5a6fd8;
        }}
        .energy-display {{
            font-size: 24px;
            font-weight: bold;
            color: {'#28a745' if binding_energy <= -10.0 else '#dc3545'};
        }}
        .smiles-box {{
            font-family: monospace;
            background: #f1f3f4;
            padding: 10px;
            border-radius: 4px;
            word-break: break-all;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 EGFR L858R Molecular Docking Result</h1>
            <p>AutoDock Vina Simulation - Job ID: {job_id}</p>
        </div>
        
        <div class="info-panel">
            <div class="info-box">
                <h3>📊 Binding Results</h3>
                <p><strong>Binding Energy:</strong> <span class="energy-display">{binding_energy:.2f} kcal/mol</span></p>
                <p><strong>Target:</strong> EGFR L858R mutant protein</p>
                <p><strong>Goal:</strong> ≤ -10.0 kcal/mol</p>
                <p><strong>Status:</strong> {'✅ Goal Achieved!' if binding_energy <= -10.0 else '📈 Needs Improvement'}</p>
            </div>
            
            <div class="info-box">
                <h3>🧪 Molecule Information</h3>
                <p><strong>SMILES:</strong></p>
                <div class="smiles-box">{molecule_smiles}</div>
                <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
            </div>
        </div>
        
        <div class="viewer-container">
            <div class="controls">
                <button class="control-btn" onclick="resetView()">🔄 Reset View</button>
                <button class="control-btn" onclick="toggleProtein()">👁️ Toggle Protein</button>
                <button class="control-btn" onclick="toggleLigand()">💊 Toggle Ligand</button>
                <button class="control-btn" onclick="showBindingSite()">🎯 Binding Site</button>
                <button class="control-btn" onclick="toggleSurface()">🌊 Surface</button>
            </div>
            <div id="viewer" style="height: 100%; width: 100%;"></div>
        </div>
    </div>

    <script>
        let viewer;
        let proteinModel, ligandModel;
        let showProtein = true, showLigand = true, showSurface = false;
        
        // Initialize 3Dmol viewer
        function initViewer() {{
            viewer = $3Dmol.createViewer("viewer", {{
                defaultcolors: $3Dmol.rasmolElementColors
            }});
            
            // Load protein
            proteinModel = viewer.addModel(`{protein_content}`, "pdbqt");
            proteinModel.setStyle({{}}, {{
                cartoon: {{color: 'lightblue', opacity: 0.8}},
                stick: {{hidden: true}}
            }});
            
            // Load ligand
            ligandModel = viewer.addModel(`{ligand_content}`, "pdbqt");
            ligandModel.setStyle({{}}, {{
                stick: {{colorscheme: 'yellowCarbon', radius: 0.3}},
                sphere: {{scale: 0.3, colorscheme: 'yellowCarbon'}}
            }});
            
            // Set initial view to show both protein and ligand
            viewer.zoomTo();
            viewer.render();
            
            // Center on the binding site (both protein and ligand)
            setTimeout(() => {{
                // First try to center on ligand
                viewer.zoomTo(ligandModel);
                // Then zoom out a bit to show context
                viewer.zoom(0.8);
                viewer.render();
            }}, 500);
        }}
        
        function resetView() {{
            // Reset to show both protein and ligand optimally
            viewer.zoomTo(ligandModel);
            viewer.zoom(0.8);
            viewer.render();
        }}
        
        function toggleProtein() {{
            showProtein = !showProtein;
            if (showProtein) {{
                proteinModel.setStyle({{}}, {{
                    cartoon: {{color: 'lightblue', opacity: 0.8}}
                }});
            }} else {{
                proteinModel.setStyle({{}}, {{hidden: true}});
            }}
            viewer.render();
        }}
        
        function toggleLigand() {{
            showLigand = !showLigand;
            if (showLigand) {{
                ligandModel.setStyle({{}}, {{
                    stick: {{colorscheme: 'yellowCarbon', radius: 0.3}},
                    sphere: {{scale: 0.3, colorscheme: 'yellowCarbon'}}
                }});
            }} else {{
                ligandModel.setStyle({{}}, {{hidden: true}});
            }}
            viewer.render();
        }}
        
        function showBindingSite() {{
            // Highlight residues within 5Å of ligand
            proteinModel.setStyle({{}}, {{cartoon: {{color: 'lightblue', opacity: 0.8}}}});
            proteinModel.setStyle({{resi: ['858', '790', '793', '745']}}, {{
                stick: {{colorscheme: 'greenCarbon'}},
                cartoon: {{color: 'green'}}
            }});
            viewer.render();
        }}
        
        function toggleSurface() {{
            showSurface = !showSurface;
            if (showSurface) {{
                viewer.addSurface($3Dmol.SurfaceType.VDW, {{
                    opacity: 0.3,
                    color: 'white'
                }}, {{model: proteinModel}});
            }} else {{
                viewer.removeSurface();
            }}
            viewer.render();
        }}
        
        // Initialize when page loads
        document.addEventListener('DOMContentLoaded', initViewer);
    </script>
</body>
</html>
"""
        
        # Write HTML file
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        return html_file

def create_visualization(protein_file: str, ligand_file: str, 
                        molecule_smiles: str, binding_energy: float,
                        job_id: str) -> str:
    """
    Convenience function to create a 3D visualization.
    
    Args:
        protein_file: Path to protein PDBQT file
        ligand_file: Path to docked ligand PDBQT file
        molecule_smiles: SMILES string of the molecule
        binding_energy: Binding energy in kcal/mol
        job_id: Unique job identifier
        
    Returns:
        Path to the created HTML file
    """
    viz_tool = VisualizationTool()
    return viz_tool.create_3d_visualization(
        protein_file, ligand_file, molecule_smiles, binding_energy, job_id
    )

if __name__ == "__main__":
    # Test with existing files
    protein_file = "compound_test/protein.pdbqt"
    ligand_files = list(Path("compound_test").glob("ligand_docked_*.pdbqt"))
    
    if ligand_files:
        ligand_file = ligand_files[0]
        test_smiles = "CCO"  # Ethanol for testing
        test_energy = -5.2
        test_job_id = "test_7551"
        
        html_file = create_visualization(
            protein_file, str(ligand_file), test_smiles, test_energy, test_job_id
        )
        
        if html_file:
            print(f"🎉 Test visualization created: {html_file}")
            print(f"📂 Open in browser: file://{Path(html_file).absolute()}")
    else:
        print("❌ No docked ligand files found for testing")