"""
Main orchestrator for the cancer research multi-agent system.
"""
import google.generativeai as genai
from typing import Dict, Any, List
try:
    from .config import (
        GOOGLE_API_KEY, MODEL_NAME, MAX_ITERATIONS, 
        GOAL_BINDING_ENERGY, PROBLEM_STATEMENT
    )
    from .agents import (
        SubAgent1Researcher, CollaborativeDiscussion, SubAgent4SolutionWriter
    )
except ImportError:
    from config import (
        GOOGLE_API_KEY, MODEL_NAME, MAX_ITERATIONS, 
        GOAL_BINDING_ENERGY, PROBLEM_STATEMENT
    )
    from agents import (
        SubAgent1Researcher, CollaborativeDiscussion, SubAgent4SolutionWriter
    )


class CancerResearchOrchestrator:
    """
    Main orchestrator that manages the multi-agent workflow for cancer research.
    """
    
    def __init__(self):
        # Configure Gemini API
        genai.configure(api_key=GOOGLE_API_KEY)
        
        # Initialize agents
        self.researcher = SubAgent1Researcher()
        self.discussion_manager = CollaborativeDiscussion()
        self.solution_writer = SubAgent4SolutionWriter()
        
        # Initialize context
        self.context = f"""
🎯 **PROBLEM STATEMENT:** {PROBLEM_STATEMENT}

🏆 **ULTIMATE GOAL:** Achieve a binding energy of {GOAL_BINDING_ENERGY} kcal/mol or lower (more negative is better)

📋 **TARGET SPECIFICATIONS:**
- Target: EGFR protein with L858R mutation
- Application: Non-small cell lung cancer treatment
- Metric: Binding energy (kcal/mol)
- Success Threshold: ≤ {GOAL_BINDING_ENERGY} kcal/mol
"""
        
        self.iteration_results = []
    
    def run_research_pipeline(self) -> Dict[str, Any]:
        """
        Runs the complete cancer research pipeline.
        
        Returns:
            Dictionary with final results and pipeline statistics
        """
        print("🚀" + "="*60)
        print("🧬 CANCER RESEARCH MULTI-AGENT PIPELINE STARTING")
        print("🎯 Goal: Design optimal EGFR L858R inhibitor")
        print("="*62)
        
        # Phase 1: Initial Research (runs once)
        print("\n📚 PHASE 1: INITIAL RESEARCH")
        research_findings = self.researcher.research(PROBLEM_STATEMENT)
        self.context += f"\n\n📊 **INITIAL RESEARCH FINDINGS:**\n{research_findings}"
        
        # Phase 2: Iterative Solution Development
        print(f"\n🔄 PHASE 2: ITERATIVE SOLUTION DEVELOPMENT (Max {MAX_ITERATIONS} iterations)")
        
        best_result = None
        
        for iteration in range(MAX_ITERATIONS):
            print(f"\n{'🔬' + '='*20} ITERATION {iteration + 1}/{MAX_ITERATIONS} {'='*20}")
            
            # Phase 2a: Collaborative Discussion
            final_molecule, discussion_history = self.discussion_manager.conduct_discussion(self.context)
            
            # Validate molecule
            if not final_molecule or len(final_molecule) < 5:
                print("❌ Discussion failed to produce valid molecule. Retrying...")
                feedback = "\n\n⚠️ **FEEDBACK:** Previous discussion failed to produce a valid SMILES string. Please propose a concrete, specific molecular structure."
                self.context += feedback
                continue
            
            # Phase 2b: AutoDock Vina Simulation
            simulation_result = self.solution_writer.run_simulation(final_molecule)
            binding_energy = simulation_result['binding_energy']
            
            # Store iteration result
            iteration_result = {
                'iteration': iteration + 1,
                'molecule_smiles': final_molecule,
                'binding_energy': binding_energy,
                'simulation_result': simulation_result,
                'discussion_history': discussion_history
            }
            self.iteration_results.append(iteration_result)
            
            # Update best result
            if best_result is None or binding_energy < best_result['binding_energy']:
                best_result = iteration_result.copy()
            
            # Phase 2c: Evaluation and Loop Control
            if binding_energy <= GOAL_BINDING_ENERGY:
                print(f"\n🎉{'*'*20} GOAL ACHIEVED! {'*'*20}")
                print(f"✅ Iteration: {iteration + 1}")
                print(f"🧬 Final Solution (SMILES): {final_molecule}")
                print(f"📊 Final Binding Energy: {binding_energy:.2f} kcal/mol")
                print(f"🎯 Target was: {GOAL_BINDING_ENERGY:.2f} kcal/mol")
                print(f"🏆 Improvement: {abs(binding_energy - GOAL_BINDING_ENERGY):.2f} kcal/mol better than target")
                
                return self._generate_final_report(iteration_result, success=True)
            
            else:
                print(f"\n📈 Goal not yet achieved in iteration {iteration + 1}")
                print(f"📊 Current: {binding_energy:.2f} kcal/mol | Target: {GOAL_BINDING_ENERGY:.2f} kcal/mol")
                print(f"📉 Gap: {binding_energy - GOAL_BINDING_ENERGY:.2f} kcal/mol")
                
                # Generate feedback for next iteration
                feedback = self._generate_iteration_feedback(iteration_result)
                self.context += feedback
        
        # Max iterations reached
        print(f"\n⏰ Maximum iterations ({MAX_ITERATIONS}) reached")
        print(f"🏆 Best result achieved:")
        print(f"🧬 Best Molecule: {best_result['molecule_smiles']}")
        print(f"📊 Best Binding Energy: {best_result['binding_energy']:.2f} kcal/mol")
        print(f"🎯 Target was: {GOAL_BINDING_ENERGY:.2f} kcal/mol")
        
        return self._generate_final_report(best_result, success=False)
    
    def _generate_iteration_feedback(self, iteration_result: Dict[str, Any]) -> str:
        """
        Generates feedback for the next iteration based on current results.
        
        Args:
            iteration_result: Results from current iteration
            
        Returns:
            Formatted feedback string
        """
        binding_energy = iteration_result['binding_energy']
        molecule = iteration_result['molecule_smiles']
        iteration_num = iteration_result['iteration']
        
        gap = binding_energy - GOAL_BINDING_ENERGY
        
        feedback = f"""

🔄 **ITERATION {iteration_num} RESULTS & FEEDBACK:**
- 🧬 **Tested Molecule:** {molecule}
- 📊 **Binding Energy:** {binding_energy:.2f} kcal/mol
- 🎯 **Target Energy:** {GOAL_BINDING_ENERGY:.2f} kcal/mol
- 📉 **Gap to Target:** {gap:.2f} kcal/mol (needs improvement)

💡 **FEEDBACK FOR NEXT ITERATION:**
The proposed structure did not achieve the target binding affinity. The binding energy needs to be more negative by {gap:.2f} kcal/mol. 

🔍 **Analysis Suggestions:**
1. Consider why this structure may have insufficient binding affinity
2. Explore modifications to increase hydrophobic interactions
3. Optimize hydrogen bonding patterns
4. Consider conformational flexibility and entropy effects
5. Address potential steric clashes or suboptimal geometry

🎯 **Next Steps:** Design a significantly improved molecule that addresses these binding limitations.
"""
        
        return feedback
    
    def _generate_final_report(self, final_result: Dict[str, Any], success: bool) -> Dict[str, Any]:
        """
        Generates comprehensive final report.
        
        Args:
            final_result: Final iteration result
            success: Whether the goal was achieved
            
        Returns:
            Comprehensive report dictionary
        """
        report = {
            'success': success,
            'goal_achieved': success,
            'target_binding_energy': GOAL_BINDING_ENERGY,
            'final_result': final_result,
            'total_iterations': len(self.iteration_results),
            'all_iterations': self.iteration_results,
            'best_binding_energy': min(r['binding_energy'] for r in self.iteration_results),
            'improvement_over_iterations': self._calculate_improvement_trend(),
            'problem_statement': PROBLEM_STATEMENT,
            'methodology': {
                'research_phase': 'ArXiv + Google Search',
                'discussion_agents': 'Medicinal Chemist + Computational Biologist',
                'simulation_tool': 'GROMACS via Inductiva API',
                'optimization_metric': 'Binding Energy (kcal/mol)'
            }
        }
        
        return report
    
    def _calculate_improvement_trend(self) -> List[float]:
        """Calculate the improvement trend across iterations."""
        return [result['binding_energy'] for result in self.iteration_results]


def main():
    """Main entry point for the cancer research pipeline."""
    try:
        orchestrator = CancerResearchOrchestrator()
        final_report = orchestrator.run_research_pipeline()
        
        print("\n📋" + "="*60)
        print("📊 FINAL PIPELINE REPORT")
        print("="*62)
        print(f"🎯 Goal Achieved: {'✅ YES' if final_report['success'] else '❌ NO'}")
        print(f"🧬 Best Molecule: {final_report['final_result']['molecule_smiles']}")
        print(f"📊 Best Binding Energy: {final_report['best_binding_energy']:.2f} kcal/mol")
        print(f"🔄 Total Iterations: {final_report['total_iterations']}")
        print("="*62)
        
        return final_report
        
    except Exception as e:
        print(f"❌ Pipeline failed with error: {e}")
        raise


if __name__ == "__main__":
    main()