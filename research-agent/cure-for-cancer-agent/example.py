#!/usr/bin/env python3
"""
Example usage of the Cancer Research Multi-Agent System.
"""

import os
import sys
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def run_example():
    """Run a simple example of the cancer research pipeline."""
    
    print("🧬 Cancer Research Multi-Agent System - Example")
    print("=" * 50)
    
    # Check if API key is available
    if not os.getenv("GOOGLE_API_KEY"):
        print("❌ GOOGLE_API_KEY not found in .env file")
        print("Please create a .env file with your Google API key:")
        print("GOOGLE_API_KEY=your_api_key_here")
        return False
    
    try:
        # Import the main orchestrator
        from main import CancerResearchOrchestrator
        
        print("🚀 Initializing Cancer Research Orchestrator...")
        orchestrator = CancerResearchOrchestrator()
        
        print("📋 Current Configuration:")
        print(f"  - Target: EGFR L858R mutation")
        print(f"  - Goal: ≤ -10.0 kcal/mol binding energy")
        print(f"  - Max Iterations: 10")
        print(f"  - Model: Gemini 2.5 Flash")
        
        print("\n🔬 Starting research pipeline...")
        print("This may take several minutes depending on API response times...")
        
        # Run the pipeline
        results = orchestrator.run_research_pipeline()
        
        # Display results
        print("\n📊 FINAL RESULTS:")
        print("=" * 30)
        print(f"🎯 Goal Achieved: {'✅ YES' if results['success'] else '❌ NO'}")
        print(f"🧬 Best Molecule: {results['final_result']['molecule_smiles']}")
        print(f"📊 Best Binding Energy: {results['best_binding_energy']:.2f} kcal/mol")
        print(f"🎯 Target Energy: {results['target_binding_energy']:.2f} kcal/mol")
        print(f"🔄 Total Iterations: {results['total_iterations']}")
        
        if results['success']:
            improvement = abs(results['best_binding_energy'] - results['target_binding_energy'])
            print(f"🏆 Exceeded target by: {improvement:.2f} kcal/mol")
        
        print("\n📈 Iteration Progress:")
        for i, iteration in enumerate(results['all_iterations']):
            status = "✅" if iteration['binding_energy'] <= results['target_binding_energy'] else "📊"
            print(f"  {status} Iteration {i+1}: {iteration['binding_energy']:.2f} kcal/mol")
        
        return True
        
    except KeyboardInterrupt:
        print("\n⏹️ Example interrupted by user")
        return False
    except Exception as e:
        print(f"\n❌ Example failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def run_quick_test():
    """Run a quick test of individual components."""
    
    print("🧪 Quick Component Test")
    print("=" * 25)
    
    try:
        # Test tool imports
        from tools import ArxivSearchTool, GromacsSimulationTool
        print("✅ Tools imported successfully")
        
        # Test mock simulation
        gromacs_tool = GromacsSimulationTool()
        test_molecule = "CCO"  # Simple ethanol
        result = gromacs_tool.run_simulation(test_molecule)
        print(f"✅ Mock simulation test: {result['binding_energy']} kcal/mol")
        
        # Test ArXiv tool (without API call)
        arxiv_tool = ArxivSearchTool()
        print("✅ ArXiv tool initialized")
        
        print("\n🎉 Quick test completed successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Quick test failed: {e}")
        return False

if __name__ == "__main__":
    print("Choose an option:")
    print("1. Run full example (requires API keys)")
    print("2. Run quick component test")
    print("3. Exit")
    
    choice = input("\nEnter choice (1-3): ").strip()
    
    if choice == "1":
        success = run_example()
    elif choice == "2":
        success = run_quick_test()
    elif choice == "3":
        print("👋 Goodbye!")
        sys.exit(0)
    else:
        print("❌ Invalid choice")
        sys.exit(1)
    
    sys.exit(0 if success else 1)