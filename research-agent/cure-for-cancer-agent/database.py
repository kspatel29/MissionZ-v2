"""
Supabase database integration for the cancer research pipeline.
"""
import os
import json
import uuid
from datetime import datetime
from typing import Dict, Any, Optional
from supabase import create_client, Client

class DatabaseManager:
    """Manages database operations for the cancer research pipeline."""
    
    def __init__(self):
        # Get Supabase credentials from environment
        self.supabase_url = os.getenv("SUPABASE_URL")
        self.supabase_key = os.getenv("SUPABASE_ANON_KEY")
        
        if not self.supabase_url or not self.supabase_key:
            raise ValueError("SUPABASE_URL and SUPABASE_ANON_KEY must be set in .env file")
        
        # Initialize Supabase client
        self.supabase: Client = create_client(self.supabase_url, self.supabase_key)
        
        # Fixed space ID as provided
        self.space_id = "c544efa7-3f84-4fc4-9150-34380e25a521"
    
    def log_agent_conversation(self, agent_name: str, agent_response: str, 
                             iteration: int = None) -> bool:
        """
        Log an agent's conversation turn to the database.
        
        Args:
            agent_name: Name of the agent (e.g., "SubAgent1Researcher")
            agent_response: The agent's response text
            iteration: Optional iteration number
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Prepare the data
            conversation_data = {
                "research_space_id": self.space_id,
                "agent_name": agent_name,
                "agent_response": {
                    "response": agent_response,
                    "iteration": iteration,
                    "timestamp": datetime.now().isoformat()
                }
            }
            
            # Insert into database
            result = self.supabase.table("research_conversation_turns").insert(conversation_data).execute()
            
            if result.data:
                return True
            else:
                print(f"❌ Failed to log conversation for {agent_name}")
                return False
                
        except Exception as e:
            print(f"❌ Database error logging conversation: {e}")
            return False
    
    def log_solution(self, molecule_smiles: str, binding_energy: float, 
                    simulation_result: Dict[str, Any], html_file_path: str = None) -> bool:
        """
        Log a solution and its results to the database.
        
        Args:
            molecule_smiles: SMILES string of the proposed molecule
            binding_energy: Binding energy result
            simulation_result: Full simulation result dictionary
            html_file_path: Path to the HTML visualization file
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Read HTML file content directly for database storage
            html_content = None
            if html_file_path:
                # Construct the full path to the HTML file
                if not os.path.isabs(html_file_path):
                    # If it's a relative path, try to find it in the compound_test directory
                    full_path = os.path.join(os.path.dirname(__file__), "compound_test", html_file_path)
                else:
                    full_path = html_file_path
                
                if os.path.exists(full_path):
                    try:
                        # Read HTML content directly
                        with open(full_path, 'r', encoding='utf-8') as f:
                            html_content = f.read()
                        print(f"✅ HTML content loaded: {len(html_content)} characters")
                    except Exception as e:
                        print(f"⚠️ Could not read HTML file: {e}")
                else:
                    print(f"⚠️ HTML file not found: {full_path}")
                    print(f"   Original path: {html_file_path}")
                    print(f"   Current directory: {os.getcwd()}")
                    print(f"   Database file location: {os.path.dirname(__file__)}")
            
            # Prepare solution data with HTML content stored directly
            solution_data = {
                "research_space_id": self.space_id,
                "solution": molecule_smiles,
                "result": json.dumps({
                    "binding_energy": binding_energy,
                    "simulation_result": simulation_result,
                    "timestamp": datetime.now().isoformat(),
                    "molecule_smiles": molecule_smiles,
                    "status": simulation_result.get('status', 'unknown'),
                    "job_id": simulation_result.get('job_id', 'unknown')
                }),
                "simulation_file": html_content  # Store HTML content directly
            }
            
            # Insert into database
            result = self.supabase.table("solutions").insert(solution_data).execute()
            
            if result.data:
                solution_id = result.data[0]['id']
                print(f"✅ Solution logged to database (ID: {solution_id})")
                if html_content:
                    print(f"✅ HTML visualization stored in database ({len(html_content)} chars)")
                return True
            else:
                print(f"❌ Failed to log solution")
                return False
                
        except Exception as e:
            print(f"❌ Database error logging solution: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    
    def get_research_space_info(self) -> Optional[Dict[str, Any]]:
        """Get information about the current research space."""
        try:
            result = self.supabase.table("research_space").select("*").eq("id", self.space_id).execute()
            
            if result.data:
                return result.data[0]
            else:
                return None
                
        except Exception as e:
            print(f"❌ Database error getting space info: {e}")
            return None
    
    def get_conversation_history(self, limit: int = 50) -> list:
        """Get recent conversation history for the research space."""
        try:
            result = self.supabase.table("research_conversation_turns")\
                .select("*")\
                .eq("research_space_id", self.space_id)\
                .order("created_at", desc=True)\
                .limit(limit)\
                .execute()
            
            return result.data if result.data else []
            
        except Exception as e:
            print(f"❌ Database error getting conversation history: {e}")
            return []
    
    def get_solutions_history(self, limit: int = 20) -> list:
        """Get recent solutions for the research space."""
        try:
            result = self.supabase.table("solutions")\
                .select("*")\
                .eq("research_space_id", self.space_id)\
                .order("created_at", desc=True)\
                .limit(limit)\
                .execute()
            
            return result.data if result.data else []
            
        except Exception as e:
            print(f"❌ Database error getting solutions history: {e}")
            return []
    
    def get_existing_solutions_summary(self) -> str:
        """
        Get a formatted summary of existing solutions to prevent duplication.
        
        Returns:
            Formatted string with existing molecules and their results
        """
        try:
            # Get all existing solutions for this research space
            result = self.supabase.table("solutions")\
                .select("solution, result, created_at")\
                .eq("research_space_id", self.space_id)\
                .order("created_at", desc=True)\
                .execute()
            
            if not result.data:
                return "No previous solutions found in database."
            
            # Format the existing solutions
            solutions_summary = "🗃️ **EXISTING SOLUTIONS (DO NOT REPEAT):**\n\n"
            
            for i, solution in enumerate(result.data, 1):
                try:
                    # Parse the result JSON
                    result_data = json.loads(solution['result'])
                    binding_energy = result_data.get('binding_energy', 'Unknown')
                    status = result_data.get('status', 'Unknown')
                    timestamp = solution['created_at'][:10]  # Just the date
                    
                    solutions_summary += f"{i}. **SMILES:** {solution['solution']}\n"
                    solutions_summary += f"   **Binding Energy:** {binding_energy} kcal/mol\n"
                    solutions_summary += f"   **Status:** {status}\n"
                    solutions_summary += f"   **Date:** {timestamp}\n\n"
                    
                except Exception as e:
                    # If JSON parsing fails, just show basic info
                    solutions_summary += f"{i}. **SMILES:** {solution['solution']}\n"
                    solutions_summary += f"   **Date:** {solution['created_at'][:10]}\n\n"
            
            solutions_summary += "⚠️ **IMPORTANT:** Do NOT propose any of the above molecules again. "
            solutions_summary += "Focus on novel, different chemical structures.\n"
            
            print(f"📋 Retrieved {len(result.data)} existing solutions from database")
            return solutions_summary
            
        except Exception as e:
            print(f"❌ Database error getting existing solutions: {e}")
            return "Could not retrieve existing solutions from database."
    
    def check_molecule_exists(self, molecule_smiles: str) -> bool:
        """
        Check if a specific molecule has already been tested.
        
        Args:
            molecule_smiles: SMILES string to check
            
        Returns:
            True if molecule exists, False otherwise
        """
        try:
            result = self.supabase.table("solutions")\
                .select("id")\
                .eq("research_space_id", self.space_id)\
                .eq("solution", molecule_smiles)\
                .execute()
            
            return len(result.data) > 0 if result.data else False
            
        except Exception as e:
            print(f"❌ Database error checking molecule existence: {e}")
            return False