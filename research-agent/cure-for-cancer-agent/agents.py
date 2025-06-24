"""
Agent implementations for the cancer research system.
"""
import google.generativeai as genai
from google.generativeai.types import Tool
from typing import Dict, Any, List
import os
try:
    from .config import GOOGLE_API_KEY, MODEL_NAME, DISCUSSION_ROUNDS, TARGET_MUTATION, TARGET_PROTEIN
    from .tools import arxiv_search_tool, autodock_vina_simulation_tool, AutoDockVinaSimulationTool
    from .database import DatabaseManager
except ImportError:
    from config import GOOGLE_API_KEY, MODEL_NAME, DISCUSSION_ROUNDS, TARGET_MUTATION, TARGET_PROTEIN
    from tools import arxiv_search_tool, autodock_vina_simulation_tool, AutoDockVinaSimulationTool
    from database import DatabaseManager



class WebSearchAgent:
    """Specialized agent for Google Search using new API."""
    
    def __init__(self):
        try:
            # Try the new Google GenAI client approach
            from google import genai
            from google.genai import types
            
            self.client = genai.Client(api_key=GOOGLE_API_KEY)
            self.grounding_tool = types.Tool(google_search=types.GoogleSearch())
            self.config = types.GenerateContentConfig(tools=[self.grounding_tool])
            self.use_new_api = True
            print("✅ Using new Google GenAI client for search")
            
        except ImportError:
            # Fallback to old API
            self.model = genai.GenerativeModel(
                model_name=MODEL_NAME,
                tools=[Tool(google_search_retrieval={})]
            )
            self.use_new_api = False
            print("⚠️ Using fallback Google search API")
    
    def search(self, query: str) -> str:
        """Perform web search and return results."""
        try:
            if self.use_new_api:
                response = self.client.models.generate_content(
                    model="gemini-2.5-flash",
                    contents=query,
                    config=self.config
                )
                return response.text if response.text else "Web search completed but no results returned"
            else:
                response = self.model.generate_content(query)
                if response.candidates and response.candidates[0].content.parts:
                    return response.text
                else:
                    return "Web search failed: No valid response parts returned"
        except Exception as e:
            return f"Web search failed: {e}"


class ArxivSearchAgent:
    """Specialized agent for ArXiv search using direct API calls."""
    
    def __init__(self):
        # Use direct ArXiv tool instead of Gemini function calling
        from tools import ArxivSearchTool
        self.arxiv_tool = ArxivSearchTool()
    
    def search(self, query: str) -> str:
        """Perform ArXiv search and return results."""
        try:
            # Use direct ArXiv API call
            results = self.arxiv_tool.search(query)
            return results
        except Exception as e:
            return f"ArXiv search failed: {e}"


class SubAgent1Researcher:
    """
    Research Agent that coordinates between specialized search agents.
    Uses separate agents to avoid tool conflicts.
    """
    
    def __init__(self, db_manager: DatabaseManager = None):
        # Initialize specialized search agents
        self.web_search_agent = WebSearchAgent()
        self.arxiv_search_agent = ArxivSearchAgent()
        
        # Coordinator model (no tools, just coordination)
        self.coordinator_model = genai.GenerativeModel(MODEL_NAME)
        
        # Database manager for logging
        self.db_manager = db_manager
    
    def research(self, problem: str) -> str:
        """
        Conducts comprehensive research on the given problem.
        
        Args:
            problem: Problem statement to research
            
        Returns:
            Comprehensive research findings
        """
        print("🔬 SUB-AGENT 1: Research Phase")
        
        # Step 1: Web search for current discoveries and solutions
        web_search_prompt = f"""
        Search for and summarize the latest (2023-2025) research discoveries and solutions for:
        - Small molecule inhibitors targeting {TARGET_MUTATION} mutated {TARGET_PROTEIN}
        - Recent breakthroughs in non-small cell lung cancer treatment
        - Novel approaches to overcome drug resistance in {TARGET_PROTEIN} inhibitors
        - Current binding affinity improvements and molecular design strategies
        
        Focus on: current discoveries, exact solutions, important findings, and research articles.
        Provide specific examples of successful molecules and their binding energies if available.
        """
        
        web_findings = self.web_search_agent.search(web_search_prompt)
        if "failed" in web_findings.lower():
            web_findings = "Web search unavailable due to API limitations."
        
        # Step 2: ArXiv search for academic research
        arxiv_search_query = f"{TARGET_PROTEIN} inhibitor {TARGET_MUTATION} mutation small molecule"
        arxiv_query = f"Search ArXiv for research papers about: {arxiv_search_query}. Focus on novel chemical structures and binding affinity improvements."
        arxiv_findings = self.arxiv_search_agent.search(arxiv_query)
        if "failed" in arxiv_findings.lower():
            arxiv_findings = "ArXiv search unavailable."
        
        # Compile comprehensive research report
        research_report = f"""
🔬 **COMPREHENSIVE RESEARCH FINDINGS**

📡 **WEB SEARCH DISCOVERIES & SOLUTIONS:**
{web_findings}

📚 **ARXIV ACADEMIC RESEARCH:**
{arxiv_findings}

🎯 **RESEARCH SUMMARY:**
The above findings provide current state-of-the-art approaches, successful molecular designs, 
and important discoveries in {TARGET_PROTEIN} {TARGET_MUTATION} inhibitor development.
"""
        
        # Log to database
        if self.db_manager:
            self.db_manager.log_agent_conversation("SubAgent1Researcher", research_report)
        
        print("✅ Research completed")
        return research_report


class SubAgent2Chemist:
    """
    Medicinal Chemist Agent - proposes and refines molecular structures.
    Has access to research_agent and web_search_agent tools.
    """
    
    def __init__(self, db_manager: DatabaseManager = None):
        self.model = genai.GenerativeModel(MODEL_NAME)
        # Add specialized search agents for additional research if needed
        self.web_search_agent = WebSearchAgent()
        self.arxiv_search_agent = ArxivSearchAgent()
        self.db_manager = db_manager
    
    def propose_solution(self, context: str, discussion_history: List[str]) -> str:
        """
        Proposes a molecular solution based on context and discussion.
        
        Args:
            context: Current context including research and feedback
            discussion_history: Previous discussion rounds
            
        Returns:
            Proposed molecular solution with reasoning
        """
        full_context = context + "\n\n" + "\n".join(discussion_history)
        
        prompt = f"""
You are a medicinal chemist specializing in {TARGET_PROTEIN} inhibitors.

CONTEXT: {full_context[-1000:]}...

TASK: Propose a novel small molecule for {TARGET_MUTATION} mutated {TARGET_PROTEIN} with improved binding affinity.

REQUIREMENTS:
1. Target the {TARGET_MUTATION} mutation specifically
2. Aim for binding energy ≤ -10.0 kcal/mol
3. Address resistance mechanisms
4. Be synthetically feasible
5. MUST provide a valid SMILES string

RESPONSE FORMAT (BE CONCISE - MAX 300 WORDS):
- Molecule name and core scaffold
- Key design features (2-3 main points)
- Expected binding improvements
- **VALID SMILES string (REQUIRED)** - Use standard drug-like molecules as templates

IMPORTANT: 
1. Base your SMILES on known EGFR inhibitors like:
   - Gefitinib: COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1
   - Erlotinib: C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1
2. DO NOT repeat any molecules from the existing solutions list above
3. Propose NOVEL, DIFFERENT chemical structures only

Focus on practical, innovative design with VALID chemistry and NO DUPLICATES.
"""
        
        response = self.model.generate_content(prompt)
        if response.candidates and response.candidates[0].content.parts:
            agent_response = response.text
            # Log to database
            if self.db_manager:
                self.db_manager.log_agent_conversation("SubAgent2Chemist", agent_response)
            return agent_response
        else:
            return "I apologize, but I cannot generate a molecular proposal at this time. Please try rephrasing your request."


class SubAgent3Biologist:
    """
    Computational Biologist Agent - validates and critiques molecular proposals.
    Has access to research_agent and web_search_agent tools.
    """
    
    def __init__(self, db_manager: DatabaseManager = None):
        self.model = genai.GenerativeModel(MODEL_NAME)
        # Add specialized search agents for additional research if needed
        self.web_search_agent = WebSearchAgent()
        self.arxiv_search_agent = ArxivSearchAgent()
        self.db_manager = db_manager
    
    def critique_solution(self, context: str, discussion_history: List[str]) -> str:
        """
        Critiques and validates the proposed molecular solution.
        
        Args:
            context: Current context including research and feedback
            discussion_history: Previous discussion rounds
            
        Returns:
            Critical analysis and suggestions for improvement
        """
        full_context = context + "\n\n" + "\n".join(discussion_history)
        
        prompt = f"""
You are a computational biologist expert in {TARGET_PROTEIN} structure.

CONTEXT: {full_context[-1000:]}...

TASK: Critically analyze the medicinal chemist's molecular proposal.

ANALYSIS REQUIREMENTS:
1. Evaluate binding viability for {TARGET_PROTEIN} {TARGET_MUTATION}
2. Assess drug-like properties and synthetic feasibility
3. Identify weaknesses and suggest improvements
4. Consider resistance potential

RESPONSE FORMAT (BE CONCISE - MAX 250 WORDS):
- Overall assessment (viable/needs work)
- Check if molecule is duplicate from existing solutions
- 2-3 main strengths
- 2-3 main weaknesses  
- Specific improvement suggestions
- Binding energy prediction

Be constructively critical and focused. REJECT any duplicate molecules.

CRITICAL INSTRUCTION FOR CONCLUDING THE DISCUSSION:

Your default behavior is to continue the discussion for further refinement. However, if you become highly confident that the current molecule is a strong candidate and further discussion is unlikely to yield significant improvements, you MUST end your response with one of the exact phrase on a new line:

CONCLUSION: Proceed to final simulation.

If you do not include the exact phrase, the discussion will automatically continue. Do not use these phrases unless you are confident.
"""
        
        response = self.model.generate_content(prompt)
        if response.candidates and response.candidates[0].content.parts:
            agent_response = response.text
            # Log to database
            if self.db_manager:
                self.db_manager.log_agent_conversation("SubAgent3Biologist", agent_response)
            return agent_response
        else:
            return "I apologize, but I cannot generate a critique at this time. Please try rephrasing your request."



    def is_confident_about_solution(self, critique_response: str) -> bool:
        """
        Analyzes the critique response to determine if the agent is confident
        about the proposed solution and wants to end the discussion early.
        
        Args:
            critique_response: The critique response from this agent
            
        Returns:
            True if confident and wants to end discussion, False otherwise
        """
        # Convert to lowercase for easier matching
        response_lower = critique_response.lower()
        
        # Keywords indicating high confidence
        confidence_indicators = [
            "excellent proposal",
            "highly promising",
            "strong candidate",
            "ready for testing",
            "confident this will",
            "optimal design",
            "no major concerns",
            "well-designed molecule",
            "comprehensive solution",
            "final recommendation"
        ]
        
        # Keywords indicating the agent wants to proceed
        proceed_indicators = [
            "Proceed to final simulation",
            "proceed with simulation",
            "test this molecule",
            "move to docking",
            "ready for autodock",
            "no further modifications needed"
        ]
        
        # Check for confidence indicators
        confidence_score = sum(1 for indicator in confidence_indicators 
                             if indicator in response_lower)
        
        # Check for proceed indicators
        proceed_score = sum(1 for indicator in proceed_indicators 
                          if indicator in response_lower)
        
        # Agent is confident if:
        # 1. High confidence score (2+ confidence indicators), OR
        # 2. Any proceed indicators, OR
        # 3. Response is short and positive (likely means satisfied)
        is_confident = (
            confidence_score >= 2 or
            proceed_score >= 1 or
            (len(response_lower) < 200 and any(word in response_lower 
                                              for word in ["good", "promising", "viable", "acceptable"]))
        )
        
        if is_confident:
            print(f"\n🎯 Agent 3 confidence analysis: {confidence_score} confidence indicators, {proceed_score} proceed indicators")
        
        return is_confident
class SubAgent4SolutionWriter:
    """
    Solution Writer Agent - converts discussion outcomes to AutoDock Vina format
    and runs molecular docking simulations.
    """
    
    def __init__(self, db_manager: DatabaseManager = None):
        self.model = genai.GenerativeModel(MODEL_NAME)
        self.simulation_tool = AutoDockVinaSimulationTool()
        self.db_manager = db_manager
    
    def extract_final_molecule(self, discussion_history: List[str]) -> str:
        """
        Extracts the final molecular structure from discussion.
        
        Args:
            discussion_history: Complete discussion between agents
            
        Returns:
            SMILES string of the final proposed molecule
        """
        full_discussion = "\n".join(discussion_history)
        
        prompt = f"""
Based on the complete discussion below, extract or formulate the FINAL proposed molecular structure.

DISCUSSION:
{full_discussion}

REQUIREMENTS:
1. Identify the most promising molecular design from the discussion
2. Convert it to a valid SMILES string format
3. If multiple molecules were discussed, choose the best one based on the analysis
4. If no specific SMILES was provided, create one based on the described features
5. Ensure the SMILES represents a drug-like molecule suitable for {TARGET_PROTEIN} inhibition

OUTPUT FORMAT:
Provide ONLY the SMILES string, nothing else. No explanations, no formatting, just the SMILES.

Example output: CC(C)C1=NC(=NC(=N1)N2CCN(CC2)C(=O)C3=CC=C(C=C3)F)N4CCOCC4
"""
        
        response = self.model.generate_content(prompt)
        if response.candidates and response.candidates[0].content.parts:
            response_text = response.text.strip()
            smiles = self._extract_valid_smiles(response_text)
        else:
            # Fallback SMILES for a generic kinase inhibitor
            smiles = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"  # Gefitinib-based
        
        # Clean up any remaining formatting
        lines = smiles.split('\n')
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('//'):
                return line
        
        return smiles
    
    def _extract_valid_smiles(self, response_text: str) -> str:
        """Extract and validate SMILES string from response text."""
        try:
            from rdkit import Chem
            
            # Look for SMILES patterns in the response
            import re
            
            # Common SMILES patterns to look for
            smiles_patterns = [
                r'SMILES[:\s]+([A-Za-z0-9\(\)\[\]@\+\-=#\.\\\\/]+)',
                r'([A-Za-z0-9\(\)\[\]@\+\-=#\.\\\\/]{20,})',  # Long chemical strings
                r'([CONSPFClBrI][A-Za-z0-9\(\)\[\]@\+\-=#\.\\\\/]{15,})'  # Starting with common atoms
            ]
            
            for pattern in smiles_patterns:
                matches = re.findall(pattern, response_text)
                for match in matches:
                    # Clean the match
                    candidate_smiles = match.strip().replace("`", "").replace("'", "").replace('"', '')
                    
                    # Validate with RDKit
                    try:
                        mol = Chem.MolFromSmiles(candidate_smiles)
                        if mol is not None:
                            # Additional validation
                            Chem.SanitizeMol(mol)
                            print(f"✅ Valid SMILES found: {candidate_smiles}")
                            return candidate_smiles
                    except:
                        continue
            
            # If no valid SMILES found, return a known good one
            print("⚠️ No valid SMILES found in response, using Gefitinib template")
            return "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
            
        except ImportError:
            # RDKit not available, use simple extraction
            import re
            match = re.search(r'([A-Za-z0-9\(\)\[\]@\+\-=#\.\\\\/]{20,})', response_text)
            if match:
                return match.group(1).strip()
            else:
                return "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
    
    def run_simulation(self, molecule_smiles: str) -> Dict[str, Any]:
        """
        Runs AutoDock Vina simulation for the proposed molecule.
        
        Args:
            molecule_smiles: SMILES string of the molecule
            
        Returns:
            Simulation results including binding energy
        """
        print(f"\n🧪 === SUB-AGENT 4: AUTODOCK VINA SIMULATION ===")
        print(f"🧬 Testing molecule: {molecule_smiles}")
        
        result = self.simulation_tool.run_simulation(molecule_smiles, f"{TARGET_PROTEIN}_{TARGET_MUTATION}")
        
        print(f"📊 Binding Energy: {result['binding_energy']} kcal/mol")
        print(f"✅ Simulation Status: {result['status']}")
        if result.get('output_file'):
            print(f"📁 3D Structure: {result['output_file']}")
        if result.get('visualization_file'):
            print(f"🎨 3D Visualization: {result['visualization_file']}")
        
        # Log solution to database
        if self.db_manager:
            self.db_manager.log_solution(
                molecule_smiles, 
                result['binding_energy'], 
                result, 
                result.get('visualization_file')
            )
        
        return result


class CollaborativeDiscussion:
    """
    Manages the collaborative discussion between Sub-Agent 2 and Sub-Agent 3.
    """
    
    def __init__(self, db_manager: DatabaseManager = None):
        self.agent2 = SubAgent2Chemist(db_manager)
        self.agent3 = SubAgent3Biologist(db_manager)
        self.agent4 = SubAgent4SolutionWriter(db_manager)
        self.db_manager = db_manager
    
    def conduct_discussion(self, context: str) -> tuple[str, List[str]]:
        """
        Conducts collaborative discussion between agents 2 and 3.
        
        Args:
            context: Research context and any feedback from previous iterations
            
        Returns:
            Tuple of (final_molecule_smiles, discussion_history)
        """
        print("💬 SUB-AGENTS 2 & 3: Collaborative Discussion")
        
        discussion_history = [f"🎯 **DISCUSSION CONTEXT:**\n{context}"]
        
        for round_num in range(DISCUSSION_ROUNDS):
            print(f"🔄 Round {round_num + 1}/{DISCUSSION_ROUNDS}")
            
            # Agent 2 (Chemist) proposes or refines
            agent2_response = self.agent2.propose_solution(context, discussion_history)
            print(f"👨‍🔬 Chemist: {agent2_response[:150]}...")
            discussion_history.append(f"\n👨‍🔬 **MEDICINAL CHEMIST (Agent 2):**\n{agent2_response}")
            
            # Agent 3 (Biologist) critiques and validates
            agent3_response = self.agent3.critique_solution(context, discussion_history)
            print(f"👩‍🔬 Biologist: {agent3_response[:150]}...")
            discussion_history.append(f"\n👩‍🔬 **COMPUTATIONAL BIOLOGIST (Agent 3):**\n{agent3_response}")

            # Check if Agent 3 is confident and wants to end early
            if self.agent3.is_confident_about_solution(agent3_response):
                print(f"✅ Discussion concluded after {round_num + 1} rounds")
                break
        
        # Extract final molecule from discussion
        print("🎯 Extracting final molecule...")
        final_molecule = self.agent4.extract_final_molecule(discussion_history)
        print(f"🧬 Final Molecule: {final_molecule}")
        
        return final_molecule, discussion_history