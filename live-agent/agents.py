"""
Agent implementations for the news agent system.
"""
import google.generativeai as genai
from typing import Dict, Any, List
import json
from datetime import datetime, timedelta

try:
    from .config import GOOGLE_API_KEY, MODEL_NAME, AGENT_NAMES, NEWS_LOOKBACK_HOURS
    from .tools import datetime_tool, google_search_agent, url_context_tool
    from .database import NewsDatabase
except ImportError:
    from config import GOOGLE_API_KEY, MODEL_NAME, AGENT_NAMES, NEWS_LOOKBACK_HOURS
    from tools import datetime_tool, google_search_agent, url_context_tool
    from database import NewsDatabase


class SubAgent1Researcher:
    """
    News Research Agent that finds new news updates from the past 48 hours.
    Gets current timeline data and identifies news that doesn't exist in the timeline.
    Uses specialized agents to avoid tool conflicts.
    """

    def __init__(self, db_manager: NewsDatabase = None):
        # Configure Gemini API
        genai.configure(api_key=GOOGLE_API_KEY)

        # Initialize specialized tool agents
        self.google_search_agent = google_search_agent
        self.url_context_tool = url_context_tool
        self.datetime_tool = datetime_tool

        # Coordinator model (no tools, just coordination)
        self.coordinator_model = genai.GenerativeModel(MODEL_NAME)

        # Database manager for logging and timeline access
        self.db_manager = db_manager
        self.agent_name = AGENT_NAMES["RESEARCHER"]
    
    def research_news(self, topic: str, space_id: str) -> Dict[str, Any]:
        """
        Research new news updates for the given topic.
        
        Args:
            topic: News topic to research (e.g., "israel vs iran war")
            space_id: UUID of the news space
            
        Returns:
            Dictionary containing the news update and metadata
        """
        print(f"🔬 SUB-AGENT 1: Researching news for topic: {topic}")
        
        # Step 1: Get current timeline data
        current_timeline = []
        if self.db_manager:
            current_timeline = self.db_manager.get_current_timeline(space_id, NEWS_LOOKBACK_HOURS)
        
        # Step 2: Prepare timeline summary for context
        timeline_summary = self._prepare_timeline_summary(current_timeline)
        
        # Step 3: Get current time
        current_time = self.datetime_tool.get_formatted_datetime()

        # Step 4: Search for recent news
        search_query = f"latest news {topic} past 48 hours {current_time}"
        search_results = self.google_search_agent.search(search_query)

        # Step 5: Analyze and synthesize findings
        research_prompt = f"""
        You are a news research agent. Based on the search results below, find ONE single news update about "{topic}" that has happened in the past 48 hours and does NOT exist in the current timeline.

        CURRENT TIME: {current_time}

        CURRENT TIMELINE SUMMARY:
        {timeline_summary}

        SEARCH RESULTS:
        {search_results}

        INSTRUCTIONS:
        1. Analyze the search results to find recent news about "{topic}"
        2. Find ONE significant news event/update that is NOT already covered in the timeline
        3. Provide an approximate event time if the exact time is not available
        4. Include source URLs for verification

        REQUIREMENTS:
        - Focus on factual, verifiable news events
        - Ensure the news is from the past 48 hours
        - Avoid duplicating existing timeline entries
        - Include source URLs for verification
        - Provide approximate event time

        Please provide your findings in a clear, structured format.
        """

        try:
            # Generate response using coordinator model
            response = self.coordinator_model.generate_content(research_prompt)

            research_result = response.text if response.text else "No research results generated"
            
            # Print detailed output for debugging
            print(f"📋 RESEARCH RESULT:")
            print(f"   Topic: {topic}")
            print(f"   Timeline entries reviewed: {len(current_timeline)}")
            print(f"   Research output: {research_result[:200]}...")

            # Log the conversation turn
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "topic": topic,
                        "research_result": research_result,
                        "timeline_entries_reviewed": len(current_timeline),
                        "timestamp": datetime.now().isoformat()
                    }
                )
            
            return {
                "agent": self.agent_name,
                "topic": topic,
                "research_result": research_result,
                "timeline_context": timeline_summary,
                "status": "completed"
            }
            
        except Exception as e:
            error_msg = f"Research failed: {e}"
            print(f"❌ {error_msg}")
            
            # Log the error
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "topic": topic,
                        "error": error_msg,
                        "timestamp": datetime.now().isoformat()
                    }
                )
            
            return {
                "agent": self.agent_name,
                "topic": topic,
                "error": error_msg,
                "status": "failed"
            }
    
    def _prepare_timeline_summary(self, timeline_entries: List[Dict[str, Any]]) -> str:
        """
        Prepare a summary of existing timeline entries.
        
        Args:
            timeline_entries: List of timeline entries from database
            
        Returns:
            Formatted summary string
        """
        if not timeline_entries:
            return "No existing timeline entries found."
        
        summary_parts = []
        for entry in timeline_entries[:10]:  # Limit to recent 10 entries
            news_data = entry.get('news', {})
            event_time = entry.get('event_time', 'Unknown time')
            news_report = news_data.get('news_report', 'No report available')
            
            # Truncate long reports
            if len(news_report) > 200:
                news_report = news_report[:200] + "..."
            
            summary_parts.append(f"- {event_time}: {news_report}")
        
        return f"Recent timeline entries ({len(timeline_entries)} total):\n" + "\n".join(summary_parts)


class SubAgent2BiasAnalyzer:
    """
    Bias Analysis Agent that analyzes news reports for bias, credibility, and source verification.
    Uses specialized agents to avoid tool conflicts.
    """

    def __init__(self, db_manager: NewsDatabase = None):
        # Configure Gemini API
        genai.configure(api_key=GOOGLE_API_KEY)

        # Initialize specialized tool agents
        self.google_search_agent = google_search_agent
        self.url_context_tool = url_context_tool
        self.datetime_tool = datetime_tool

        # Coordinator model (no tools, just coordination)
        self.coordinator_model = genai.GenerativeModel(MODEL_NAME)

        # Database manager for logging
        self.db_manager = db_manager
        self.agent_name = AGENT_NAMES["BIAS_ANALYZER"]

    def analyze_bias(self, initial_report: str, space_id: str, discussion_history: List[str] = None) -> Dict[str, Any]:
        """
        Analyze the initial news report for bias, credibility, and source verification.

        Args:
            initial_report: The initial news report from Sub-Agent 1
            space_id: UUID of the news space
            discussion_history: Previous discussion context

        Returns:
            Dictionary containing bias analysis results
        """
        print(f"🔍 SUB-AGENT 2: Analyzing bias in news report")

        # Prepare discussion context
        context = ""
        if discussion_history:
            context = "\n".join(discussion_history[-3:])  # Last 3 exchanges

        # Step 1: Search for verification sources
        # Extract key terms from initial report for search
        key_terms = initial_report[:100].replace('\n', ' ')  # First 100 chars as key terms
        verification_query = f"fact check verify {key_terms} news sources credibility"
        verification_results = self.google_search_agent.search(verification_query)

        # Step 2: Analyze with coordinator model
        analysis_prompt = f"""
        You are a bias analysis expert. Analyze the following news report for bias, credibility, and source verification.

        INITIAL NEWS REPORT:
        {initial_report}

        PREVIOUS DISCUSSION CONTEXT:
        {context}

        VERIFICATION SEARCH RESULTS:
        {verification_results}

        ANALYSIS REQUIREMENTS:
        1. Identify potential bias indicators:
           - Language bias (emotional, loaded terms)
           - Source bias (political leaning, credibility)
           - Selection bias (what's included/excluded)
           - Confirmation bias (cherry-picking facts)

        2. Assess credibility factors:
           - Source reputation and track record
           - Fact-checking and verification
           - Multiple source confirmation
           - Primary vs secondary sources

        3. Provide specific analysis on:
           - Bias percentage (0-100%)
           - Bias direction (if any)
           - Credibility score (0-100%)
           - Source quality assessment
           - Specific bias indicators found
           - Recommendations for verification

        Focus on being thorough and objective in your analysis.
        """

        try:
            # Generate response using coordinator model
            response = self.coordinator_model.generate_content(analysis_prompt)

            analysis_result = response.text if response.text else "No analysis results generated"

            # Print detailed output for debugging
            print(f"📋 BIAS ANALYSIS RESULT:")
            print(f"   Initial report length: {len(initial_report)} chars")
            print(f"   Analysis output: {analysis_result[:200]}...")

            # Log the conversation turn
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "analysis_result": analysis_result,
                        "initial_report_length": len(initial_report),
                        "timestamp": datetime.now().isoformat()
                    }
                )

            return {
                "agent": self.agent_name,
                "analysis_result": analysis_result,
                "status": "completed"
            }

        except Exception as e:
            error_msg = f"Bias analysis failed: {e}"
            print(f"❌ {error_msg}")

            # Log the error
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "error": error_msg,
                        "timestamp": datetime.now().isoformat()
                    }
                )

            return {
                "agent": self.agent_name,
                "error": error_msg,
                "status": "failed"
            }


class SubAgent3CounterClaims:
    """
    Counter-Claims Agent that provides counter-claims, verifies sources, and manages discussion confidence.
    Uses specialized agents to avoid tool conflicts.
    """

    def __init__(self, db_manager: NewsDatabase = None):
        # Configure Gemini API
        genai.configure(api_key=GOOGLE_API_KEY)

        # Initialize specialized tool agents
        self.google_search_agent = google_search_agent
        self.url_context_tool = url_context_tool
        self.datetime_tool = datetime_tool

        # Coordinator model (no tools, just coordination)
        self.coordinator_model = genai.GenerativeModel(MODEL_NAME)

        # Database manager for logging
        self.db_manager = db_manager
        self.agent_name = AGENT_NAMES["COUNTER_CLAIMS"]

    def provide_counter_claims(self, bias_analysis: str, space_id: str, discussion_history: List[str] = None) -> Dict[str, Any]:
        """
        Provide counter-claims and alternative perspectives on the news analysis.

        Args:
            bias_analysis: The bias analysis from Sub-Agent 2
            space_id: UUID of the news space
            discussion_history: Previous discussion context

        Returns:
            Dictionary containing counter-claims and verification results
        """
        print(f"⚖️ SUB-AGENT 3: Providing counter-claims and verification")

        # Prepare discussion context
        context = ""
        if discussion_history:
            context = "\n".join(discussion_history[-5:])  # Last 5 exchanges

        # Step 1: Search for alternative viewpoints
        key_terms = bias_analysis[:100].replace('\n', ' ')  # Extract key terms
        alternative_query = f"alternative viewpoint counter-argument {key_terms}"
        alternative_results = self.google_search_agent.search(alternative_query)

        # Step 2: Generate counter-claims analysis
        counter_claims_prompt = f"""
        You are a counter-claims verification expert. Your role is to provide alternative perspectives, counter-claims, and thorough verification.

        BIAS ANALYSIS FROM SUB-AGENT 2:
        {bias_analysis}

        DISCUSSION HISTORY:
        {context}

        ALTERNATIVE VIEWPOINTS SEARCH RESULTS:
        {alternative_results}

        VERIFICATION REQUIREMENTS:
        1. Provide counter-claims with sources:
           - Alternative interpretations of events
           - Contradictory evidence or reports
           - Different expert opinions
           - Historical context that challenges the narrative

        2. Verify factual claims:
           - Cross-reference with multiple sources
           - Check primary sources when possible
           - Identify unverified claims
           - Note conflicting information

        3. Assess discussion confidence:
           - How well-verified are the claims?
           - Are there significant contradictions?
           - Quality of sources on both sides
           - Overall confidence in the analysis (0.0-1.0)

        4. Determine if discussion should continue:
           - If confidence > 0.8, recommend ending discussion
           - If confidence < 0.8, recommend continuing with specific areas to investigate

        Provide thorough counter-analysis with specific sources and evidence.
        Include a clear confidence score (0.0-1.0) in your response.
        """

        try:
            # Generate response using coordinator model
            response = self.coordinator_model.generate_content(counter_claims_prompt)

            counter_claims_result = response.text if response.text else "No counter-claims results generated"

            # Extract confidence score from the response (simple pattern matching)
            confidence_score = self._extract_confidence_score(counter_claims_result)

            # Print detailed output for debugging
            print(f"📋 COUNTER-CLAIMS RESULT:")
            print(f"   Confidence score extracted: {confidence_score}")
            print(f"   Counter-claims output: {counter_claims_result[:200]}...")

            # Log the conversation turn
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "counter_claims_result": counter_claims_result,
                        "confidence_score": confidence_score,
                        "timestamp": datetime.now().isoformat()
                    }
                )

            return {
                "agent": self.agent_name,
                "counter_claims_result": counter_claims_result,
                "confidence_score": confidence_score,
                "status": "completed"
            }

        except Exception as e:
            error_msg = f"Counter-claims analysis failed: {e}"
            print(f"❌ {error_msg}")

            # Log the error
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "error": error_msg,
                        "timestamp": datetime.now().isoformat()
                    }
                )

            return {
                "agent": self.agent_name,
                "error": error_msg,
                "status": "failed"
            }

    def _extract_confidence_score(self, text: str) -> float:
        """
        Extract confidence score from the response text.

        Args:
            text: Response text to analyze

        Returns:
            Confidence score between 0.0 and 1.0
        """
        import re

        # Look for confidence patterns
        patterns = [
            r'confidence[:\s]+([0-9]*\.?[0-9]+)',
            r'confidence score[:\s]+([0-9]*\.?[0-9]+)',
            r'overall confidence[:\s]+([0-9]*\.?[0-9]+)'
        ]

        for pattern in patterns:
            match = re.search(pattern, text.lower())
            if match:
                try:
                    score = float(match.group(1))
                    # Normalize to 0-1 range if needed
                    if score > 1.0:
                        score = score / 100.0
                    return min(max(score, 0.0), 1.0)
                except ValueError:
                    continue

        # Default confidence if not found
        return 0.5


class SubAgent4Formatter:
    """
    Report Formatter Agent that creates structured JSON output from the discussion results.
    """

    def __init__(self, db_manager: NewsDatabase = None):
        # Configure Gemini API
        genai.configure(api_key=GOOGLE_API_KEY)

        # Initialize coordinator model (no tools needed for formatting)
        self.coordinator_model = genai.GenerativeModel(MODEL_NAME)
        self.datetime_tool = datetime_tool

        # Database manager for logging
        self.db_manager = db_manager
        self.agent_name = AGENT_NAMES["FORMATTER"]

    def format_final_report(self, discussion_results: Dict[str, Any], space_id: str) -> Dict[str, Any]:
        """
        Create a structured JSON output from all the discussion results.

        Args:
            discussion_results: Combined results from all previous agents
            space_id: UUID of the news space

        Returns:
            Dictionary containing the final formatted report
        """
        print(f"📝 SUB-AGENT 4: Formatting final report")

        # Extract data from discussion results
        initial_report = discussion_results.get("initial_report", "")
        bias_analysis = discussion_results.get("bias_analysis", "")
        counter_claims = discussion_results.get("counter_claims", "")
        confidence_score = discussion_results.get("confidence_score", 0.5)

        # Get current time for reference
        current_time = self.datetime_tool.get_current_datetime()

        formatting_prompt = f"""
        You are a report formatting expert. Create a structured JSON output from the following discussion results.

        CURRENT TIME: {current_time}

        INITIAL RESEARCH REPORT:
        {initial_report}

        BIAS ANALYSIS:
        {bias_analysis}

        COUNTER-CLAIMS AND VERIFICATION:
        {counter_claims}

        CONFIDENCE SCORE: {confidence_score}

        FORMAT REQUIREMENTS:
        Create a JSON object with the following structure:
        {{
            "title": "Brief, clear title for the news event",
            "event_time": "ISO timestamp of when the NEWS EVENT occurred (NOT current time - extract from news content)",
            "news_report": "Clear, factual summary of the verified news event",
            "bias": "Description of identified bias (if any)",
            "bias_percentage": "Numeric bias percentage (0-100)",
            "credibility_score": "Numeric credibility score (0-100)",
            "sources": ["array", "of", "source", "URLs"],
            "verification_status": "verified/partially_verified/unverified",
            "confidence_level": "high/medium/low based on confidence score",
            "claims_verified": ["list", "of", "verified", "claims"],
            "claims_disputed": ["list", "of", "disputed", "claims"],
            "additional_context": "Any important context or caveats"
        }}

        FORMATTING GUIDELINES:
        1. Extract the actual event time from the news content (when the event happened, NOT current time)
        2. Create a clear, concise title for the news event
        3. Highlight unverified claims with words like "claim", "alleged", "reported"
        4. Clearly distinguish between verified facts and claims
        5. Include all relevant source URLs
        6. Provide approximate event time if exact time is unavailable
        7. Ensure bias percentage reflects the evidence (hard evidence = low bias)
        8. Make the news_report clear and factual

        Return ONLY the JSON object, no additional text or markdown formatting.
        """

        try:
            # Generate response using coordinator model
            response = self.coordinator_model.generate_content(formatting_prompt)

            formatted_result = response.text if response.text else "{}"

            # Try to parse as JSON to validate
            try:
                # Clean the response - remove markdown formatting if present
                clean_result = formatted_result.strip()
                if clean_result.startswith('```json'):
                    clean_result = clean_result[7:]
                if clean_result.endswith('```'):
                    clean_result = clean_result[:-3]
                clean_result = clean_result.strip()

                json_result = json.loads(clean_result)

                # Ensure required fields exist
                required_fields = ["title", "event_time", "news_report", "bias", "bias_percentage",
                                 "credibility_score", "sources", "verification_status",
                                 "confidence_level", "claims_verified", "claims_disputed", "additional_context"]

                for field in required_fields:
                    if field not in json_result:
                        if field == "title":
                            json_result[field] = "News Update"
                        elif field == "event_time":
                            json_result[field] = datetime.now().isoformat()
                        elif field in ["bias_percentage", "credibility_score"]:
                            json_result[field] = 50
                        elif field in ["sources", "claims_verified", "claims_disputed"]:
                            json_result[field] = []
                        else:
                            json_result[field] = "Unknown"

            except json.JSONDecodeError as e:
                print(f"⚠️ JSON parsing failed: {e}")
                print(f"Raw response: {formatted_result[:500]}...")

                # If JSON parsing fails, create a basic structure
                json_result = {
                    "title": "News Update - Formatting Failed",
                    "event_time": datetime.now().isoformat(),
                    "news_report": "Report formatting failed - raw content available",
                    "bias": "Unknown",
                    "bias_percentage": 50,
                    "credibility_score": 50,
                    "sources": [],
                    "verification_status": "unverified",
                    "confidence_level": "low",
                    "claims_verified": [],
                    "claims_disputed": [],
                    "additional_context": f"Raw formatting result: {formatted_result}"
                }

            # Print detailed output for debugging
            print(f"📋 FORMATTER RESULT:")
            print(f"   Title: {json_result.get('title', 'N/A')}")
            print(f"   Event Time: {json_result.get('event_time', 'N/A')}")
            print(f"   Bias: {json_result.get('bias_percentage', 'N/A')}%")
            print(f"   Credibility: {json_result.get('credibility_score', 'N/A')}/100")
            print(f"   Verification: {json_result.get('verification_status', 'N/A')}")
            print(f"   Sources: {len(json_result.get('sources', []))} sources")

            # Log the conversation turn
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "formatted_result": json_result,
                        "confidence_score": confidence_score,
                        "timestamp": datetime.now().isoformat()
                    }
                )

            return {
                "agent": self.agent_name,
                "formatted_result": json_result,
                "status": "completed"
            }

        except Exception as e:
            error_msg = f"Report formatting failed: {e}"
            print(f"❌ {error_msg}")

            # Create fallback JSON structure
            fallback_result = {
                "title": "News Update - Error",
                "event_time": datetime.now().isoformat(),
                "news_report": "Report generation failed",
                "bias": "Unknown",
                "bias_percentage": 50,
                "credibility_score": 0,
                "sources": [],
                "verification_status": "unverified",
                "confidence_level": "low",
                "claims_verified": [],
                "claims_disputed": [],
                "additional_context": f"Error: {error_msg}"
            }

            # Log the error
            if self.db_manager:
                self.db_manager.log_conversation_turn(
                    space_id,
                    self.agent_name,
                    {
                        "error": error_msg,
                        "fallback_result": fallback_result,
                        "timestamp": datetime.now().isoformat()
                    }
                )

            return {
                "agent": self.agent_name,
                "formatted_result": fallback_result,
                "error": error_msg,
                "status": "failed"
            }


class DiscussionManager:
    """
    Manages the discussion loop between Sub-Agent 2 and Sub-Agent 3.
    """

    def __init__(self, db_manager: NewsDatabase = None):
        self.agent2 = SubAgent2BiasAnalyzer(db_manager)
        self.agent3 = SubAgent3CounterClaims(db_manager)
        self.agent4 = SubAgent4Formatter(db_manager)
        self.db_manager = db_manager

    def conduct_discussion(self, initial_report: str, space_id: str, max_rounds: int = 3, confidence_threshold: float = 0.8) -> Dict[str, Any]:
        """
        Conduct discussion between agents 2 and 3 until confidence threshold is reached.

        Args:
            initial_report: Initial news report from Sub-Agent 1
            space_id: UUID of the news space
            max_rounds: Maximum number of discussion rounds
            confidence_threshold: Confidence threshold to exit discussion

        Returns:
            Dictionary containing final discussion results
        """
        print(f"💬 DISCUSSION MANAGER: Starting discussion (max {max_rounds} rounds, threshold {confidence_threshold})")

        discussion_history = []
        confidence_score = 0.0
        round_num = 0

        # Add initial report to history
        discussion_history.append(f"INITIAL REPORT:\n{initial_report}")

        while round_num < max_rounds and confidence_score < confidence_threshold:
            round_num += 1
            print(f"🔄 Discussion Round {round_num}/{max_rounds}")

            # Agent 2 (Bias Analyzer) analyzes
            bias_result = self.agent2.analyze_bias(initial_report, space_id, discussion_history)
            if bias_result["status"] == "completed":
                bias_analysis = bias_result["analysis_result"]
                discussion_history.append(f"\nROUND {round_num} - BIAS ANALYZER:\n{bias_analysis}")
                print(f"✅ Bias analysis completed")
            else:
                print(f"❌ Bias analysis failed: {bias_result.get('error', 'Unknown error')}")
                break

            # Agent 3 (Counter-Claims) responds
            counter_result = self.agent3.provide_counter_claims(bias_analysis, space_id, discussion_history)
            if counter_result["status"] == "completed":
                counter_claims = counter_result["counter_claims_result"]
                confidence_score = counter_result["confidence_score"]
                discussion_history.append(f"\nROUND {round_num} - COUNTER-CLAIMS:\n{counter_claims}")
                print(f"✅ Counter-claims completed (confidence: {confidence_score:.2f})")
            else:
                print(f"❌ Counter-claims failed: {counter_result.get('error', 'Unknown error')}")
                break

            # Check if confidence threshold is reached
            if confidence_score >= confidence_threshold:
                print(f"🎯 Confidence threshold reached ({confidence_score:.2f} >= {confidence_threshold})")
                break
            else:
                print(f"🔄 Continuing discussion (confidence: {confidence_score:.2f} < {confidence_threshold})")

        # Prepare final discussion results
        final_results = {
            "initial_report": initial_report,
            "bias_analysis": bias_result.get("analysis_result", "") if 'bias_result' in locals() else "",
            "counter_claims": counter_result.get("counter_claims_result", "") if 'counter_result' in locals() else "",
            "confidence_score": confidence_score,
            "discussion_rounds": round_num,
            "discussion_history": discussion_history,
            "threshold_reached": confidence_score >= confidence_threshold
        }

        print(f"📊 Discussion completed: {round_num} rounds, final confidence: {confidence_score:.2f}")

        return final_results
