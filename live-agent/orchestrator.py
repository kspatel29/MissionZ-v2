"""
Main orchestrator for the news agent system.
"""
import time
import uuid
from typing import Dict, Any
from datetime import datetime

try:
    from .config import (
        GOOGLE_API_KEY, MODEL_NAME, DISCUSSION_ROUNDS, 
        CONFIDENCE_THRESHOLD, NEWS_LOOKBACK_HOURS
    )
    from .agents import (
        SubAgent1Researcher, DiscussionManager, SubAgent4Formatter
    )
    from .database import NewsDatabase
except ImportError:
    from config import (
        GOOGLE_API_KEY, MODEL_NAME, DISCUSSION_ROUNDS, 
        CONFIDENCE_THRESHOLD, NEWS_LOOKBACK_HOURS
    )
    from agents import (
        SubAgent1Researcher, DiscussionManager, SubAgent4Formatter
    )
    from database import NewsDatabase


class NewsAgentOrchestrator:
    """
    Main orchestrator that manages the multi-agent workflow for news processing.
    """
    
    def __init__(self):
        # Initialize database manager
        try:
            self.db_manager = NewsDatabase()
            print("✅ Database connection established")
        except Exception as e:
            print(f"⚠️ Database connection failed: {e}")
            self.db_manager = None
        
        # Initialize agents
        self.researcher = SubAgent1Researcher(self.db_manager)
        self.discussion_manager = DiscussionManager(self.db_manager)
        self.formatter = SubAgent4Formatter(self.db_manager)
        
        print("✅ News Agent Orchestrator initialized")
    
    def process_news_topic(self, topic: str, space_id: str) -> Dict[str, Any]:
        """
        Process a single news topic through the complete agent workflow.
        
        Args:
            topic: News topic to research (e.g., "israel vs iran war")
            space_id: UUID of the news space
            
        Returns:
            Dictionary containing the complete processing results
        """
        print(f"\n🚀 PROCESSING NEWS TOPIC: {topic}")
        print(f"📍 Space ID: {space_id}")
        print("=" * 80)
        
        try:
            # Step 1: Research new news (Sub-Agent 1)
            print("\n📰 STEP 1: News Research")
            research_result = self.researcher.research_news(topic, space_id)
            
            if research_result["status"] != "completed":
                return {
                    "status": "failed",
                    "step": "research",
                    "error": research_result.get("error", "Research failed"),
                    "topic": topic,
                    "space_id": space_id
                }
            
            initial_report = research_result["research_result"]
            print(f"✅ Research completed: {len(initial_report)} characters")
            
            # Step 2: Discussion between agents 2 and 3
            print("\n💬 STEP 2: Bias Analysis & Counter-Claims Discussion")
            discussion_results = self.discussion_manager.conduct_discussion(
                initial_report, 
                space_id, 
                max_rounds=DISCUSSION_ROUNDS,
                confidence_threshold=CONFIDENCE_THRESHOLD
            )
            
            print(f"✅ Discussion completed: {discussion_results['discussion_rounds']} rounds")
            
            # Step 3: Format final report (Sub-Agent 4)
            print("\n📝 STEP 3: Final Report Formatting")
            format_result = self.formatter.format_final_report(discussion_results, space_id)
            
            if format_result["status"] != "completed":
                return {
                    "status": "failed",
                    "step": "formatting",
                    "error": format_result.get("error", "Formatting failed"),
                    "topic": topic,
                    "space_id": space_id
                }
            
            final_report = format_result["formatted_result"]
            print(f"✅ Report formatted successfully")
            
            # Step 4: Save to database
            print("\n💾 STEP 4: Saving to Database")
            if self.db_manager:
                save_success = self.db_manager.save_news_update(space_id, final_report)
                if save_success:
                    print("✅ News update saved to timeline")
                else:
                    print("❌ Failed to save news update")
            else:
                print("⚠️ Database not available - skipping save")
            
            # Return complete results
            return {
                "status": "completed",
                "topic": topic,
                "space_id": space_id,
                "research_result": research_result,
                "discussion_results": discussion_results,
                "final_report": final_report,
                "processing_time": datetime.now().isoformat(),
                "saved_to_db": self.db_manager is not None
            }
            
        except Exception as e:
            error_msg = f"Orchestrator error: {e}"
            print(f"❌ {error_msg}")
            
            return {
                "status": "failed",
                "step": "orchestrator",
                "error": error_msg,
                "topic": topic,
                "space_id": space_id
            }
    
    def run_continuous_loop(self, topic: str, space_id: str, interval_minutes: int = 30):
        """
        Run the news processing in a continuous loop.
        
        Args:
            topic: News topic to monitor
            space_id: UUID of the news space
            interval_minutes: Minutes between processing cycles
        """
        print(f"\n🔄 STARTING CONTINUOUS NEWS MONITORING")
        print(f"📰 Topic: {topic}")
        print(f"📍 Space ID: {space_id}")
        print(f"⏰ Interval: {interval_minutes} minutes")
        print("=" * 80)
        
        cycle_count = 0
        
        try:
            while True:
                cycle_count += 1
                print(f"\n🔄 CYCLE {cycle_count} - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                
                # Process the news topic
                result = self.process_news_topic(topic, space_id)
                
                if result["status"] == "completed":
                    print(f"✅ Cycle {cycle_count} completed successfully")
                    
                    # Print summary of the final report
                    final_report = result["final_report"]
                    print(f"📊 Summary:")
                    print(f"   - Title: {final_report.get('title', 'Unknown')}")
                    print(f"   - Event Time: {final_report.get('event_time', 'Unknown')}")
                    print(f"   - Bias: {final_report.get('bias_percentage', 'Unknown')}%")
                    print(f"   - Credibility: {final_report.get('credibility_score', 'Unknown')}/100")
                    print(f"   - Verification: {final_report.get('verification_status', 'Unknown')}")
                else:
                    print(f"❌ Cycle {cycle_count} failed: {result.get('error', 'Unknown error')}")
                
                # Wait for next cycle
                print(f"⏳ Waiting {interval_minutes} minutes until next cycle...")
                time.sleep(interval_minutes * 60)
                
        except KeyboardInterrupt:
            print(f"\n🛑 Continuous monitoring stopped by user after {cycle_count} cycles")
        except Exception as e:
            print(f"\n❌ Continuous monitoring failed: {e}")
    
    def get_space_summary(self, space_id: str) -> Dict[str, Any]:
        """
        Get a summary of the news space and recent activity.
        
        Args:
            space_id: UUID of the news space
            
        Returns:
            Dictionary containing space summary
        """
        if not self.db_manager:
            return {"error": "Database not available"}
        
        try:
            # Get space info
            space_info = self.db_manager.get_news_space_info(space_id)
            
            # Get recent timeline
            timeline = self.db_manager.get_current_timeline(space_id, NEWS_LOOKBACK_HOURS)
            
            # Get recent conversations
            conversations = self.db_manager.get_conversation_history(space_id, 20)
            
            return {
                "space_info": space_info,
                "timeline_entries": len(timeline),
                "recent_conversations": len(conversations),
                "last_update": timeline[0]["created_at"] if timeline else None,
                "status": "active" if timeline else "inactive"
            }
            
        except Exception as e:
            return {"error": f"Failed to get space summary: {e}"}
