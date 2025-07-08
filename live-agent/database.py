"""
Supabase database integration for the news agent system.
"""
import os
import json
import uuid
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List
from supabase import create_client, Client
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

class NewsDatabase:
    """Manages database operations for the news agent system."""
    
    def __init__(self):
        # Get Supabase credentials from environment
        self.supabase_url = os.getenv("SUPABASE_URL")
        self.supabase_key = os.getenv("SUPABASE_KEY")
        
        if not self.supabase_url or not self.supabase_key:
            raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
        
        # Initialize Supabase client
        self.supabase: Client = create_client(self.supabase_url, self.supabase_key)
        print("✅ Database connection established")
    
    def log_conversation_turn(self, space_id: str, agent_name: str, agent_response: Dict[str, Any]) -> bool:
        """
        Log a conversation turn to the news_conversation_turns table.
        
        Args:
            space_id: UUID of the news space
            agent_name: Name of the agent
            agent_response: Agent's response data
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Prepare the data
            conversation_data = {
                "space_id": space_id,
                "agent_name": agent_name,
                "agent_response": agent_response
            }
            
            # Insert into database
            result = self.supabase.table("news_conversation_turns").insert(conversation_data).execute()
            
            if result.data:
                print(f"✅ Logged conversation turn for {agent_name}")
                return True
            else:
                print(f"❌ Failed to log conversation for {agent_name}")
                return False
                
        except Exception as e:
            print(f"❌ Database error logging conversation: {e}")
            return False
    
    def get_current_timeline(self, space_id: str, hours_back: int = 48) -> List[Dict[str, Any]]:
        """
        Get current timeline data for a space from the past specified hours.
        
        Args:
            space_id: UUID of the news space
            hours_back: Number of hours to look back (default 48)
            
        Returns:
            List of timeline entries
        """
        try:
            # Calculate the cutoff time
            cutoff_time = datetime.now() - timedelta(hours=hours_back)
            
            result = self.supabase.table("news_timeline")\
                .select("*")\
                .eq("space_id", space_id)\
                .gte("created_at", cutoff_time.isoformat())\
                .order("created_at", desc=True)\
                .execute()
            
            return result.data if result.data else []
            
        except Exception as e:
            print(f"❌ Database error getting timeline: {e}")
            return []
    
    def save_news_update(self, space_id: str, news_data: Dict[str, Any]) -> bool:
        """
        Save a new news update to the news_timeline table.
        
        Args:
            space_id: UUID of the news space
            news_data: Dictionary containing event_time, news_report, bias, etc.
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Prepare the data
            timeline_data = {
                "space_id": space_id,
                "news": news_data,
                "event_time": news_data.get("event_time")
            }
            
            # Insert into database
            result = self.supabase.table("news_timeline").insert(timeline_data).execute()
            
            if result.data:
                news_id = result.data[0]['id']
                print(f"✅ News update saved to timeline (ID: {news_id})")
                return True
            else:
                print(f"❌ Failed to save news update")
                return False
                
        except Exception as e:
            print(f"❌ Database error saving news update: {e}")
            return False
    
    def get_news_space_info(self, space_id: str) -> Optional[Dict[str, Any]]:
        """Get information about a specific news space."""
        try:
            result = self.supabase.table("news_space").select("*").eq("id", space_id).execute()
            
            if result.data:
                return result.data[0]
            else:
                return None
                
        except Exception as e:
            print(f"❌ Database error getting space info: {e}")
            return None
    
    def get_conversation_history(self, space_id: str, limit: int = 50) -> List[Dict[str, Any]]:
        """Get recent conversation history for a news space."""
        try:
            result = self.supabase.table("news_conversation_turns")\
                .select("*")\
                .eq("space_id", space_id)\
                .order("created_at", desc=True)\
                .limit(limit)\
                .execute()
            
            return result.data if result.data else []
            
        except Exception as e:
            print(f"❌ Database error getting conversation history: {e}")
            return []
