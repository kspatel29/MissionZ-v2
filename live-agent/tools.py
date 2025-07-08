"""
Tools for the news agent system.
"""
import requests
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from typing import Dict, Any
import google.generativeai as genai
from google.generativeai.types import Tool

try:
    from .config import GOOGLE_API_KEY, MODEL_NAME
except ImportError:
    from config import GOOGLE_API_KEY, MODEL_NAME


class DateTimeTool:
    """Tool for getting current date and time."""
    
    def get_current_datetime(self) -> str:
        """
        Get the current date and time in ISO format.
        
        Returns:
            Current datetime string in ISO format with timezone
        """
        current_time = datetime.now(timezone.utc)
        return current_time.isoformat()
    
    def get_formatted_datetime(self) -> str:
        """
        Get the current date and time in a human-readable format.
        
        Returns:
            Formatted datetime string
        """
        current_time = datetime.now(timezone.utc)
        return current_time.strftime("%Y-%m-%d %H:%M:%S UTC")


class GoogleSearchAgent:
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


class URLContextTool:
    """Tool for fetching and analyzing content from URLs."""
    
    def __init__(self):
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
    
    def fetch_url_content(self, url: str) -> str:
        """
        Fetch content from a URL and return text content.
        
        Args:
            url: URL to fetch content from
            
        Returns:
            Text content from the URL or error message
        """
        try:
            response = requests.get(url, headers=self.headers, timeout=30)
            response.raise_for_status()
            
            # Basic text extraction (you might want to use BeautifulSoup for better parsing)
            content = response.text
            
            # Simple text extraction - remove HTML tags
            import re
            text_content = re.sub(r'<[^>]+>', '', content)
            text_content = re.sub(r'\s+', ' ', text_content).strip()
            
            # Limit content length
            if len(text_content) > 5000:
                text_content = text_content[:5000] + "... [content truncated]"
            
            return text_content
            
        except requests.RequestException as e:
            return f"Error fetching URL content: {e}"
        except Exception as e:
            return f"Error processing URL content: {e}"
    
    def analyze_url_content(self, url: str, analysis_prompt: str = "") -> str:
        """
        Fetch URL content and analyze it using Gemini.
        
        Args:
            url: URL to analyze
            analysis_prompt: Specific analysis instructions
            
        Returns:
            Analysis results
        """
        content = self.fetch_url_content(url)
        
        if content.startswith("Error"):
            return content
        
        try:
            model = genai.GenerativeModel(MODEL_NAME)
            
            prompt = f"""
            Analyze the following web content from URL: {url}
            
            {analysis_prompt if analysis_prompt else "Provide a summary of the key information, focusing on factual content, sources, and credibility indicators."}
            
            Content:
            {content}
            """
            
            response = model.generate_content(prompt)
            return response.text if response.text else "Analysis failed: No response generated"
            
        except Exception as e:
            return f"Analysis failed: {e}"


# Tool instances that can be imported
datetime_tool = DateTimeTool()
google_search_agent = GoogleSearchAgent()
url_context_tool = URLContextTool()
