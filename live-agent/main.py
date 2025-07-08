#!/usr/bin/env python3
"""
Main entry point for the news agent system.
"""
import sys
import uuid
import argparse
from datetime import datetime

try:
    from .orchestrator import NewsAgentOrchestrator
    from .config import DISCUSSION_ROUNDS, CONFIDENCE_THRESHOLD, NEWS_LOOKBACK_HOURS
except ImportError:
    from orchestrator import NewsAgentOrchestrator
    from config import DISCUSSION_ROUNDS, CONFIDENCE_THRESHOLD, NEWS_LOOKBACK_HOURS


def main():
    """Main function to run the news agent system."""
    parser = argparse.ArgumentParser(description="News Agent System")
    parser.add_argument("--topic", type=str, required=True, help="News topic to research")
    parser.add_argument("--space-id", type=str, help="News space UUID (will generate if not provided)")
    parser.add_argument("--continuous", action="store_true", help="Run in continuous mode")
    parser.add_argument("--interval", type=int, default=30, help="Interval in minutes for continuous mode")
    parser.add_argument("--single", action="store_true", help="Run single processing cycle")
    parser.add_argument("--summary", action="store_true", help="Show space summary")
    
    args = parser.parse_args()
    
    # Generate space ID if not provided
    if not args.space_id:
        args.space_id = str(uuid.uuid4())
        print(f"📍 Generated new space ID: {args.space_id}")
    
    # Initialize orchestrator
    print("🚀 Initializing News Agent System...")
    orchestrator = NewsAgentOrchestrator()
    
    # Show configuration
    print(f"⚙️ Configuration:")
    print(f"   - Topic: {args.topic}")
    print(f"   - Space ID: {args.space_id}")
    print(f"   - Discussion rounds: {DISCUSSION_ROUNDS}")
    print(f"   - Confidence threshold: {CONFIDENCE_THRESHOLD}")
    print(f"   - News lookback: {NEWS_LOOKBACK_HOURS} hours")
    
    try:
        if args.summary:
            # Show space summary
            print("\n📊 SPACE SUMMARY")
            summary = orchestrator.get_space_summary(args.space_id)
            if "error" in summary:
                print(f"❌ {summary['error']}")
            else:
                print(f"   - Timeline entries: {summary['timeline_entries']}")
                print(f"   - Recent conversations: {summary['recent_conversations']}")
                print(f"   - Last update: {summary['last_update']}")
                print(f"   - Status: {summary['status']}")
        
        elif args.continuous:
            # Run continuous monitoring
            orchestrator.run_continuous_loop(args.topic, args.space_id, args.interval)
        
        elif args.single:
            # Run single processing cycle
            print("\n🔄 SINGLE PROCESSING CYCLE")
            result = orchestrator.process_news_topic(args.topic, args.space_id)
            
            if result["status"] == "completed":
                print("\n✅ PROCESSING COMPLETED SUCCESSFULLY")
                final_report = result["final_report"]
                print(f"\n📊 FINAL REPORT SUMMARY:")
                print(f"   - Title: {final_report.get('title', 'Unknown')}")
                print(f"   - Event Time: {final_report.get('event_time', 'Unknown')}")
                print(f"   - News Report: {final_report.get('news_report', 'No report')[:100]}...")
                print(f"   - Bias: {final_report.get('bias_percentage', 'Unknown')}%")
                print(f"   - Credibility: {final_report.get('credibility_score', 'Unknown')}/100")
                print(f"   - Verification: {final_report.get('verification_status', 'Unknown')}")
                print(f"   - Sources: {len(final_report.get('sources', []))} sources")
                
                # Show full JSON if requested
                import json
                print(f"\n📄 FULL JSON OUTPUT:")
                print(json.dumps(final_report, indent=2))
            else:
                print(f"\n❌ PROCESSING FAILED: {result.get('error', 'Unknown error')}")
                sys.exit(1)
        
        else:
            # Default: show help and run single cycle
            print("\n💡 No mode specified. Running single processing cycle...")
            print("   Use --continuous for continuous monitoring")
            print("   Use --summary to see space summary")
            print("   Use --help for all options")
            
            result = orchestrator.process_news_topic(args.topic, args.space_id)
            
            if result["status"] == "completed":
                print("\n✅ PROCESSING COMPLETED SUCCESSFULLY")
                final_report = result["final_report"]
                print(f"\n📊 SUMMARY:")
                print(f"   - Title: {final_report.get('title', 'Unknown')}")
                print(f"   - Event Time: {final_report.get('event_time', 'Unknown')}")
                print(f"   - Bias: {final_report.get('bias_percentage', 'Unknown')}%")
                print(f"   - Credibility: {final_report.get('credibility_score', 'Unknown')}/100")
                print(f"   - Verification: {final_report.get('verification_status', 'Unknown')}")
            else:
                print(f"\n❌ PROCESSING FAILED: {result.get('error', 'Unknown error')}")
                sys.exit(1)
    
    except KeyboardInterrupt:
        print("\n🛑 Stopped by user")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


def test_workflow():
    """Test function to verify the complete workflow."""
    print("🧪 TESTING NEWS AGENT WORKFLOW")
    print("=" * 50)
    
    # Test topics
    test_topics = [
        "israel vs iran war",
        "ukraine russia conflict",
        "artificial intelligence news"
    ]
    
    orchestrator = NewsAgentOrchestrator()
    
    for topic in test_topics:
        print(f"\n🧪 Testing topic: {topic}")
        space_id = str(uuid.uuid4())
        
        try:
            result = orchestrator.process_news_topic(topic, space_id)
            
            if result["status"] == "completed":
                print(f"✅ Test passed for: {topic}")
                final_report = result["final_report"]
                print(f"   - Credibility: {final_report.get('credibility_score', 'Unknown')}/100")
                print(f"   - Verification: {final_report.get('verification_status', 'Unknown')}")
            else:
                print(f"❌ Test failed for: {topic}")
                print(f"   Error: {result.get('error', 'Unknown')}")
        
        except Exception as e:
            print(f"❌ Test exception for {topic}: {e}")
    
    print("\n🧪 Testing completed")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        test_workflow()
    else:
        main()
