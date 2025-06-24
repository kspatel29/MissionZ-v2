#!/usr/bin/env python3
"""
Test script to verify database integration and existing solutions retrieval.
"""

import os
import sys

def test_database_connection():
    """Test basic database connection."""
    try:
        from database import DatabaseManager
        
        print("🧪 Testing Database Connection")
        print("=" * 40)
        
        # Initialize database manager
        db_manager = DatabaseManager()
        print("✅ Database manager initialized")
        
        # Test getting research space info
        space_info = db_manager.get_research_space_info()
        if space_info:
            print(f"✅ Research space found: {space_info.get('space_title', 'Untitled')}")
        else:
            print("⚠️ Research space not found (this is okay for testing)")
        
        return db_manager
        
    except Exception as e:
        print(f"❌ Database connection failed: {e}")
        return None

def test_existing_solutions(db_manager):
    """Test retrieving existing solutions."""
    if not db_manager:
        return
    
    print("\n🗃️ Testing Existing Solutions Retrieval")
    print("=" * 45)
    
    try:
        # Get existing solutions summary
        solutions_summary = db_manager.get_existing_solutions_summary()
        print("✅ Solutions summary retrieved:")
        print(solutions_summary[:500] + "..." if len(solutions_summary) > 500 else solutions_summary)
        
        # Test checking if a specific molecule exists
        test_molecule = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"  # Gefitinib
        exists = db_manager.check_molecule_exists(test_molecule)
        print(f"\n✅ Molecule existence check (Gefitinib): {'Found' if exists else 'Not found'}")
        
        # Get solutions history
        history = db_manager.get_solutions_history(limit=5)
        print(f"✅ Recent solutions count: {len(history)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error testing existing solutions: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_logging(db_manager):
    """Test logging functionality."""
    if not db_manager:
        return
    
    print("\n📝 Testing Logging Functionality")
    print("=" * 35)
    
    try:
        # Test conversation logging
        success = db_manager.log_agent_conversation(
            "TestAgent", 
            "This is a test conversation from the test script.",
            iteration=0
        )
        print(f"✅ Conversation logging: {'Success' if success else 'Failed'}")
        
        # Test solution logging (without HTML file)
        test_result = {
            "binding_energy": -8.5,
            "status": "test_simulation",
            "job_id": "test_123"
        }
        
        success = db_manager.log_solution(
            "CCO",  # Simple ethanol for testing
            -8.5,
            test_result
        )
        print(f"✅ Solution logging: {'Success' if success else 'Failed'}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error testing logging: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all database tests."""
    print("🧬 Database Integration Test")
    print("=" * 50)
    
    # Check environment variables
    if not os.getenv("SUPABASE_URL") or not os.getenv("SUPABASE_ANON_KEY"):
        print("❌ Missing Supabase credentials in .env file")
        print("Please add SUPABASE_URL and SUPABASE_ANON_KEY to your .env file")
        return False
    
    tests_passed = 0
    total_tests = 3
    
    # Test 1: Database connection
    db_manager = test_database_connection()
    if db_manager:
        tests_passed += 1
    
    # Test 2: Existing solutions
    if test_existing_solutions(db_manager):
        tests_passed += 1
    
    # Test 3: Logging functionality
    if test_logging(db_manager):
        tests_passed += 1
    
    print(f"\n📊 Test Results: {tests_passed}/{total_tests} passed")
    
    if tests_passed == total_tests:
        print("🎉 All database tests passed!")
        print("✅ Database integration is working correctly")
        return True
    else:
        print("❌ Some database tests failed")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)