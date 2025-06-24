#!/usr/bin/env python3
"""
Health check script for the cancer research pipeline.
"""

import os
import sys
import subprocess
from datetime import datetime

def check_python_imports():
    """Check if all required Python packages are available."""
    required_packages = [
        'google.generativeai',
        'supabase',
        'rdkit',
        'vina',
        'requests'
    ]
    
    print("🐍 Checking Python packages...")
    failed_imports = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"  ✅ {package}")
        except ImportError:
            print(f"  ❌ {package}")
            failed_imports.append(package)
    
    return len(failed_imports) == 0

def check_system_dependencies():
    """Check if system dependencies are available."""
    print("\n🔧 Checking system dependencies...")
    
    dependencies = {
        'obabel': 'Open Babel',
        'python': 'Python'
    }
    
    failed_deps = []
    
    for cmd, name in dependencies.items():
        try:
            result = subprocess.run([cmd, '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                print(f"  ✅ {name}")
            else:
                print(f"  ❌ {name} (command failed)")
                failed_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            print(f"  ❌ {name} (not found)")
            failed_deps.append(name)
    
    return len(failed_deps) == 0

def check_environment_variables():
    """Check if required environment variables are set."""
    print("\n🔑 Checking environment variables...")
    
    required_vars = [
        'GOOGLE_API_KEY',
        'SUPABASE_URL',
        'SUPABASE_ANON_KEY'
    ]
    
    missing_vars = []
    
    for var in required_vars:
        value = os.getenv(var)
        if value:
            print(f"  ✅ {var} (length: {len(value)})")
        else:
            print(f"  ❌ {var} (not set)")
            missing_vars.append(var)
    
    return len(missing_vars) == 0

def check_file_structure():
    """Check if required files and directories exist."""
    print("\n📁 Checking file structure...")
    
    required_paths = [
        'compound_test',
        'compound_test/protein.pdbqt',
        'compound_test/visualizations',
        'config.py',
        'agents.py',
        'tools.py',
        'database.py',
        'main.py'
    ]
    
    missing_paths = []
    
    for path in required_paths:
        if os.path.exists(path):
            print(f"  ✅ {path}")
        else:
            print(f"  ❌ {path} (missing)")
            missing_paths.append(path)
    
    return len(missing_paths) == 0

def check_database_connection():
    """Check database connectivity."""
    print("\n🗄️ Checking database connection...")
    
    try:
        from database import DatabaseManager
        
        db_manager = DatabaseManager()
        space_info = db_manager.get_research_space_info()
        
        print("  ✅ Database connection successful")
        if space_info:
            print(f"  ✅ Research space found: {space_info.get('space_title', 'Untitled')}")
        else:
            print("  ⚠️ Research space not found (this may be normal)")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Database connection failed: {e}")
        return False

def check_autodock_vina():
    """Check AutoDock Vina functionality."""
    print("\n🧪 Checking AutoDock Vina...")
    
    try:
        from tools import AutoDockVinaSimulationTool
        
        vina_tool = AutoDockVinaSimulationTool()
        
        if vina_tool.use_real_vina:
            print("  ✅ AutoDock Vina dependencies available")
            print("  ✅ Real simulations will be used")
        else:
            print("  ⚠️ AutoDock Vina dependencies missing")
            print("  ⚠️ Mock simulations will be used")
        
        return True
        
    except Exception as e:
        print(f"  ❌ AutoDock Vina check failed: {e}")
        return False

def generate_health_report():
    """Generate a comprehensive health report."""
    print("🧬 Cancer Research Pipeline - Health Check")
    print("=" * 50)
    print(f"📅 Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"🐳 Container: {os.getenv('HOSTNAME', 'Unknown')}")
    print(f"👤 User: {os.getenv('USER', 'Unknown')}")
    print(f"📂 Working Directory: {os.getcwd()}")
    
    checks = [
        ("Python Packages", check_python_imports),
        ("System Dependencies", check_system_dependencies),
        ("Environment Variables", check_environment_variables),
        ("File Structure", check_file_structure),
        ("Database Connection", check_database_connection),
        ("AutoDock Vina", check_autodock_vina)
    ]
    
    passed_checks = 0
    total_checks = len(checks)
    
    for check_name, check_func in checks:
        try:
            if check_func():
                passed_checks += 1
        except Exception as e:
            print(f"\n❌ Error in {check_name}: {e}")
    
    print(f"\n📊 Health Check Summary")
    print("=" * 30)
    print(f"✅ Passed: {passed_checks}/{total_checks}")
    print(f"❌ Failed: {total_checks - passed_checks}/{total_checks}")
    
    if passed_checks == total_checks:
        print("\n🎉 All health checks passed!")
        print("✅ System is ready to run the cancer research pipeline")
        return True
    else:
        print(f"\n⚠️ {total_checks - passed_checks} health check(s) failed")
        print("❌ Please fix the issues above before running the pipeline")
        return False

def main():
    """Main health check function."""
    try:
        success = generate_health_report()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n⏹️ Health check interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 Health check crashed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()