#!/bin/bash

# Cancer Research Agent Deployment Script
# This script sets up and runs the cancer research pipeline on a VM

set -e  # Exit on any error

echo "🧬 Cancer Research Agent Deployment"
echo "=================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ️ $1${NC}"
}

# Check if Docker is installed
check_docker() {
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        echo "Visit: https://docs.docker.com/get-docker/"
        exit 1
    fi
    print_status "Docker is installed"
}

# Check if Docker Compose is installed
check_docker_compose() {
    if ! command -v docker-compose &> /dev/null; then
        print_error "Docker Compose is not installed. Please install Docker Compose first."
        echo "Visit: https://docs.docker.com/compose/install/"
        exit 1
    fi
    print_status "Docker Compose is installed"
}

# Check if .env file exists
check_env_file() {
    if [ ! -f .env ]; then
        print_warning ".env file not found. Creating from template..."
        if [ -f .env.example ]; then
            cp .env.example .env
            print_info "Please edit .env file with your API keys:"
            print_info "  - GOOGLE_API_KEY=your_google_api_key"
            print_info "  - SUPABASE_URL=your_supabase_url"
            print_info "  - SUPABASE_ANON_KEY=your_supabase_anon_key"
            echo ""
            read -p "Press Enter after you've updated the .env file..."
        else
            print_error ".env.example file not found. Please create .env file manually."
            exit 1
        fi
    fi
    print_status ".env file exists"
}

# Validate environment variables
validate_env() {
    source .env
    
    if [ -z "$GOOGLE_API_KEY" ]; then
        print_error "GOOGLE_API_KEY is not set in .env file"
        exit 1
    fi
    
    if [ -z "$SUPABASE_URL" ]; then
        print_error "SUPABASE_URL is not set in .env file"
        exit 1
    fi
    
    if [ -z "$SUPABASE_ANON_KEY" ]; then
        print_error "SUPABASE_ANON_KEY is not set in .env file"
        exit 1
    fi
    
    print_status "Environment variables validated"
}

# Build Docker image
build_image() {
    print_info "Building Docker image..."
    docker-compose build
    print_status "Docker image built successfully"
}

# Run setup (protein download)
run_setup() {
    print_info "Running initial setup (downloading protein structure)..."
    docker-compose --profile setup run --rm cancer-research-setup
    print_status "Setup completed"
}

# Test the system
run_tests() {
    print_info "Running system tests..."
    
    # Test database connection
    print_info "Testing database connection..."
    docker-compose --profile test run --rm cancer-research-test
    
    # Test AutoDock Vina
    print_info "Testing AutoDock Vina..."
    docker-compose --profile test run --rm cancer-research-vina-test
    
    print_status "All tests completed"
}

# Run the main pipeline
run_pipeline() {
    print_info "Starting cancer research pipeline..."
    print_warning "This will run indefinitely until the binding energy goal is achieved."
    print_warning "Press Ctrl+C to stop the pipeline."
    echo ""
    
    # Run the main pipeline
    docker-compose up cancer-research
}

# Interactive development mode
run_dev_mode() {
    print_info "Starting development mode..."
    docker-compose --profile dev run --rm cancer-research-dev
}

# Clean up containers and volumes
cleanup() {
    print_info "Cleaning up containers and volumes..."
    docker-compose down -v
    docker system prune -f
    print_status "Cleanup completed"
}

# Show logs
show_logs() {
    docker-compose logs -f cancer-research
}

# Main menu
show_menu() {
    echo ""
    echo "🧬 Cancer Research Agent - Deployment Options"
    echo "============================================="
    echo "1. 🚀 Full Setup and Run (Recommended for first time)"
    echo "2. 🔬 Run Pipeline Only"
    echo "3. 🧪 Run Tests Only"
    echo "4. 👨‍💻 Development Mode (Interactive Shell)"
    echo "5. 📋 Show Logs"
    echo "6. 🧹 Cleanup"
    echo "7. ❌ Exit"
    echo ""
}

# Main execution
main() {
    # Initial checks
    check_docker
    check_docker_compose
    check_env_file
    validate_env
    
    # Show menu if no arguments provided
    if [ $# -eq 0 ]; then
        while true; do
            show_menu
            read -p "Choose an option (1-7): " choice
            
            case $choice in
                1)
                    build_image
                    run_setup
                    run_tests
                    run_pipeline
                    break
                    ;;
                2)
                    run_pipeline
                    break
                    ;;
                3)
                    run_tests
                    break
                    ;;
                4)
                    run_dev_mode
                    break
                    ;;
                5)
                    show_logs
                    break
                    ;;
                6)
                    cleanup
                    break
                    ;;
                7)
                    print_info "Goodbye!"
                    exit 0
                    ;;
                *)
                    print_error "Invalid option. Please choose 1-7."
                    ;;
            esac
        done
    else
        # Handle command line arguments
        case $1 in
            "setup")
                build_image
                run_setup
                ;;
            "test")
                run_tests
                ;;
            "run")
                run_pipeline
                ;;
            "dev")
                run_dev_mode
                ;;
            "logs")
                show_logs
                ;;
            "cleanup")
                cleanup
                ;;
            "full")
                build_image
                run_setup
                run_tests
                run_pipeline
                ;;
            *)
                echo "Usage: $0 [setup|test|run|dev|logs|cleanup|full]"
                echo "Or run without arguments for interactive menu"
                exit 1
                ;;
        esac
    fi
}

# Run main function
main "$@"