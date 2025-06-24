# 🧬 Cancer Research Agent - VM Deployment Guide

This guide will help you deploy the Cancer Research Agent on your VM using Docker.

## 📋 Prerequisites

### System Requirements
- **OS**: Ubuntu 20.04+ / CentOS 8+ / Any Docker-compatible Linux
- **RAM**: Minimum 4GB, Recommended 8GB+
- **CPU**: Minimum 2 cores, Recommended 4+ cores
- **Storage**: Minimum 10GB free space
- **Network**: Internet connection for API calls

### Required Software
- **Docker** (20.10+)
- **Docker Compose** (1.29+)
- **Git** (for cloning the repository)

## 🚀 Quick Start

### 1. Clone the Repository
```bash
git clone <your-repository-url>
cd research-agent/cure-for-cancer-agent
```

### 2. Set Up Environment Variables
```bash
# Copy the example environment file
cp .env.example .env

# Edit with your API credentials
nano .env
```

Add your credentials:
```bash
GOOGLE_API_KEY=your_google_gemini_api_key_here
SUPABASE_URL=your_supabase_project_url
SUPABASE_ANON_KEY=your_supabase_anon_key
```

### 3. Run the Deployment Script
```bash
./deploy.sh
```

Choose option **1** for full setup and run.

## 🔧 Deployment Options

### Interactive Menu
```bash
./deploy.sh
```

### Command Line Options
```bash
# Full setup and run (recommended for first time)
./deploy.sh full

# Run pipeline only
./deploy.sh run

# Run tests only
./deploy.sh test

# Development mode (interactive shell)
./deploy.sh dev

# Show logs
./deploy.sh logs

# Cleanup containers and volumes
./deploy.sh cleanup
```

## 🐳 Docker Services

### Main Services

1. **cancer-research** - Main pipeline
   ```bash
   docker-compose up cancer-research
   ```

2. **cancer-research-test** - Database tests
   ```bash
   docker-compose --profile test run cancer-research-test
   ```

3. **cancer-research-vina-test** - AutoDock Vina tests
   ```bash
   docker-compose --profile test run cancer-research-vina-test
   ```

4. **cancer-research-dev** - Development shell
   ```bash
   docker-compose --profile dev run cancer-research-dev
   ```

### Manual Docker Commands

```bash
# Build the image
docker-compose build

# Run setup (download protein structure)
docker-compose --profile setup run --rm cancer-research-setup

# Run main pipeline
docker-compose up cancer-research

# Run in background
docker-compose up -d cancer-research

# View logs
docker-compose logs -f cancer-research

# Stop services
docker-compose down

# Clean up everything
docker-compose down -v
docker system prune -f
```

## 📁 Volume Mounts

The following directories are mounted for persistence:

- `./compound_test` → `/app/compound_test` (Protein files and docking results)
- `./logs` → `/app/logs` (Application logs)
- `./.env` → `/app/.env` (Environment variables)
- `cancer_research_data` → `/app/data` (Persistent data volume)

## 🔍 Monitoring and Logs

### View Real-time Logs
```bash
docker-compose logs -f cancer-research
```

### Check Container Status
```bash
docker-compose ps
```

### Access Container Shell
```bash
docker-compose exec cancer-research bash
```

### View Resource Usage
```bash
docker stats cancer-research-agent
```

## 🧪 Testing

### Run All Tests
```bash
./deploy.sh test
```

### Individual Tests
```bash
# Database connectivity
docker-compose --profile test run --rm cancer-research-test

# AutoDock Vina functionality
docker-compose --profile test run --rm cancer-research-vina-test
```

## 🔧 Troubleshooting

### Common Issues

1. **Permission Denied on deploy.sh**
   ```bash
   chmod +x deploy.sh
   ```

2. **Docker not found**
   ```bash
   # Install Docker
   curl -fsSL https://get.docker.com -o get-docker.sh
   sudo sh get-docker.sh
   sudo usermod -aG docker $USER
   # Log out and back in
   ```

3. **Docker Compose not found**
   ```bash
   # Install Docker Compose
   sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
   sudo chmod +x /usr/local/bin/docker-compose
   ```

4. **Environment Variables Not Set**
   - Check `.env` file exists and has correct values
   - Ensure no quotes around values in `.env`

5. **Database Connection Issues**
   - Verify Supabase credentials
   - Check network connectivity
   - Run database test: `./deploy.sh test`

6. **AutoDock Vina Issues**
   - Container includes all dependencies
   - Run Vina test: `docker-compose --profile test run cancer-research-vina-test`

### Debug Mode

Run in development mode for debugging:
```bash
./deploy.sh dev
```

This gives you an interactive shell inside the container.

## 📊 Expected Output

When running successfully, you should see:
```
🧬 SUB-AGENT 1: Research Phase
✅ Research completed
💬 SUB-AGENTS 2 & 3: Collaborative Discussion
🔄 Round 1/2
👨‍🔬 Chemist: Proposing novel EGFR inhibitor...
👩‍🔬 Biologist: Analyzing molecular proposal...
🎯 Extracting final molecule...
🧬 Final Molecule: [SMILES string]
🧪 SUB-AGENT 4: AUTODOCK VINA SIMULATION
📊 Binding Energy: -X.XX kcal/mol
✅ Database connection established
✅ Solution logged to database
```

## 🛑 Stopping the Pipeline

### Graceful Stop
```bash
# In the terminal running the pipeline
Ctrl+C

# Or from another terminal
docker-compose stop cancer-research
```

### Force Stop
```bash
docker-compose down
```

### Complete Cleanup
```bash
./deploy.sh cleanup
```

## 🔄 Updates and Maintenance

### Update the Application
```bash
git pull origin main
docker-compose build
docker-compose up cancer-research
```

### View Database Results
The results are stored in your Supabase database:
- **research_conversation_turns**: All agent conversations
- **solutions**: All tested molecules with binding energies and HTML visualizations

## 📞 Support

If you encounter issues:
1. Check the logs: `docker-compose logs cancer-research`
2. Run tests: `./deploy.sh test`
3. Try development mode: `./deploy.sh dev`
4. Check this troubleshooting guide

## 🎯 Success Criteria

The pipeline is working correctly when:
- ✅ Database connection established
- ✅ AutoDock Vina simulations running
- ✅ Agent conversations logged to database
- ✅ Solutions with binding energies stored
- ✅ HTML visualizations generated and saved

The pipeline will run indefinitely until it achieves a binding energy ≤ -10.0 kcal/mol!