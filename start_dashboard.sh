#!/bin/bash

# NC Multiverse Dashboard Starter Script
# This script starts the interactive dashboard

echo "🚀 Starting NC Multiverse Dashboard..."
echo "📊 Loading data from nc_multiverse.csv..."
echo "🌐 Dashboard will be available at: http://localhost:8050"
echo "⏹️  Press Ctrl+C to stop the dashboard"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed or not in PATH"
    echo "Please install Python 3 and try again"
    exit 1
fi

# Check if the data file exists
if [ ! -f "nc_multiverse.csv" ]; then
    echo "❌ nc_multiverse.csv not found in current directory"
    echo "Please ensure you're in the correct directory"
    exit 1
fi

# Check if main.py exists
if [ ! -f "main.py" ]; then
    echo "❌ main.py not found in current directory"
    echo "Please ensure you're in the correct directory"
    exit 1
fi

echo "✅ All checks passed. Starting dashboard..."
echo ""

# Start the dashboard
python3 main.py
