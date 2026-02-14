#!/bin/bash
# Setup script for mission_planner with Python 3.10

echo "Setting up Python 3.10 virtual environment for mission_planner..."

# Check if python3.10 is available
if ! command -v python3.10 &> /dev/null; then
    echo "Error: Python 3.10 not found. Please install it first:"
    echo "  Ubuntu/Debian: sudo apt install python3.10 python3.10-venv"
    echo "  Or use pyenv: pyenv install 3.10.13"
    exit 1
fi

# Create virtual environment
python3.10 -m venv venv

# Activate and install dependencies
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

echo ""
echo "✓ Setup complete!"
echo ""
echo "To activate the environment:"
echo "  source venv/bin/activate"
echo ""
echo "To run the app:"
echo "  streamlit run app.py"
