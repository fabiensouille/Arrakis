@echo off
uv run --project "%~dp0..\.." python "%~dp0..\arrakis\clean_data.py" %*
