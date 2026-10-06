@echo off
uv run --project "%~dp0..\.." python "%~dp0..\arrakis\validate_arrakis.py" %*
