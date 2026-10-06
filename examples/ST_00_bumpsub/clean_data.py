#!/usr/bin/env python3
"""
Cross-platform script to clean result data files

How to execute: python clean_data.py
"""
import os
import glob

def clean_data():
    """
    Remove all .dat files in RESU directory
    """
    # Check if RESU directory exists
    if not os.path.exists("RESU"):
        print("RESU directory not found. Nothing to clean.")
        return
    
    # Find and remove all .dat files
    dat_files = glob.glob(os.path.join("RESU", "*.dat"))
    
    if not dat_files:
        print("No .dat files found in RESU directory.")
        return
    
    removed_count = 0
    for file in dat_files:
        try:
            if os.path.isfile(file):
                os.remove(file)
                removed_count += 1
        except OSError as e:
            print(f"Error removing {file}: {e}")
    
    print(f"Cleaned {removed_count} .dat file(s) from RESU directory.")

if __name__ == '__main__':
    clean_data()
