import os
import argparse
from argparse import Namespace

def create_output_dirs(samples_dir:str, output_dir:str) -> None:
    # Get all the subdirectory names in the samples folder
    sample_folders:list[str] = [d for d in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, d))]

    # Create the corresponding directories in first_analysis and second_analysis
    for sample in sample_folders:
        first_analysis_path:str = os.path.join(output_dir, "first_analysis", sample)
        second_analysis_path:str = os.path.join(output_dir, "second_analysis", sample)
        
        # Create directories if they don't exist
        os.makedirs(first_analysis_path, exist_ok=True)
        os.makedirs(second_analysis_path, exist_ok=True)
        
        print(f"Directories created: {first_analysis_path}, {second_analysis_path}")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create output directories for samples.")
    
    parser.add_argument('--samples_dir', type=str, required=True, help='Path to the samples directory')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to the output directory')

    args = parser.parse_args()

    create_output_dirs(args.samples_dir, args.output_dir)

