#!/usr/bin/env python3
"""
Script to check multiple genome files in a directory for chromosomes over 32 Mbp
Each file represents a complete genome with multiple chromosomes
"""

from Bio import SeqIO
import os
import glob

def check_genome_files(directory_path, threshold_mbp=32):
    """
    Check all genome FASTA files in a directory for chromosomes above threshold
    Each file represents a complete genome with multiple chromosomes
    
    Args:
        directory_path (str): Path to directory containing genome files
        threshold_mbp (int): Threshold in megabase pairs (default: 32)
    """
    
    threshold_bp = threshold_mbp * 1000000  # Convert Mbp to bp
    
    # Find all FASTA files in the directory
    fasta_patterns = ["*.fna", "*.fasta", "*.fa"]
    fasta_files = []
    
    for pattern in fasta_patterns:
        fasta_files.extend(glob.glob(os.path.join(directory_path, pattern)))
    
    if not fasta_files:
        print(f"No genome FASTA files found in directory: {directory_path}")
        return
    
    print(f"Found {len(fasta_files)} genome file(s) in directory")
    print(f"Checking each genome for chromosomes above {threshold_mbp} Mbp ({threshold_bp:,} bp)")
    print("=" * 80)
    
    genomes_with_large_chromosomes = []
    
    # Process each genome file
    for file_path in sorted(fasta_files):
        filename = os.path.basename(file_path)
        print(f"\n GENOME: {filename}")
        print("-" * 60)
        
        try:
            large_chromosomes = []
            total_chromosomes = 0
            
            # Check each chromosome in the genome
            for record in SeqIO.parse(file_path, "fasta"):
                total_chromosomes += 1
                seq_length = len(record.seq)
                seq_length_mbp = seq_length / 1000000
                
                if seq_length > threshold_bp:
                    large_chromosomes.append({
                        'id': record.id,
                        'description': record.description,
                        'length_bp': seq_length,
                        'length_mbp': seq_length_mbp
                    })
                    print(f"   {record.id}: {seq_length:,} bp ({seq_length_mbp:.2f} Mbp) - ABOVE THRESHOLD")
                else:
                    print(f"     {record.id}: {seq_length:,} bp ({seq_length_mbp:.2f} Mbp)")
            
            # Summary for this genome
            if large_chromosomes:
                genomes_with_large_chromosomes.append({
                    'genome_name': filename,
                    'filepath': file_path,
                    'large_chromosomes': large_chromosomes,
                    'total_chromosomes': total_chromosomes
                })
                print(f"   RESULT: {len(large_chromosomes)} out of {total_chromosomes} chromosomes above {threshold_mbp} Mbp")
            else:
                print(f"   RESULT: No chromosomes above {threshold_mbp} Mbp (total: {total_chromosomes} chromosomes)")
                
        except Exception as e:
            print(f"   ERROR: Could not process {filename}: {e}")
    
    # Final summary
    print("\n" + "=" * 80)
    print(" FINAL SUMMARY - GENOMES WITH LARGE CHROMOSOMES")
    print("=" * 80)
    
    if genomes_with_large_chromosomes:
        print(f"\n GENOMES with chromosomes ABOVE {threshold_mbp} Mbp:")
        print("-" * 60)
        
        for genome_info in genomes_with_large_chromosomes:
            print(f"\n {genome_info['genome_name']}")
            print(f"    Path: {genome_info['filepath']}")
            print(f"    Large chromosomes: {len(genome_info['large_chromosomes'])} out of {genome_info['total_chromosomes']}")
            
            for chr_info in genome_info['large_chromosomes']:
                print(f"       {chr_info['id']}: {chr_info['length_mbp']:.2f} Mbp")
        
        print(f"\n TOTAL RESULT: {len(genomes_with_large_chromosomes)} genome(s) have chromosomes above {threshold_mbp} Mbp")
        
        # Create a simple list of genome names for easy reference
        print(f"\n GENOME NAMES WITH LARGE CHROMOSOMES:")
        for genome_info in genomes_with_large_chromosomes:
            print(f"    {genome_info['genome_name']}")
            
    else:
        print(f"\n RESULT: No genomes found with chromosomes above {threshold_mbp} Mbp")
        print("All genomes have chromosomes smaller than the threshold.")

def main():
    """
    Main function - set your directory path here
    """
    # MODIFY THIS PATH to point to your directory containing chromosome files
    directory_path = "/data/home/bt24090/viralrecall/Porifera/chromosomes_only"
    
    # Optional: change the threshold (default is 32 Mbp)
    threshold_mbp = 32
    
    # Check if directory exists
    if not os.path.exists(directory_path):
        print(f"Error: Directory '{directory_path}' not found!")
        print("Please modify the 'directory_path' variable in the main() function.")
        return
    
    if not os.path.isdir(directory_path):
        print(f"Error: '{directory_path}' is not a directory!")
        return
    
    # Run the analysis
    check_genome_files(directory_path, threshold_mbp)

if __name__ == "__main__":
    main()
