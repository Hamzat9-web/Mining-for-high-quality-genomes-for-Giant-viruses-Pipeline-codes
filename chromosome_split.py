
from Bio import SeqIO
import os
import glob
import shutil
import math

def get_species_name(filename):
    """
    Extract species name from filename
    
    Args:
        filename (str): Genome filename
        
    Returns:
        str: Species name for directory creation
    """
    # Remove file extension
    base_name = os.path.splitext(filename)[0]
    
    # Try to extract species name (assuming format: Genus_species_...)
    parts = base_name.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    else:
        return base_name

def calculate_split_parts(length_bp, target_size_mbp=30):
    """
    Calculate how many parts to split a chromosome into
    
    Args:
        length_bp (int): Chromosome length in base pairs
        target_size_mbp (int): Target size per part in Mbp
        
    Returns:
        int: Number of parts to split into
    """
    target_size_bp = target_size_mbp * 1000000
    parts = math.ceil(length_bp / target_size_bp)
    return max(2, parts)  # Minimum 2 parts

def extract_chromosome_number(record):
    """
    Extract chromosome number from record description
    
    Args:
        record: BioPython SeqRecord
        
    Returns:
        str: Chromosome number (e.g., "1", "15") or original ID if not found
    """
    description = record.description.lower()
    
    # Look for "chromosome: X" pattern
    if "chromosome:" in description:
        parts = description.split("chromosome:")
        if len(parts) > 1:
            # Get the part after "chromosome:" and extract the number
            chr_part = parts[1].strip().split()[0]  # Get first word after colon
            # Remove any non-digit characters and return just the number
            chr_number = ''.join(filter(str.isdigit, chr_part))
            if chr_number:
                return f"chr{chr_number}"
    
    # Fallback to original ID if chromosome number not found
    return record.id
    """
    Calculate how many parts to split a chromosome into
    
    Args:
        length_bp (int): Chromosome length in base pairs
        target_size_mbp (int): Target size per part in Mbp
        
    Returns:
        int: Number of parts to split into
    """
    target_size_bp = target_size_mbp * 1000000
    parts = math.ceil(length_bp / target_size_bp)
    return max(2, parts)  # Minimum 2 parts

def split_chromosome(record, num_parts):
    """
    Split a chromosome record into specified number of parts
    
    Args:
        record: BioPython SeqRecord
        num_parts (int): Number of parts to split into
        
    Returns:
        list: List of split SeqRecord objects
    """
    seq_length = len(record.seq)
    part_size = seq_length // num_parts
    split_records = []
    
    for i in range(num_parts):
        start = i * part_size
        if i == num_parts - 1:  # Last part gets remainder
            end = seq_length
        else:
            end = (i + 1) * part_size
        
        # Create new record for this part
        part_record = record[start:end]
        part_record.id = f"{record.id}_part{i+1}"
        part_record.description = f"{record.description} (part {i+1} of {num_parts})"
        split_records.append(part_record)
    
    return split_records

def process_genome_file(file_path, output_base_dir, threshold_mbp=32, target_size_mbp=30):
    """
    Process a single genome file - identify and split large chromosomes into separate files
    
    Args:
        file_path (str): Path to genome file
        output_base_dir (str): Base directory for outputs
        threshold_mbp (int): Threshold for splitting (default: 32 Mbp)
        target_size_mbp (int): Target size for split parts (default: 30 Mbp)
        
    Returns:
        dict: Processing results
    """
    filename = os.path.basename(file_path)
    species_name = get_species_name(filename)
    threshold_bp = threshold_mbp * 1000000
    
    print(f"\n PROCESSING: {filename}")
    print(f"   Species: {species_name}")
    print("-" * 80)
    
    # Create species-specific output directory
    species_output_dir = os.path.join(output_base_dir, species_name)
    os.makedirs(species_output_dir, exist_ok=True)
    
    try:
        # Analyze genome and identify chromosomes needing splitting
        large_chromosomes = []
        small_chromosomes = []
        split_operations = []
        output_files_created = []
        
        for record in SeqIO.parse(file_path, "fasta"):
            seq_length = len(record.seq)
            seq_length_mbp = seq_length / 1000000
            
            if seq_length > threshold_bp:
                # Calculate how many parts needed
                num_parts = calculate_split_parts(seq_length, target_size_mbp)
                part_size_mbp = seq_length_mbp / num_parts
                
                large_chromosomes.append({
                    'record': record,
                    'length_mbp': seq_length_mbp,
                    'num_parts': num_parts,
                    'part_size_mbp': part_size_mbp
                })
                
                print(f"    {record.id}: {seq_length:,} bp ({seq_length_mbp:.2f} Mbp)")
                print(f"      → Will split into {num_parts} parts (~{part_size_mbp:.2f} Mbp each)")
                
                # Split the chromosome and save each part as separate file
                split_records = split_chromosome(record, num_parts)
                
                # Extract chromosome number for cleaner naming
                chr_identifier = extract_chromosome_number(record)
                
                for i, split_record in enumerate(split_records):
                    # New naming convention: [species_name]_[chr_number]_[part_number]
                    part_filename = f"{species_name}_{chr_identifier}_part{i+1}.fasta"
                    part_filepath = os.path.join(species_output_dir, part_filename)
                    SeqIO.write([split_record], part_filepath, "fasta")
                    output_files_created.append(part_filename)
                    
                    part_size_mbp = len(split_record.seq) / 1000000
                    print(f"       Created: {part_filename} ({part_size_mbp:.2f} Mbp)")
                
                split_operations.append({
                    'original_id': record.id,
                    'original_size_mbp': seq_length_mbp,
                    'num_parts': num_parts,
                    'split_files': [f"{species_name}_{chr_identifier}_part{i+1}.fasta" for i in range(num_parts)]
                })
                
            else:
                # Keep chromosome as-is for other_chromosomes file
                small_chromosomes.append(record)
                print(f"   {record.id}: {seq_length:,} bp ({seq_length_mbp:.2f} Mbp) - keeping intact")
        
        if not large_chromosomes:
            print("   No chromosomes above threshold - skipping this genome")
            # Remove empty directory
            os.rmdir(species_output_dir)
            return {
                'filename': filename,
                'species': species_name,
                'processed': False,
                'reason': 'No large chromosomes',
                'large_chromosomes': 0
            }
        
        # Create other_chromosomes.fasta file with all small chromosomes
        if small_chromosomes:
            other_chromosomes_file = "other_chromosomes.fasta"
            other_chromosomes_path = os.path.join(species_output_dir, other_chromosomes_file)
            SeqIO.write(small_chromosomes, other_chromosomes_path, "fasta")
            output_files_created.append(other_chromosomes_file)
            print(f"   Created: {other_chromosomes_file} ({len(small_chromosomes)} chromosomes)")
        
        # Create summary report
        report_path = os.path.join(species_output_dir, f"{species_name}_splitting_report.txt")
        with open(report_path, 'w') as report:
            report.write(f"Chromosome Splitting Report for {species_name}\n")
            report.write("=" * 60 + "\n\n")
            report.write(f"Original file: {filename}\n")
            report.write(f"Threshold: {threshold_mbp} Mbp\n")
            report.write(f"Target split size: ~{target_size_mbp} Mbp\n\n")
            
            report.write("Output Files Created:\n")
            report.write("-" * 30 + "\n")
            for file in sorted(output_files_created):
                report.write(f"• {file}\n")
            
            report.write(f"\nSplit Operations:\n")
            report.write("-" * 20 + "\n")
            for op in split_operations:
                report.write(f"\n• {op['original_id']}: {op['original_size_mbp']:.2f} Mbp\n")
                report.write(f"  → Split into {op['num_parts']} files:\n")
                for file in op['split_files']:
                    report.write(f"    - {file}\n")
            
            if small_chromosomes:
                report.write(f"\nIntact Chromosomes (< {threshold_mbp} Mbp):\n")
                report.write("-" * 35 + "\n")
                for record in small_chromosomes:
                    size_mbp = len(record.seq) / 1000000
                    report.write(f"• {record.id}: {size_mbp:.2f} Mbp\n")
            
            report.write(f"\nSummary:\n")
            report.write("-" * 10 + "\n")
            report.write(f"Total output files: {len(output_files_created)}\n")
            report.write(f"Split chromosome files: {sum(len(op['split_files']) for op in split_operations)}\n")
            report.write(f"Other chromosomes file: {'Yes' if small_chromosomes else 'No'}\n")
            report.write(f"Large chromosomes split: {len(large_chromosomes)}\n")
            report.write(f"Small chromosomes kept intact: {len(small_chromosomes)}\n")
        
        print(f"   Successfully processed - output saved to: {species_output_dir}/")
        print(f"      Total files created: {len(output_files_created)}")
        print(f"      Report: {species_name}_splitting_report.txt")
        
        return {
            'filename': filename,
            'species': species_name,
            'processed': True,
            'output_dir': species_output_dir,
            'output_files': output_files_created,
            'large_chromosomes': len(large_chromosomes),
            'split_operations': len(split_operations),
            'small_chromosomes': len(small_chromosomes),
            'total_output_files': len(output_files_created)
        }
        
    except Exception as e:
        print(f"   ERROR processing {filename}: {e}")
        # Clean up on error
        if os.path.exists(species_output_dir):
            shutil.rmtree(species_output_dir)
        
        return {
            'filename': filename,
            'species': species_name,
            'processed': False,
            'reason': f"Error: {e}",
            'large_chromosomes': 0
        }

def adaptive_chromosome_splitter(input_directory, output_base_directory="split_genomes", 
                                threshold_mbp=32, target_size_mbp=30, delete_originals=True):
    """
    Main function to process all genome files in a directory
    
    Args:
        input_directory (str): Directory containing genome files
        output_base_directory (str): Base directory for outputs
        threshold_mbp (int): Threshold for splitting chromosomes
        target_size_mbp (int): Target size for split parts
        delete_originals (bool): Whether to delete original files after processing
    """
    
    print(" ADAPTIVE CHROMOSOME SPLITTER")
    print("=" * 80)
    print(f"Input directory: {input_directory}")
    print(f"Output directory: {output_base_directory}")
    print(f"Split threshold: {threshold_mbp} Mbp")
    print(f"Target split size: ~{target_size_mbp} Mbp")
    print(f"Delete originals: {'Yes' if delete_originals else 'No'}")
    print("=" * 80)
    
    # Find all genome files
    fasta_patterns = ["*.fna", "*.fasta", "*.fa"]
    fasta_files = []
    
    for pattern in fasta_patterns:
        fasta_files.extend(glob.glob(os.path.join(input_directory, pattern)))
    
    if not fasta_files:
        print(f" No genome files found in directory: {input_directory}")
        return
    
    print(f"Found {len(fasta_files)} genome file(s) to process")
    
    # Create output base directory
    os.makedirs(output_base_directory, exist_ok=True)
    
    # Process each genome file
    results = []
    successful_deletions = []
    
    for file_path in sorted(fasta_files):
        result = process_genome_file(file_path, output_base_directory, threshold_mbp, target_size_mbp)
        results.append(result)
        
        # Delete original file if processing was successful
        if result['processed'] and delete_originals:
            try:
                os.remove(file_path)
                successful_deletions.append(result['filename'])
                print(f"    Deleted original file: {result['filename']}")
            except Exception as e:
                print(f"    Could not delete {result['filename']}: {e}")
    
    # Final Summary
    print("\n" + "=" * 80)
    print(" FINAL PROCESSING SUMMARY")
    print("=" * 80)
    
    processed_genomes = [r for r in results if r['processed']]
    failed_genomes = [r for r in results if not r['processed']]
    
    if processed_genomes:
        print(f"\n SUCCESSFULLY PROCESSED GENOMES: {len(processed_genomes)}")
        print("-" * 60)
        
        total_splits = 0
        for result in processed_genomes:
            total_splits += result['split_operations']
            print(f" {result['species']}")
            print(f"    Output: {result['output_dir']}")
            print(f"    Large chromosomes split: {result['split_operations']}")
            print(f"    Total output files: {result['total_output_files']}")
            print(f"    Small chromosomes kept: {result['small_chromosomes']}")
        
        print(f"\n STATISTICS:")
        print(f"   • Total genomes processed: {len(processed_genomes)}")
        print(f"   • Total chromosomes split: {total_splits}")
        print(f"   • Original files deleted: {len(successful_deletions)}")
    
    if failed_genomes:
        print(f"\n FAILED TO PROCESS: {len(failed_genomes)}")
        for result in failed_genomes:
            print(f"   • {result['filename']}: {result['reason']}")
    
    if not processed_genomes:
        print(" No genomes were processed - no chromosomes above threshold found")

def main():
    """
    Main function - configure your settings here
    """
    # CONFIGURE THESE SETTINGS
    input_directory = ""
    output_directory = ""
    threshold_mbp = 32          # Split chromosomes above this size
    target_size_mbp = 30        # Target size for split parts
    delete_originals = True     # Set to False to keep original files
    
    # Validate input directory
    if not os.path.exists(input_directory):
        print(f" Error: Input directory '{input_directory}' not found!")
        print("Please modify the 'input_directory' variable in the main() function.")
        return
    
    if not os.path.isdir(input_directory):
        print(f" Error: '{input_directory}' is not a directory!")
        return
    
    # Run the adaptive splitter
    adaptive_chromosome_splitter(
        input_directory=input_directory,
        output_base_directory=output_directory,
        threshold_mbp=threshold_mbp,
        target_size_mbp=target_size_mbp,
        delete_originals=delete_originals
    )

if __name__ == "__main__":
    main()
