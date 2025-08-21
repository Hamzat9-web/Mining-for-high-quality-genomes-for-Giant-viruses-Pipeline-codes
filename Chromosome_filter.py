import os
import re
from Bio import SeqIO
import glob
import shutil

def extract_chromosomes_and_rename(input_fasta, output_dir):
    """Extract chromosome sequences and rename with species names."""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the accession from the parent directory name
    parent_dir = os.path.dirname(input_fasta)
    accession = os.path.basename(parent_dir)  # This will be like GCA_947044365.1 or GCF_947044365.1
    
    chromosome_records = []
    species_name = None
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Check if this is a chromosome sequence (look for various chromosome patterns)
        if re.search(r'chromosome\s*[:\s]\s*\d+|chromosome\s*[:\s]', record.description, re.IGNORECASE):
            
            # Extract species name from the description
            # Looking for patterns like "Petrosia ficiformis genome assembly" or "Whitmania pigra isolate"
            species_match = re.search(r'(\w+\s+\w+)\s+(?:genome\s+assembly|isolate)', record.description)
            if species_match:
                species_name = species_match.group(1).replace(' ', '_')
            
            chromosome_records.append(record)
    
    if chromosome_records:
        # Create output filename with species name
        if species_name:
            output_filename = f"{species_name}_{accession}_chromosomes.fna"
        else:
            output_filename = f"{accession}_chromosomes.fna"
        
        output_path = os.path.join(output_dir, output_filename)
        
        # Write chromosome sequences to new file
        SeqIO.write(chromosome_records, output_path, "fasta")
        print(f"✓ {species_name or accession}: {len(chromosome_records)} chromosome sequences")
        
        # Delete the entire directory after successful processing
        try:
            shutil.rmtree(parent_dir)
            print(f"  → Deleted entire directory: {accession}")
        except OSError as e:
            print(f"  ⚠ Warning: Could not delete directory {parent_dir}: {e}")
        
        return output_path, species_name or accession, len(chromosome_records)
    else:
        print(f"✗ {accession}: No chromosome sequences found")
        # Don't delete the directory if no chromosomes were found
        return None, accession, 0

def process_all_porifera_genomes(data_dir, output_dir):
    """Process all Porifera genome files (both GCA and GCF)."""
    
    # Find all GCA and GCF directories
    gca_dirs = glob.glob(os.path.join(data_dir, "GCA_*"))
    gcf_dirs = glob.glob(os.path.join(data_dir, "GCF_*"))
    
    # Combine both lists
    all_dirs = gca_dirs + gcf_dirs
    
    print(f"Found {len(gca_dirs)} GCA directories")
    print(f"Found {len(gcf_dirs)} GCF directories")
    print(f"Total: {len(all_dirs)} directories to process")
    print("=" * 60)
    
    processed_files = []
    summary_stats = []
    
    for genome_dir in sorted(all_dirs):
        # Find the genomic.fna file in this directory
        genome_files = glob.glob(os.path.join(genome_dir, "*_genomic.fna"))
        
        if genome_files:
            genome_file = genome_files[0]  # Should only be one
            result = extract_chromosomes_and_rename(genome_file, output_dir)
            
            if result[0]:  # If file was created
                processed_files.append(result[0])
            
            summary_stats.append({
                'accession': os.path.basename(genome_dir),
                'species': result[1],
                'chromosomes': result[2],
                'processed': result[0] is not None,
                'type': 'GCA' if genome_dir.startswith(os.path.join(data_dir, 'GCA_')) else 'GCF'
            })
        else:
            print(f"✗ {os.path.basename(genome_dir)}: No genomic.fna file found")
    
    # Print summary
    print("\n" + "=" * 60)
    print("PROCESSING SUMMARY")
    print("=" * 60)
    
    successful = len([s for s in summary_stats if s['processed']])
    total = len(summary_stats)
    gca_processed = len([s for s in summary_stats if s['processed'] and s['type'] == 'GCA'])
    gcf_processed = len([s for s in summary_stats if s['processed'] and s['type'] == 'GCF'])
    
    print(f"Successfully processed: {successful}/{total} genomes")
    print(f"  - GCA genomes: {gca_processed}")
    print(f"  - GCF genomes: {gcf_processed}")
    print(f"Total chromosome sequences extracted: {sum(s['chromosomes'] for s in summary_stats)}")
    
    print("\nProcessed genomes:")
    for stat in summary_stats:
        if stat['processed']:
            print(f"  ✓ {stat['species']} ({stat['accession']}, {stat['type']}): {stat['chromosomes']} chromosomes")
    
    if any(not s['processed'] for s in summary_stats):
        print("\nSkipped (no chromosomes found):")
        for stat in summary_stats:
            if not stat['processed']:
                print(f"  ✗ {stat['species']} ({stat['accession']}, {stat['type']})")
    
    return processed_files

# Run the processing
if __name__ == "__main__":
    # Set your paths
    data_directory = "/data/home/bt24090/viralrecall/Ctenophora/ncbi_dataset/data"
    output_directory = "/data/home/bt24090/viralrecall/Ctenophora/chromosomes_only"
    
    # Process all genomes
    processed_files = process_all_porifera_genomes(data_directory, output_directory)
    
    print(f"\n🎉 All done! Check your results in: {output_directory}")
