import os
import subprocess
import argparse
from Bio import SeqIO
import glob
import pandas as pd

def run_blast_pass1(input_file, output_file, db_path):
    """Run BLAST+ Pass 1: Get top hit only to identify viral candidates"""
    cmd = [
        "blastp",
        "-query", input_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6",
        "-max_target_seqs", "1",  # Only top hit
        "-evalue", "1e-5"
    ]
    subprocess.run(cmd, check=True)

def run_blast_pass2(input_file, output_file, db_path):
    """Run BLAST+ Pass 2: Get all significant hits for validation"""
    cmd = [
        "blastp",
        "-query", input_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "6",
        "-max_target_seqs", "1000",  # Large number to get all significant hits
        "-evalue", "1e-5"
    ]
    subprocess.run(cmd, check=True)

def identify_viral_candidates(blast_tsv):
    """Pass 1: Identify sequences with viral top hits"""
    viral_candidates = set()
    
    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        return viral_candidates
    
    try:
        blast_df = pd.read_csv(blast_tsv, sep='\t', header=None,
                              names=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
    except Exception as e:
        print(f"    Error reading Pass 1 BLAST results: {e}")
        return viral_candidates
    
    for _, row in blast_df.iterrows():
        if row['sseqid'].startswith('GV_'):
            viral_candidates.add(row['qseqid'])
    
    return viral_candidates

def validate_viral_candidates(blast_tsv, viral_candidates):
    """Pass 2: Validate viral candidates using relative scoring"""
    sequences_to_keep = {}
    
    # Statistics
    pass1_candidates = len(viral_candidates)
    validated_sequences = 0
    failed_euk_competition = 0
    no_viral_hit = 0
    
    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        print(f"    Pass 2: No BLAST hits found")
        return sequences_to_keep, pass1_candidates, validated_sequences, failed_euk_competition, no_viral_hit
    
    try:
        blast_df = pd.read_csv(blast_tsv, sep='\t', header=None,
                              names=['qseqid', 'sseqid', 'pident', 'length',
                                    'mismatch', 'gapopen', 'qstart', 'qend',
                                    'sstart', 'send', 'evalue', 'bitscore'])
    except Exception as e:
        print(f"    Error reading Pass 2 BLAST results: {e}")
        return sequences_to_keep, pass1_candidates, validated_sequences, failed_euk_competition, no_viral_hit
    
    # Group hits by query sequence
    for query_id in viral_candidates:
        query_hits = blast_df[blast_df['qseqid'] == query_id]
        
        if len(query_hits) == 0:
            no_viral_hit += 1
            continue
        
        # Find best viral and eukaryotic hits
        viral_hits = query_hits[query_hits['sseqid'].str.startswith('GV_')]
        euk_hits = query_hits[
            query_hits['sseqid'].str.startswith('EUK_') |
            (query_hits['sseqid'].str.contains(r'\[', regex=True) & query_hits['sseqid'].str.contains(r'\]', regex=True))
        ]
        
        if len(viral_hits) == 0:
            no_viral_hit += 1
            continue
        
        best_viral_hit = viral_hits.iloc[0]  # Already sorted by e-value
        best_viral_evalue = best_viral_hit['evalue']
        
        # Check if there are competing eukaryotic hits
        if len(euk_hits) > 0:
            best_euk_evalue = euk_hits.iloc[0]['evalue']
            
            # Relative scoring: viral hit must be 1000x better (2 orders of magnitude)
            if best_euk_evalue <= 1e-5:  # Eukaryotic hit is significant
                if best_viral_evalue > (best_euk_evalue / 1000):  # Viral not 100x better
                    failed_euk_competition += 1
                    continue
        
        # Sequence passes validation
        validated_sequences += 1
        virus_id = best_viral_hit['sseqid'].replace('GV_', '', 1)
        sequences_to_keep[query_id] = {
            'virus_id': virus_id,
            'viral_evalue': best_viral_evalue,
            'best_euk_evalue': euk_hits.iloc[0]['evalue'] if len(euk_hits) > 0 else None
        }
    
    return sequences_to_keep, pass1_candidates, validated_sequences, failed_euk_competition, no_viral_hit

def create_viral_candidate_fasta(input_fasta, viral_candidates, temp_fasta):
    """Create temporary FASTA file with only viral candidates for Pass 2"""
    sequences = SeqIO.parse(input_fasta, "fasta-blast")
    viral_seqs = []
    
    for seq in sequences:
        if seq.id in viral_candidates:
            viral_seqs.append(seq)
    
    SeqIO.write(viral_seqs, temp_fasta, "fasta")
    return len(viral_seqs)

def process_blast_results_two_pass(pass1_tsv, pass2_tsv, input_fasta, output_fasta, temp_fasta):
    """Two-pass BLAST processing with validation"""
    
    # Pass 1: Identify viral candidates
    print(f"    Pass 1: Identifying viral candidates...")
    viral_candidates = identify_viral_candidates(pass1_tsv)
    
    if not viral_candidates:
        print(f"    Pass 1: No viral candidates found")
        with open(output_fasta, 'w') as f:
            pass
        return 0
    
    print(f"    Pass 1: Found {len(viral_candidates)} viral candidates")
    
    # Create temporary FASTA with viral candidates
    temp_count = create_viral_candidate_fasta(input_fasta, viral_candidates, temp_fasta)
    print(f"    Created temporary FASTA with {temp_count} sequences")
    
    # Pass 2: Validate candidates with relative scoring
    print(f"    Pass 2: Validating candidates with relative scoring...")
    sequences_to_keep, pass1_count, validated_count, failed_euk, no_viral = \
        validate_viral_candidates(pass2_tsv, viral_candidates)
    
    # Print detailed statistics
    print(f"    Validation Statistics:")
    print(f"      Pass 1 candidates: {pass1_count}")
    print(f"      Validated sequences: {validated_count}")
    print(f"      Failed eukaryotic competition: {failed_euk}")
    print(f"      Lost viral hit in Pass 2: {no_viral}")
    
    # Write final filtered sequences
    sequences = SeqIO.parse(input_fasta, "fasta-blast")
    filtered_sequences = []
    
    for seq in sequences:
        if seq.id in sequences_to_keep:
            seq_info = sequences_to_keep[seq.id]
            virus_id = seq_info['virus_id']
            viral_eval = seq_info['viral_evalue']
            euk_eval = seq_info['best_euk_evalue']
            
            # Enhanced description with validation info
            if euk_eval is not None:
                seq.description = f"{seq.description} {virus_id} [V:{viral_eval:.2e},E:{euk_eval:.2e}]"
            else:
                seq.description = f"{seq.description} {virus_id} [V:{viral_eval:.2e}]"
            
            filtered_sequences.append(seq)
    
    SeqIO.write(filtered_sequences, output_fasta, "fasta")
    return len(filtered_sequences)

def process_species_directory(species_dir, db_path):
    """Process a single species directory with two-pass validation"""
    species_name = os.path.basename(species_dir)
    print(f"\n{'='*60}")
    print(f"Processing species: {species_name}")
    print(f"{'='*60}")
    
    # Define paths
    query_sequences_dir = os.path.join(species_dir, "query_sequences")
    filter_results_dir = os.path.join(species_dir, "blast_filter_results")
    blast_results_dir = os.path.join(species_dir, "blast_results")
    temp_dir = os.path.join(species_dir, "temp")
    
    # Check if query_sequences directory exists
    if not os.path.exists(query_sequences_dir):
        print(f"  Warning: query_sequences directory not found in {species_dir}")
        return []
    
    # Create directories
    os.makedirs(filter_results_dir, exist_ok=True)
    os.makedirs(blast_results_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Find FASTA files
    fasta_pattern = os.path.join(query_sequences_dir, "*.fasta")
    input_files = glob.glob(fasta_pattern)
    
    if not input_files:
        print(f"  No .fasta files found in {query_sequences_dir}")
        return []
    
    print(f"  Found {len(input_files)} FASTA files to process")
    
    species_summary = []
    
    for input_file in input_files:
        filename = os.path.basename(input_file)
        base_name = os.path.splitext(filename)[0]
        
        # Output files
        pass1_output = os.path.join(blast_results_dir, f"{base_name}_pass1.tsv")
        pass2_output = os.path.join(blast_results_dir, f"{base_name}_pass2.tsv")
        temp_fasta = os.path.join(temp_dir, f"{base_name}_candidates.fasta")
        filtered_output = os.path.join(filter_results_dir, f"filtered_{filename}")
        
        print(f"\n  Processing {filename}...")
        
        # Count input sequences
        try:
            input_count = sum(1 for _ in SeqIO.parse(input_file, "fasta-blast"))
        except Exception as e:
            print(f"    Error counting sequences: {e}")
            input_count = 0
        
        if input_count == 0:
            print(f"    Warning: {filename} contains no sequences. Skipping.")
            species_summary.append((filename, 0, 0))
            continue
        
        print(f"    Input sequences: {input_count}")
        
        # Pass 1: Run BLAST to identify viral candidates
        try:
            run_blast_pass1(input_file, pass1_output, db_path)
        except subprocess.CalledProcessError as e:
            print(f"    Error in BLAST Pass 1: {e}")
            species_summary.append((filename, input_count, 0))
            continue
        
        # Pass 2: Run BLAST on viral candidates only
        viral_candidates = identify_viral_candidates(pass1_output)
        if viral_candidates:
            create_viral_candidate_fasta(input_file, viral_candidates, temp_fasta)
            try:
                run_blast_pass2(temp_fasta, pass2_output, db_path)
            except subprocess.CalledProcessError as e:
                print(f"    Error in BLAST Pass 2: {e}")
                species_summary.append((filename, input_count, 0))
                continue
        else:
            # Create empty Pass 2 file
            with open(pass2_output, 'w') as f:
                pass
        
        # Process results with two-pass validation
        try:
            kept_count = process_blast_results_two_pass(
                pass1_tsv=pass1_output,
                pass2_tsv=pass2_output,
                input_fasta=input_file,
                output_fasta=filtered_output,
                temp_fasta=temp_fasta
            )
            print(f"    Final result: {kept_count} out of {input_count} sequences")
        except Exception as e:
            print(f"    Error processing results: {e}")
            kept_count = 0
        
        # Clean up temp file
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        
        species_summary.append((filename, input_count, kept_count))
    
    # Clean up temp directory if empty
    try:
        os.rmdir(temp_dir)
    except OSError:
        pass  # Directory not empty or other error
    
    # Print species summary
    print(f"\n  Species Summary for {species_name}:")
    print(f"  {'-'*50}")
    print(f"  {'File':<25} {'Input':<10} {'Kept':<10}")
    print(f"  {'-'*50}")
    for filename, input_count, kept_count in species_summary:
        print(f"  {filename:<25} {input_count:<10} {kept_count:<10}")
    
    # Calculate species totals
    total_input = sum(i for _, i, _ in species_summary)
    total_kept = sum(k for _, _, k in species_summary)
    kept_percentage = (total_kept / total_input * 100) if total_input > 0 else 0
    
    print(f"  {'-'*50}")
    print(f"  {'SPECIES TOTAL':<25} {total_input:<10} {total_kept:<10} ({kept_percentage:.1f}%)")
    
    return species_summary

def main():
    parser = argparse.ArgumentParser(description='Filter FASTA files for giant virus sequences using two-pass BLAST validation')
    parser.add_argument('--screening_results_dir', required=True,
                        help='Path to screening_results directory containing species subdirectories')
    parser.add_argument('--db_path', required=True, help='Path to BLAST database')
    
    args = parser.parse_args()
    
    # Check if screening_results directory exists
    if not os.path.exists(args.screening_results_dir):
        print(f"Error: screening_results directory not found: {args.screening_results_dir}")
        return
    
    # Find all species directories
    species_dirs = []
    for item in os.listdir(args.screening_results_dir):
        item_path = os.path.join(args.screening_results_dir, item)
        if os.path.isdir(item_path):
            query_seq_path = os.path.join(item_path, "query_sequences")
            if os.path.exists(query_seq_path):
                species_dirs.append(item_path)
    
    if not species_dirs:
        print(f"No species directories with query_sequences subdirectories found in {args.screening_results_dir}")
        return
    
    print(f"Found {len(species_dirs)} species directories to process")
    print(f"Using two-pass BLAST validation with relative scoring")
    print(f"Pass 1: Identify viral candidates (top hit only)")
    print(f"Pass 2: Validate with relative scoring (viral must be 100x better than eukaryotic)")
    
    # Process each species directory
    global_summary = []
    
    for species_dir in sorted(species_dirs):
        species_summary = process_species_directory(species_dir, args.db_path)
        
        species_name = os.path.basename(species_dir)
        for filename, input_count, kept_count in species_summary:
            global_summary.append((species_name, filename, input_count, kept_count))
    
    # Print global summary
    print(f"\n\n{'='*80}")
    print(f"GLOBAL SUMMARY - TWO-PASS BLAST VALIDATION")
    print(f"{'='*80}")
    print(f"{'Species':<30} {'File':<25} {'Input':<10} {'Kept':<10}")
    print(f"{'-'*80}")
    
    for species_name, filename, input_count, kept_count in global_summary:
        print(f"{species_name:<30} {filename:<25} {input_count:<10} {kept_count:<10}")
    
    # Calculate global totals
    total_input = sum(i for _, _, i, _ in global_summary)
    total_kept = sum(k for _, _, _, k in global_summary)
    kept_percentage = (total_kept / total_input * 100) if total_input > 0 else 0
    
    print(f"{'-'*80}")
    print(f"{'GLOBAL TOTAL':<30} {'':<25} {total_input:<10} {total_kept:<10} ({kept_percentage:.1f}%)")
    
    print(f"\nTwo-pass BLAST processing complete!")
    print(f"Processed {len(species_dirs)} species directories with maximum false positive reduction.")

if __name__ == "__main__":
    main()
