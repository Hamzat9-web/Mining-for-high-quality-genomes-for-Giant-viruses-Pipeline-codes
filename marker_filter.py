import os
import glob
from collections import defaultdict

def read_markerout_mapping(markerout_file):
    """Read target name to query name mapping from markerout file."""
    target_to_query = {}
    
    try:
        with open(markerout_file, 'r') as f:
            lines = f.readlines()
            
        # Skip all header/comment lines and process only data lines
        for line in lines:
            line = line.strip()
            # Skip empty lines, comment lines, and header separator lines
            if (line and 
                not line.startswith('#') and 
                not line.startswith('-') and
                not line.startswith('=')):
                
                parts = line.split()
                if parts and len(parts) >= 3:  # Ensure we have at least target, accession, query
                    target_name = parts[0]  # 1st column: target name
                    query_name = parts[2]   # 3rd column: query name
                    
                    # Basic validation
                    if '.' in target_name and '_' in target_name:
                        target_to_query[target_name] = query_name
                    
    except Exception as e:
        print(f"Error reading {markerout_file}: {e}")
        return {}
    
    return target_to_query

def group_sequences_by_query(filtered_faa_file, target_to_query, output_dir, base_name):
    """Group sequences by query name and write separate files."""
    
    # Create query_sequences directory if it doesn't exist
    query_sequences_dir = os.path.join(output_dir, "query_sequences")
    try:
        os.makedirs(query_sequences_dir, exist_ok=True)
        print(f"    Created/verified directory: {os.path.basename(query_sequences_dir)}")
    except Exception as e:
        print(f"    Error creating directory {query_sequences_dir}: {e}")
        return {}
    
    # Dictionary to store sequences for each query
    query_sequences = defaultdict(list)
    
    try:
        with open(filtered_faa_file, 'r') as f:
            current_header = ""
            current_sequence = []
            
            for line in f:
                if line.startswith('>'):
                    # Process previous sequence if exists
                    if current_header and current_sequence:
                        process_sequence(current_header, current_sequence, target_to_query, query_sequences)
                    
                    # Start new sequence
                    current_header = line.strip()
                    current_sequence = []
                else:
                    current_sequence.append(line.strip())
            
            # Process last sequence
            if current_header and current_sequence:
                process_sequence(current_header, current_sequence, target_to_query, query_sequences)
                
    except Exception as e:
        print(f"Error reading {filtered_faa_file}: {e}")
        return {}
    
    # Write separate files for each query in the query_sequences directory
    query_counts = {}
    for query_name, sequences in query_sequences.items():
        if sequences:  # Only create files for queries that have sequences
            output_file = os.path.join(query_sequences_dir, f"{base_name}_{query_name}.fasta")
            
            try:
                with open(output_file, 'w') as f:
                    for header, seq_lines in sequences:
                        f.write(f"{header}\n")
                        f.write(f"{''.join(seq_lines)}\n")
                
                query_counts[query_name] = len(sequences)
                print(f"    Created query_sequences/{os.path.basename(output_file)} with {len(sequences)} sequences")
                
            except Exception as e:
                print(f"    Error writing {output_file}: {e}")
    
    return query_counts

def process_sequence(header, sequence_lines, target_to_query, query_sequences):
    """Process a single sequence and assign it to the correct query group."""
    
    # Extract target name from header (first part after >)
    header_parts = header[1:].split()  # Remove > and split
    if header_parts:
        target_name = header_parts[0]
        
        # Find the query name for this target
        if target_name in target_to_query:
            query_name = target_to_query[target_name]
            query_sequences[query_name].append((header, sequence_lines))
        else:
            print(f"    Warning: Target '{target_name}' not found in markerout mapping")

def process_directory(dir_path):
    """Process a single directory to group sequences by query name."""
    
    dir_name = os.path.basename(dir_path)
    print(f"\nProcessing directory: {dir_name}")
    
    # Find markerout and filtered faa files
    markerout_files = glob.glob(os.path.join(dir_path, "*.markerout"))
    filtered_faa_files = glob.glob(os.path.join(dir_path, "*_filtered_sequences.faa"))
    
    if not markerout_files:
        print(f"  No .markerout files found")
        return False
    
    if not filtered_faa_files:
        print(f"  No *_filtered_sequences.faa files found")
        return False
    
    success = False
    
    # Process each markerout file with corresponding filtered faa file
    for markerout_file in markerout_files:
        markerout_base = os.path.basename(markerout_file).replace('.markerout', '')
        
        # Find corresponding filtered faa file
        corresponding_faa = None
        for faa_file in filtered_faa_files:
            if markerout_base in os.path.basename(faa_file):
                corresponding_faa = faa_file
                break
        
        if not corresponding_faa:
            print(f"  No corresponding filtered FAA file found for {os.path.basename(markerout_file)}")
            continue
        
        print(f"  Processing: {os.path.basename(markerout_file)} + {os.path.basename(corresponding_faa)}")
        
        # Read target to query mapping
        target_to_query = read_markerout_mapping(markerout_file)
        if not target_to_query:
            print(f"    No valid mappings found in markerout file")
            continue
        
        print(f"    Found mappings for {len(target_to_query)} targets")
        
        # Get base name for output files (use species name part)
        # Extract species name from directory name
        species_name = dir_name.split('_GCA_')[0] if '_GCA_' in dir_name else dir_name
        
        # Group sequences by query
        query_counts = group_sequences_by_query(corresponding_faa, target_to_query, dir_path, species_name)
        
        if query_counts:
            print(f"    Created {len(query_counts)} query-specific files:")
            for query, count in sorted(query_counts.items()):
                print(f"      {query}: {count} sequences")
            success = True
        else:
            print(f"     No sequences were grouped")
    
    return success

def main():
    """Main function to process all directories."""
    
    base_dir = ""
    
    if not os.path.exists(base_dir):
        print(f"Error: Directory {base_dir} does not exist!")
        return
    
    print(f"Scanning directory: {base_dir}")
    
    # Find all subdirectories
    subdirectories = [d for d in os.listdir(base_dir) 
                     if os.path.isdir(os.path.join(base_dir, d))]
    
    if not subdirectories:
        print("No subdirectories found!")
        return
    
    print(f"Found {len(subdirectories)} subdirectories to process")
    
    processed_count = 0
    failed_count = 0
    
    # Process each subdirectory
    for subdir in sorted(subdirectories):
        subdir_path = os.path.join(base_dir, subdir)
        try:
            if process_directory(subdir_path):
                processed_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"  Error processing {subdir}: {e}")
            failed_count += 1
    
    print(f"\n" + "="*60)
    print(f"GROUPING COMPLETE!")
    print(f"Successfully processed: {processed_count} directories")
    print(f"Failed to process: {failed_count} directories")
    print(f"Total directories: {len(subdirectories)}")
    
    if processed_count > 0:
        print(f"\nQuery-specific FASTA files created in 'query_sequences' subdirectories")
        print(f"File naming pattern: [SpeciesName]_[QueryName].fasta")

if __name__ == "__main__":
    main()
