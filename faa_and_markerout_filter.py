import os
import glob

def read_markerout_targets(markerout_file):
    """Read target names from markerout file, skipping all header lines."""
    targets = set()
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
                if parts and len(parts) > 1:  # Ensure it's a data line with multiple columns
                    target = parts[0]
                    # Additional validation: target should look like a sequence ID
                    if '.' in target and '_' in target:  # Basic format check
                        targets.add(target)
                    
    except Exception as e:
        print(f"Error reading {markerout_file}: {e}")
        return set()
    
    return targets

def extract_sequences(faa_file, target_names, output_file):
    """Extract sequences matching target names from faa file."""
    sequences_found = 0
    total_sequences = 0
    
    try:
        with open(faa_file, 'r') as f_in, open(output_file, 'w') as f_out:
            write_sequence = False
            
            for line in f_in:
                if line.startswith('>'):
                    total_sequences += 1
                    
                    # Reset write flag
                    write_sequence = False
                    
                    # Check if any target name appears at the start of header (after >)
                    for target in target_names:
                        if line.strip().startswith(f'>{target} ') or line.strip() == f'>{target}':
                            write_sequence = True
                            sequences_found += 1
                            break
                    
                    if write_sequence:
                        f_out.write(line)
                
                elif write_sequence:
                    f_out.write(line)
                    
    except Exception as e:
        print(f"Error processing {faa_file}: {e}")
        return 0, 0
    
    return sequences_found, total_sequences

def process_directory(dir_path):
    """Process a single directory containing .markerout and .faa files."""
    print(f"\nProcessing directory: {os.path.basename(dir_path)}")
    
    # Find .markerout and .faa files (exclude already filtered files)
    markerout_files = glob.glob(os.path.join(dir_path, "*.markerout"))
    faa_files = [f for f in glob.glob(os.path.join(dir_path, "*.faa")) 
                 if not f.endswith("_filtered_sequences.faa")]
    
    if not markerout_files:
        print(f"  No .markerout files found")
        return False
    
    if not faa_files:
        print(f"  No .faa files found")
        return False
    
    success = False
    
    # Process each combination of markerout and faa files
    for markerout_file in markerout_files:
        for faa_file in faa_files:
            print(f"  Processing: {os.path.basename(markerout_file)} + {os.path.basename(faa_file)}")
            
            # Read target names
            target_names = read_markerout_targets(markerout_file)
            if not target_names:
                print(f"    No valid target names found")
                continue
            
            print(f"    Found {len(target_names)} target names")
            
            # Count total sequences in FAA
            with open(faa_file, 'r') as f:
                total_faa_sequences = sum(1 for line in f if line.startswith('>'))
            
            # Create output filename
            base_name = os.path.basename(markerout_file).replace('.markerout', '')
            output_file = os.path.join(dir_path, f"{base_name}_filtered_sequences.faa")
            
            # Extract sequences
            sequences_found, sequences_processed = extract_sequences(faa_file, target_names, output_file)
            
            print(f"    Total sequences in FAA: {total_faa_sequences}")
            print(f"    Sequences extracted: {sequences_found}")
            
            if sequences_found > 0:
                filtering_ratio = (sequences_found / total_faa_sequences) * 100 if total_faa_sequences > 0 else 0
                print(f"    Filtering ratio: {sequences_found}/{total_faa_sequences} = {filtering_ratio:.1f}%")
                
                if sequences_found == total_faa_sequences:
                    print(f"    ⚠️  WARNING: All sequences matched - check if filtering is working correctly")
                else:
                    print(f"    ✅ Success: Filtering working correctly!")
                    print(f"    Output saved to: {os.path.basename(output_file)}")
                success = True
            else:
                print(f"    ❌ No matching sequences found")
                # Remove empty output file
                if os.path.exists(output_file):
                    os.remove(output_file)
    
    return success

def main():
    # Base directory path
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
    print(f"BATCH PROCESSING COMPLETE!")
    print(f"Successfully processed: {processed_count} directories")
    print(f"Failed to process: {failed_count} directories")
    print(f"Total directories: {len(subdirectories)}")
    
    if processed_count > 0:
        print(f"\nFiltered sequence files have been created with suffix '_filtered_sequences.faa'")

if __name__ == "__main__":
    main()
