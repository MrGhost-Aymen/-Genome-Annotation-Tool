import os
import json
from io import StringIO
from flask import Flask, render_template, request, jsonify, send_file, redirect, url_for, send_from_directory, Response, session
from werkzeug.utils import secure_filename
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from pathlib import Path
import uuid
import gzip
import shutil
import struct
import io
import time
import subprocess
import traceback
import mimetypes

print("DEBUG: App script started successfully.") # Added for early confirmation

# --- CONFIGURATION: SET PATHS TO EXTERNAL EXECUTABLES ---
# Make sure these paths are correct for your system!
BGZIP_PATH = r"C:\msys64\mingw64\bin\bgzip.exe" # <-- Your bgzip.exe path
TABIX_PATH = r"C:\msys64\mingw64\bin\tabix.exe" # <-- Your tabix.exe path
SORT_PATH = r"C:\msys64\usr\bin\sort.exe"    # <-- Your sort.exe path
# -----------------------------------------------------------------------

# --- FLAG: SET TO True TO BYPASS BGZIP/TABIX AND SERVE UNCOMPRESSED GFF3 ---
# Set to True for the current working solution that avoids "invalid block type"
USE_UNCOMPRESSED_GFF = True 
# ----------------------------------------------------------------------------

mimetypes.add_type('application/x-gzip', '.gz')
mimetypes.add_type('application/octet-stream', '.tbi')

app = Flask(__name__, static_folder='static', static_url_path='/static')
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
app.secret_key = os.urandom(24) # THIS IS CRUCIAL FOR SESSIONS!
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Helper for escaping GFF3 values (Moved to higher scope)
def _escape_gff3_value(value):
    """Escapes special characters in a GFF3 attribute value."""
    return str(value).replace('%', '%25') \
                     .replace(';', '%3B') \
                     .replace('=', '%3D') \
                     .replace(',', '%2C') \
                     .replace('&', '%26') \
                     .replace('\t', '%09') \
                     .replace('\n', '%0A') \
                     .replace('\r', '%0D')

# Helper for unescaping GFF3 values (Moved to higher scope)
def _unescape_gff3_value(value):
    """Unescapes special characters in a GFF3 attribute value."""
    return str(value).replace('%0D', '\r') \
                     .replace('%0A', '\n') \
                     .replace('%09', '\t') \
                     .replace('%26', '&') \
                     .replace('%2C', ',') \
                     .replace('%3D', '=') \
                     .replace('%3B', ';') \
                     .replace('%25', '%')


def load_codon_tables():
    try:
        # Load from static/codon_tables.json
        with open(os.path.join(app.static_folder, 'codon_tables.json'), 'r', encoding='utf-8') as f:
            raw_tables = json.load(f)
            codon_tables = []
            for table_key, table_data in raw_tables.items():
                codon_tables.append({
                    'id': str(table_data.get('id', table_key)),
                    'name': table_data['name'],
                    'start_codons': table_data.get('start', []),
                    'stop_codons': [
                        codon for codon, aa in table_data['table'].items()
                        if aa == '*'
                    ],
                    'table': table_data['table']
                })
            return codon_tables
    except Exception as e:
        print(f"ERROR: Error loading codon tables: {str(e)}")
        # Provide a minimal default if loading fails to prevent crash
        return [{'id': '1', 'name': 'Standard', 'start_codons': ['ATG'], 'stop_codons': ['TAA', 'TAG', 'TGA'], 'table': {}}]

CODON_TABLES = load_codon_tables()

@app.route('/api/codon_tables')
def get_codon_tables():
    return jsonify(CODON_TABLES)

def setup_igv_data_and_config(seq_record, gff_feature_lines_content, session_id):
    print(f"DEBUG: Setting up IGV data for session {session_id}")
    igv_data_base_dir = Path("static/igvjs/data")
    session_data_dir = igv_data_base_dir / session_id
    session_data_dir.mkdir(parents=True, exist_ok=True)
    session_data_dir_abs = session_data_dir.resolve()
    print(f"DEBUG: Session data directory created at {session_data_dir_abs}")

    fasta_data_filename = "sequence.fasta"
    fasta_filepath_fs = session_data_dir_abs / fasta_data_filename
    try:
        with open(fasta_filepath_fs, 'w') as f:
            SeqIO.write(seq_record, f, 'fasta')
        print(f"DEBUG: FASTA file written to {fasta_filepath_fs}")
    except Exception as e:
        print(f"ERROR: Failed to write FASTA file: {e}")
        raise

    fai_filepath_fs = fasta_filepath_fs.with_suffix('.fasta.fai')
    try:
        with open(fai_filepath_fs, "w") as f:
            seq_len = len(seq_record.seq)
            line_length = 80 # Standard FASTA line length
            bytes_per_line = line_length + 1 # +1 for newline character
            # Calculate offset to the first base
            offset = len(f">{seq_record.id}\n") 
            
            f.write(f"{seq_record.id}\t{seq_len}\t{offset}\t{line_length}\t{bytes_per_line}\n")
        print(f"DEBUG: FASTA .fai file written to {fai_filepath_fs}")
    except Exception as e:
        print(f"ERROR: Failed to write FASTA .fai file: {e}")
        raise

    gff_data_filename_base = "annotations.gff3"
    
    header_lines = [
        "##gff-version 3",
        f"##sequence-region {seq_record.id} 1 {len(seq_record.seq)}"
    ]
    
    # Write the GFF data for initial setup
    _write_gff_data_to_file_system(
        session_data_dir_abs, 
        gff_data_filename_base, 
        header_lines, 
        gff_feature_lines_content # This content is already formatted as lines
    )
    print(f"DEBUG: GFF data written to file system for IGV.js")

    time.sleep(0.5) # Small delay to ensure files are written, especially for FASTA.fai

    # --- Generate IGV.js Configuration ---
    igv_config = {
        "genome": {
            "id": seq_record.id,
            "fastaURL": f"/static/igvjs/data/{session_id}/{fasta_data_filename}",
            "indexURL": f"/static/igvjs/data/{session_id}/{fasta_data_filename}.fai"
        },
        "tracks": [
            {
                "name": "Genome Annotations",
                "format": "gff3",
                "url": f"/static/igvjs/data/{session_id}/{gff_data_filename_base}.gz" if not USE_UNCOMPRESSED_GFF else f"/static/igvjs/data/{session_id}/{gff_data_filename_base}",
                "indexURL": f"/static/igvjs/data/{session_id}/{gff_data_filename_base}.gz.tbi" if not USE_UNCOMPRESSED_GFF else None,
                "displayMode": "EXPANDED",
                "visibilityWindow": 1000000,
                "height": 200
            }
        ],
        "showNavigation": True,
        "showRuler": True,
        "minimumBases": 40,
        "genomeList": []
    }
    
    initial_locus_string = f"{seq_record.id}:1-{min(len(seq_record.seq), 5000)}" 
    print(f"DEBUG: IGV.js configuration generated.")
    
    return {
        'session_id': session_id,
        'sequence_id': seq_record.id,
        'sequence_length': len(seq_record.seq),
        'features': [feature_to_dict(f) for f in seq_record.features], 
        'igv_config': igv_config,
        'initial_locus_string': initial_locus_string
    }

def _write_gff_data_to_file_system(session_data_dir_abs, gff_data_filename_base, header_lines, feature_gff_content):
    """
    Helper function to write GFF data to the file system, handling sorting and compression.
    `feature_gff_content` is expected to be a single string containing all feature lines.
    """
    print(f"DEBUG: _write_gff_data_to_file_system started for {gff_data_filename_base}")
    gff_gz_filepath_fs = session_data_dir_abs / (gff_data_filename_base + ".gz")
    gff_tbi_filepath_fs = session_data_dir_abs / (gff_data_filename_base + ".gz.tbi")
    gff_uncompressed_filepath_fs = session_data_dir_abs / gff_data_filename_base

    temp_feature_file = session_data_dir_abs / "temp_features_unsorted.gff3"
    temp_gff_sorted_filepath = session_data_dir_abs / "temp_features_sorted.gff3"
    
    try:
        with open(temp_feature_file, 'w', encoding='utf-8', newline='\n') as f:
            f.write(feature_gff_content)
        print(f"DEBUG: Feature lines written to temporary file: {temp_feature_file}")
    except Exception as e:
        print(f"ERROR: Failed to write temporary feature file: {e}")
        raise
    
    try:
        env = os.environ.copy()
        env['LC_ALL'] = 'C' # Use C locale for consistent sort order

        sort_command = [
            SORT_PATH,
            '-k1,1',  # Sort by sequence ID (column 1)
            '-k4,4n', # Then by start position (column 4, numeric)
            str(temp_feature_file)
        ]
        
        print(f"DEBUG: Running sort command: {' '.join(sort_command)}")
        
        process = subprocess.run(
            sort_command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            encoding='utf-8', # Specify encoding for subprocess output
            errors='replace' # Handle decoding errors
        )
        
        if process.stderr:
            print(f"SORT STDERR: {process.stderr.strip()}")

        with open(temp_gff_sorted_filepath, 'w', encoding='utf-8', newline='\n') as outfile:
            outfile.write(process.stdout)
        
        print(f"DEBUG: Sorted features into {temp_gff_sorted_filepath}")
        
        with open(gff_uncompressed_filepath_fs, 'w', encoding='utf-8', newline='\n') as outfile:
            outfile.write('\n'.join(header_lines) + '\n')
            with open(temp_gff_sorted_filepath, 'r', encoding='utf-8', newline='\n') as f_sorted:
                outfile.write(f_sorted.read())
                
        print(f"DEBUG: Combined headers and sorted features into {gff_uncompressed_filepath_fs}")
        
    except FileNotFoundError:
        print(f"ERROR: sort executable not found at {SORT_PATH}. Please check the path configuration.")
        raise RuntimeError(f"sort executable not found at {SORT_PATH}. Please check the path configuration.")
    except subprocess.CalledProcessError as e:
        stdout_output = e.stdout.strip() if e.stdout else ""
        stderr_output = e.stderr.strip() if e.stderr else ""
        print(f"ERROR: sort failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
        raise RuntimeError(f"sort failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
    except Exception as e:
        print(f"ERROR: An unexpected error occurred during GFF3 sorting: {e}")
        traceback.print_exc()
        raise RuntimeError(f"An unexpected error occurred during GFF3 sorting: {e}")
    finally:
        for temp_file in [temp_feature_file, temp_gff_sorted_filepath]:
            try:
                if temp_file.exists():
                    temp_file.unlink()
                    print(f"DEBUG: Removed temporary file: {temp_file}")
            except Exception as e:
                print(f"WARNING: Could not remove {temp_file}: {e}")

    if USE_UNCOMPRESSED_GFF:
        print(f"DEBUG: Using uncompressed GFF: {gff_uncompressed_filepath_fs}")
        if gff_gz_filepath_fs.exists():
            try:
                gff_gz_filepath_fs.unlink()
                print(f"DEBUG: Removed lingering .gz file: {gff_gz_filepath_fs}")
            except Exception as e:
                print(f"WARNING: Could not remove lingering .gz file: {e}")
        if gff_tbi_filepath_fs.exists():
            try:
                gff_tbi_filepath_fs.unlink()
                print(f"DEBUG: Removed lingering .tbi file: {gff_tbi_filepath_fs}")
            except Exception as e:
                print(f"WARNING: Could not remove lingering .tbi file: {e}")
    else:
        # BGZIP and TABIX steps (currently bypassed by USE_UNCOMPRESSED_GFF = True)
        try:
            abs_gff_uncompressed_filepath = str(gff_uncompressed_filepath_fs)
            abs_gff_gz_filepath_fs = str(gff_gz_filepath_fs)

            print("DEBUG: Using external bgzip.exe for compression.")
            
            bgzip_command = [BGZIP_PATH, '-@', '2', '-c', abs_gff_uncompressed_filepath]

            print(f"DEBUG: Running bgzip command: {' '.join(bgzip_command)}")

            process = subprocess.run(
                bgzip_command,
                check=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            
            if process.stderr:
                print(f"BGZIP STDERR: {process.stderr.decode('utf-8', errors='replace').strip()}")

            with open(abs_gff_gz_filepath_fs, 'wb') as f_gz:
                f_gz.write(process.stdout)

            if not gff_gz_filepath_fs.exists() or gff_gz_filepath_fs.stat().st_size == 0:
                raise RuntimeError(f"bgzip failed: Output file {gff_gz_filepath_fs} was not created or is empty.")
            print(f"DEBUG: BGZIP compression successful to {gff_gz_filepath_fs}")

            try:
                with gzip.open(abs_gff_gz_filepath_fs, 'rt', encoding='utf-8') as test_f:
                    first_line = test_f.readline()
                    if not first_line.startswith('##gff-version 3'):
                        raise RuntimeError("Decompressed content doesn't start with GFF header")
                    
                print("DEBUG: BGZF validation passed (via Python's gzip)")
            except Exception as e:
                if os.path.exists(abs_gff_gz_filepath_fs):
                    os.remove(abs_gff_gz_filepath_fs)
                print(f"ERROR: BGZF validation failed (via Python's gzip): {str(e)}")
                raise RuntimeError(f"BGZF validation failed (via Python's gzip): {str(e)}")
                
        except FileNotFoundError:
            print(f"ERROR: bgzip executable not found at {BGZIP_PATH}. Please check the path configuration.")
            raise RuntimeError(f"bgzip executable not found at {BGZIP_PATH}. Please check the path configuration.")
        except subprocess.CalledProcessError as e:
            stdout_output = e.stdout.decode('utf-8', errors='replace').strip() if e.stdout else ""
            stderr_output = e.stderr.decode('utf-8', errors='replace').strip() if e.stderr else ""
            print(f"ERROR: bgzip failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
            raise RuntimeError(f"bgzip failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
        except Exception as e:
            print(f"ERROR: An unexpected error occurred during bgzip compression: {e}")
            traceback.print_exc()
            raise RuntimeError(f"An unexpected error occurred during bgzip compression: {e}")

        try:
            abs_gff_gz_filepath_fs = str(gff_gz_filepath_fs)
            abs_gff_tbi_filepath_fs = str(gff_tbi_filepath_fs)

            env = os.environ.copy()
            env['LC_ALL'] = 'C'

            tabix_command = [TABIX_PATH, '-f', '-p', 'gff', abs_gff_gz_filepath_fs] 
            print(f"DEBUG: Running tabix command: {' '.join(tabix_command)}")
            process = subprocess.run(
                tabix_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                encoding='utf-8',
                errors='replace'
            )
            print(f"DEBUG: tabix stdout: {process.stdout.strip()}")
            print(f"DEBUG: tabix stderr: {process.stderr.strip()}")
            
            if not gff_tbi_filepath_fs.exists() or gff_tbi_filepath_fs.stat().st_size == 0:
                error_message = f"Tabix indexing failed: .tbi file was not created or is empty ({abs_gff_tbi_filepath_fs}). "
                error_message += f"Tabix stdout: {process.stdout.strip()}\nTabix stderr: {process.stderr.strip()}"
                print(f"ERROR: {error_message}")
                raise RuntimeError(error_message)
            print(f"DEBUG: Tabix indexing successful to {gff_tbi_filepath_fs}")

        except FileNotFoundError:
            print(f"ERROR: tabix executable not found at {TABIX_PATH}. Please check the path configuration.")
            raise RuntimeError(f"tabix executable not found at {TABIX_PATH}. Please check the path configuration.")
        except subprocess.CalledProcessError as e:
            stdout_output = e.stdout.strip() if e.stdout else ""
            stderr_output = e.stderr.strip() if e.stderr else ""
            print(f"ERROR: tabix failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
            raise RuntimeError(f"tabix failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
        except Exception as e:
            print(f"ERROR: An unexpected error occurred during tabix indexing: {e}")
            traceback.print_exc()
            raise RuntimeError(f"An unexpected error occurred during tabix indexing: {e}")

@app.route('/update_gff_file', methods=['POST'])
def update_gff_file():
    print("DEBUG: /update_gff_file endpoint hit.")
    try:
        data = request.get_json()
        session_id = data.get('session_id')
        updated_features_data = data.get('features')

        if not session_id or not updated_features_data:
            print("ERROR: Missing session_id or features data in /update_gff_file request.")
            return jsonify({'error': 'Missing session_id or features data'}), 400

        session_data_dir_abs = (Path("static/igvjs/data") / session_id).resolve()
        
        sequence_id = session.get(f'{session_id}_sequence_id')
        sequence_length = session.get(f'{session_id}_sequence_length')

        if not sequence_id or not sequence_length:
            print("ERROR: Sequence metadata not found in session for GFF update.")
            return jsonify({'error': 'Sequence metadata not found for this session. Please re-upload files.'}), 400

        print(f"DEBUG: Converting {len(updated_features_data)} features from dict to Bio.SeqFeature objects.")
        dummy_seq_record = SeqRecord(Seq('N' * sequence_length), id=sequence_id) # Sequence content not needed here
        bio_features = [dict_to_feature(f, sequence_length) for f in updated_features_data]
        dummy_seq_record.features = bio_features
        print(f"DEBUG: Features converted successfully.")

        # Prepare features for frontend display (re-aggregate exons if necessary)
        features_for_frontend = _prepare_frontend_features(dummy_seq_record)
        print(f"DEBUG: Re-prepared {len(features_for_frontend)} features for frontend after update.")

        print(f"DEBUG: Generating GFF3 content from {len(bio_features)} Bio.SeqFeature objects.")
        gff_feature_lines_content = write_gff3(dummy_seq_record)
        print(f"DEBUG: GFF3 content generated successfully (length: {len(gff_feature_lines_content)}).")
        print("--- GFF3 Preview Start (first 10 lines) ---")
        for i, line in enumerate(gff_feature_lines_content.splitlines()):
            if i >= 10: break
            print(line)
        print("--- GFF3 Preview End ---")


        header_lines = [
            "##gff-version 3",
            f"##sequence-region {sequence_id} 1 {sequence_length}"
        ]

        gff_data_filename_base = "annotations.gff3"
        
        _write_gff_data_to_file_system(
            session_data_dir_abs, 
            gff_data_filename_base, 
            header_lines, 
            gff_feature_lines_content
        )
        print("DEBUG: GFF file updated successfully on server file system.")

        return jsonify({'success': True, 'message': 'GFF file updated successfully on server.', 'features': features_for_frontend})

    except Exception as e:
        print(f"ERROR: Server error during GFF update: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

def cleanup_igv(session_id):
    print(f"DEBUG: Cleaning up IGV data for session {session_id}")
    igv_data_dir = Path(app.static_folder) / 'igvjs' / 'data' / session_id
    if igv_data_dir.exists():
        try:
            shutil.rmtree(igv_data_dir, ignore_errors=True)
            print(f"DEBUG: Removed directory {igv_data_dir}")
        except Exception as e:
            print(f"WARNING: Could not remove directory {igv_data_dir}: {e}")
    # Removed session[f'{session_id}_sequence_data'] from here, as it's no longer stored
    session.pop(f'{session_id}_sequence_id', None)
    session.pop(f'{session_id}_sequence_length', None)
    print(f"DEBUG: Session variables cleaned for {session_id}.")


# --- REVISED parse_files function ---
def parse_files(gff_file_stream, fasta_file_stream):
    print("DEBUG: parse_files function started.")
    fasta_content = fasta_file_stream.read().decode('utf-8')
    try:
        seq_record = SeqIO.read(StringIO(fasta_content), 'fasta')
        print(f"DEBUG: FASTA file parsed. Sequence ID: {seq_record.id}, length: {len(seq_record.seq)}")
    except Exception as e:
        print(f"ERROR: FASTA parsing error: {e}")
        raise ValueError(f"FASTA parsing error: {str(e)}")
    
    seq_length = len(seq_record.seq)
    if seq_length == 0:
        print("ERROR: Empty sequence in FASTA file.")
        raise ValueError("Empty sequence in FASTA file")
    
    gff_content = gff_file_stream.read().decode('utf-8')
    print(f"DEBUG: GFF content read (length: {len(gff_content)}).")
    
    # Temporarily store raw parsed data grouped by feature ID
    temp_features_by_id = {} 
    
    for line_num, line in enumerate(gff_content.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        try:
            parts = line.split('\t')
            if len(parts) < 9:
                print(f"WARNING: Skipping malformed GFF3 line {line_num} (not enough columns): {line}")
                continue

            seqid, source, type_, start, end, score, strand_str, phase_str, attributes_str = parts

            start_pos = int(start) - 1
            end_pos = int(end)

            if start_pos < 0 or end_pos > seq_length or start_pos >= end_pos:
                print(f"WARNING: Invalid coordinates ({start}-{end}) or outside sequence length ({seq_length}) at line {line_num}: {line}")
                continue
                
            strand_val = 0
            if strand_str == '+':
                strand_val = 1
            elif strand_str == '-':
                strand_val = -1
            
            attr_dict = {}
            for attr_pair in attributes_str.split(';'):
                if '=' in attr_pair:
                    key, val_str = attr_pair.split('=', 1)
                    val = _unescape_gff3_value(val_str) 
                    attr_dict.setdefault(key, []).extend(val.split(',')) 
            
            if type_.lower() == 'cds' and phase_str != '.':
                try:
                    attr_dict.setdefault('codon_start', []).append(str(int(phase_str) + 1))
                except ValueError:
                    print(f"WARNING: Could not parse phase '{phase_str}' to int for CDS at line {line_num}. Setting codon_start to '1'.")
                    attr_dict.setdefault('codon_start', ['1']) 

            feature_id = attr_dict.get('ID', [''])[0]
            if not feature_id: 
                feature_id = attr_dict.get('Name', [''])[0]
            if not feature_id:
                 feature_id = f"{type_}-{start}-{end}-{str(uuid.uuid4())[:8]}"
                 print(f"DEBUG: Generated temporary ID '{feature_id}' for feature at line {line_num}.")
                 attr_dict['ID'] = [feature_id] 
            
            # Store raw parsed data for this line
            temp_features_by_id.setdefault(feature_id, []).append({
                'type': type_,
                'location': FeatureLocation(start_pos, end_pos, strand=strand_val),
                'qualifiers': attr_dict
            })
            
        except ValueError as ve:
            print(f"ERROR: Error parsing line {line_num} due to data validation: {ve}. Line: '{line}'")
        except Exception as e:
            print(f"ERROR: General error parsing line {line_num}: {e}. Line: '{line}'")
            traceback.print_exc() 

    final_consolidated_features = []
    print(f"DEBUG: Consolidating {len(temp_features_by_id)} unique IDs into Bio.SeqFeatures.")
    
    # Define type priority for deciding the primary type if multiple types share an ID
    type_priority_for_consolidation = {
        'gene': 0, 'mrna': 1, 'transcript': 1, 'ncrna': 1, 'rrna': 1, 'trna': 1, 
        'cds': 2, 'exon': 3, 'utr': 4, 'intron': 5
    }

    for feature_id, raw_parts_data in temp_features_by_id.items():
        if not raw_parts_data:
            continue

        # Sort raw parts to determine primary type and consistent location order
        raw_parts_data.sort(key=lambda x: (type_priority_for_consolidation.get(x['type'].lower(), 99), x['location'].start, x['location'].end))

        primary_feature_type = raw_parts_data[0]['type'] 
        
        all_locations_for_feature = []
        consolidated_qualifiers = {}

        # Merge qualifiers from all parts for this feature_id
        for part in raw_parts_data:
            all_locations_for_feature.append(part['location'])
            
            for key, values in part['qualifiers'].items():
                if key not in consolidated_qualifiers:
                    consolidated_qualifiers[key] = []
                for val in values:
                    if val not in consolidated_qualifiers[key]:
                        consolidated_qualifiers[key].append(val)
        
        # Ensure 'ID' is correctly set in consolidated qualifiers
        if 'ID' not in consolidated_qualifiers or not consolidated_qualifiers['ID']:
            consolidated_qualifiers['ID'] = [feature_id] 

        # Sort locations for CompoundLocation construction
        all_locations_for_feature.sort(key=lambda loc: (loc.start, loc.end))

        final_location = None
        if len(all_locations_for_feature) > 1:
            final_location = CompoundLocation(all_locations_for_feature)
            print(f"DEBUG: Consolidated feature {feature_id} into CompoundLocation with {len(all_locations_for_feature)} parts.")
        elif all_locations_for_feature:
            final_location = all_locations_for_feature[0]
            print(f"DEBUG: Consolidated feature {feature_id} has a single location.")
        else:
            print(f"WARNING: No valid locations for consolidated feature ID {feature_id}. Setting dummy location.")
            final_location = FeatureLocation(0, 0, strand=0)
        
        # Create the final consolidated Bio.SeqFeature
        new_feature = SeqFeature(
            location=final_location,
            type=primary_feature_type, 
            qualifiers=consolidated_qualifiers
        )
        new_feature.id = feature_id 
        final_consolidated_features.append(new_feature)
    
    seq_record.features = final_consolidated_features
    print(f"DEBUG: parse_files function finished. {len(final_consolidated_features)} consolidated features extracted.")
    return seq_record


# --- REVISED write_gff3 function ---
def write_gff3(seq_record):
    print("DEBUG: write_gff3 function started (hierarchical version).")
    gff_output_lines = []
    
    # 1. Build a map of all features by ID and parent-child relationships
    print(f"DEBUG: Total features received by write_gff3: {len(seq_record.features)}")
    if not seq_record.features:
        print("DEBUG: No features to process in write_gff3. Returning empty GFF.")
        return ""

    feature_map = {f.id: f for f in seq_record.features}
    print(f"DEBUG: Feature map IDs: {list(feature_map.keys())}")

    parent_to_children_map = {}
    top_level_features = []

    # Filter out 'region' features or others that shouldn't be printed as main lines
    # Also, ensure all features have an ID.
    printable_features = [f for f in seq_record.features if f.type.lower() != 'region' and f.id]
    print(f"DEBUG: Printable features (after filtering 'region' and no ID): {len(printable_features)}")

    for feature in printable_features:
        parent_ids = feature.qualifiers.get('Parent', [])
        
        is_top_level = True
        if parent_ids:
            found_internal_parent = False # Flag to check if any parent is *within* the currently parsed features
            for p_id in parent_ids:
                if p_id in feature_map: # If parent exists in our current set of features
                    parent_to_children_map.setdefault(p_id, []).append(feature)
                    found_internal_parent = True
                    # If a feature has multiple parents, we assume the first found internal parent is the primary one for this tree
                    break 
            if found_internal_parent:
                is_top_level = False
        
        if is_top_level:
            top_level_features.append(feature)

    print(f"DEBUG: Discovered {len(top_level_features)} top-level features.")
    print(f"DEBUG: Parent-to-children map: { {pid: [f.id for f in children] for pid, children in parent_to_children_map.items()} }")

    # Define a consistent order for feature types for sorting
    type_order = {
        'gene': 0, 'mrna': 1, 'transcript': 1, 'ncrna': 1, 'rrna': 1, 'trna': 1, 
        'utr': 2, 'intron': 2, 'cds': 3, 'exon': 4 
    }

    def get_sort_key_for_gff_output(f):
        return (type_order.get(f.type.lower(), 99), f.location.start, f.location.end)

    # Sort top-level features
    top_level_features.sort(key=get_sort_key_for_gff_output)

    # Recursive helper to write features and their children
    def _write_feature_and_children(feature_obj):
        # Ensure feature_obj has a valid ID before proceeding
        feature_id = feature_obj.id if feature_obj.id else feature_obj.qualifiers.get('ID', [''])[0]
        if not feature_id:
            print(f"WARNING: Skipping feature without a valid ID: {feature_obj.type} {feature_obj.location}. Cannot write to GFF.")
            return

        base_qualifiers = feature_obj.qualifiers.copy()
        
        start_coord = int(feature_obj.location.start) + 1
        end_coord = int(feature_obj.location.end)
        strand_str = '+' if feature_obj.location.strand == 1 else ('-' if feature_obj.location.strand == -1 else '.')
        
        phase = '.'
        if feature_obj.type.lower() == 'cds':
            try:
                codon_start_str = base_qualifiers.get('codon_start', ['1'])[0]
                phase = str(int(codon_start_str) - 1)
            except ValueError:
                pass # Default to '.'

        attributes = []
        attributes.append(f"ID={_escape_gff3_value(feature_id)}")

        # Add Parent attribute
        if 'Parent' in base_qualifiers and base_qualifiers['Parent']:
            parent_values = [str(val).strip() for val in base_qualifiers['Parent']]
            attributes.append(f"Parent={','.join(parent_values)}")
        
        # Add Name (preferred over gene/product for display)
        display_name = base_qualifiers.get('Name', [''])[0]
        if not display_name and base_qualifiers.get('gene'):
            display_name = base_qualifiers['gene'][0]
        if not display_name and base_qualifiers.get('product'):
            display_name = base_qualifiers['product'][0]
        if display_name:
            attributes.append(f"Name={_escape_gff3_value(display_name)}")

        # Add all other qualifiers, escaping values
        for key in sorted([k for k in base_qualifiers.keys() if k not in ['ID', 'Parent', 'Name', 'codon_start', 'part']]):
            values = base_qualifiers[key]
            processed_values = [_escape_gff3_value(str(v)) for v in values]
            if processed_values:
                attributes.append(f"{key}={','.join(processed_values)}")

        main_fields = [
            seq_record.id,
            'GenomeAnnotationTool',
            feature_obj.type,
            str(start_coord),
            str(end_coord),
            '.', # Score
            strand_str,
            phase
        ]
        gff_output_lines.append('\t'.join(main_fields) + '\t' + ';'.join(attributes))
        print(f"DEBUG: Wrote GFF line for {feature_obj.type}: ID={feature_id}, Loc={start_coord}-{end_coord}")

        # 2. If this feature has a CompoundLocation and is an mRNA or CDS, explicitly write 'exon' parts
        # This handles cases where exons are not explicitly provided as separate features but implied by compound locations
        # Also ensures exons are parented correctly.
        if isinstance(feature_obj.location, CompoundLocation) and feature_obj.type.lower() in ['mrna', 'transcript', 'cds']:
            locations_to_write_as_exons = list(feature_obj.location.parts)
            locations_to_write_as_exons.sort(key=lambda x: x.start)

            cumulative_len_for_phase = 0 
            if feature_obj.type.lower() == 'cds':
                try:
                    codon_start_val = int(base_qualifiers.get('codon_start', ['1'])[0]) - 1
                except ValueError:
                    codon_start_val = 0 # Default if codon_start is bad
            else:
                codon_start_val = 0 # Not applicable for non-CDS features

            for i, loc_part in enumerate(locations_to_write_as_exons):
                exon_start = int(loc_part.start) + 1
                exon_end = int(loc_part.end)
                exon_strand_str = '+' if loc_part.strand == 1 else ('-' if loc_part.strand == -1 else '.')
                
                exon_attributes = []
                exon_id = f"exon-{feature_id}.{i+1}" # A consistent ID for the exon part
                exon_attributes.append(f"ID={_escape_gff3_value(exon_id)}")
                # Parent the exon to the mRNA (or CDS if no mRNA parent)
                parent_for_exon = feature_id # Parent it to the feature that defines its segments (mRNA or CDS)
                exon_attributes.append(f"Parent={_escape_gff3_value(parent_for_exon)}") 

                # Inherit Name, gene, product, locus_tag, etc., from the parent feature for exons
                for key in ['Name', 'gene', 'product', 'locus_tag', 'Dbxref', 'Note', 'transcript_id', 'protein_id']:
                     if key in base_qualifiers and base_qualifiers[key]:
                        values = base_qualifiers[key]
                        processed_values = [_escape_gff3_value(str(v)) for v in values]
                        if processed_values:
                            exon_attributes.append(f"{key}={','.join(processed_values)}")
                
                exon_attributes.append(f"part={i+1}") 

                exon_phase = '.'
                if feature_obj.type.lower() == 'cds':
                    if len(loc_part) == 0:
                        print(f"WARNING: Zero-length exon part ({loc_part.start}-{loc_part.end}) for CDS {feature_id}. Skipping phase calculation for this part.")
                        exon_phase = '.' 
                    else:
                        calculated_phase_for_exon = (codon_start_val + cumulative_len_for_phase) % 3
                        exon_phase = str(calculated_phase_for_exon)
                        cumulative_len_for_phase += len(loc_part) 

                exon_fields = [
                    seq_record.id,
                    'GenomeAnnotationTool',
                    'exon', # Explicitly use 'exon' type for segments
                    str(exon_start),
                    str(exon_end),
                    '.', # Score
                    exon_strand_str,
                    exon_phase 
                ]
                gff_output_lines.append('\t'.join(exon_fields) + '\t' + ';'.join(exon_attributes))
                print(f"DEBUG: Wrote generated exon line for {feature_obj.type}: ID={exon_id}, Parent={parent_for_exon}, Loc={exon_start}-{exon_end}")

        # 3. Recursively write children
        children = parent_to_children_map.get(feature_id, [])
        children.sort(key=get_sort_key_for_gff_output) # Sort children for consistent output order

        for child_feature in children:
            _write_feature_and_children(child_feature)

    # Start the recursive writing from top-level features
    if not top_level_features:
        print("WARNING: No top-level features found. This might mean all features have parents, or input GFF is not hierarchical.")
        # Fallback: if no top-level features are found by Parent attribute,
        # iterate through all features and treat them as top-level for printing
        # This prevents empty GFF if Parent attributes are missing or mislinked.
        print("DEBUG: Attempting to print all printable features as top-level due to no detected hierarchy.")
        for feature in printable_features:
            # Only print if it hasn't been printed as a child already
            # This is a heuristic and might lead to duplicate if hierarchy is partially broken
            # A more robust solution might involve explicitly tracking printed IDs.
            is_already_child = False
            for child_list in parent_to_children_map.values():
                if feature in child_list:
                    is_already_child = True
                    break
            if not is_already_child:
                _write_feature_and_children(feature)
    else:
        for feature in top_level_features:
            _write_feature_and_children(feature)

    print(f"DEBUG: write_gff3 function finished. Total lines generated: {len(gff_output_lines)}")
    return '\n'.join(gff_output_lines)


def _prepare_frontend_features(seq_record):
    """
    Prepares a list of feature dictionaries for the frontend.
    This relies on parse_files having already consolidated multi-segment features
    into single Bio.SeqFeature objects with CompoundLocation.
    """
    print(f"DEBUG: _prepare_frontend_features started. Processing {len(seq_record.features)} consolidated features.")
    
    features_for_frontend = []
    
    # We now expect seq_record.features to contain consolidated Bio.SeqFeatures,
    # where multi-segment features (like spliced CDS) will have CompoundLocation.
    for feature in seq_record.features:
        # Skip 'region' features or any that are clearly not main annotations
        if feature.type.lower() == 'region':
            print(f"DEBUG: Skipping 'region' feature {feature.id} for frontend list.")
            continue 

        # Create the basic feature dictionary
        feature_dict = {
            'id': feature.id, 
            'type': feature.type,
            'start': int(feature.location.start) + 1, 
            'end': int(feature.location.end),       
            'strand': int(feature.location.strand if feature.location.strand is not None else 0), 
            
            'gene': feature.qualifiers.get('gene', [''])[0],
            'product': feature.qualifiers.get('product', [''])[0],
            'note': feature.qualifiers.get('Note', [''])[0],
            'codon_start': feature.qualifiers.get('codon_start', ['1'])[0], 
            'transcript_id': feature.qualifiers.get('transcript_id', [''])[0],
            'protein_id': feature.qualifiers.get('protein_id', [''])[0],
            'locus_tag': feature.qualifiers.get('locus_tag', [''])[0],
            'db_xref': ','.join(feature.qualifiers.get('Dbxref', [])) 
        }

        parent_ids = feature.qualifiers.get('Parent', [])
        feature_dict['parentId'] = parent_ids[0] if parent_ids else ''

        # Populate the 'exons' list based on the feature's location
        exons_list = []
        if isinstance(feature.location, CompoundLocation):
            for loc_part in sorted(feature.location.parts, key=lambda x: x.start):
                exons_list.append({
                    'start': int(loc_part.start) + 1,
                    'end': int(loc_part.end),
                    'strand': int(loc_part.strand if loc_part.strand is not None else 0) 
                })
            print(f"DEBUG: Feature {feature.id} (type: {feature.type}) had CompoundLocation with {len(exons_list)} parts for frontend 'exons'.")
        else: # Single location feature (e.g., non-spliced CDS, gene, tRNA, rRNA)
            exons_list.append({
                'start': int(feature.location.start) + 1,
                'end': int(feature.location.end),
                'strand': int(feature.location.strand if feature.location.strand is not None else 0) 
            })
            print(f"DEBUG: Feature {feature.id} (type: {feature.type}) is single-segment, added its single location to frontend 'exons'.")
        
        feature_dict['exons'] = exons_list
        features_for_frontend.append(feature_dict)

    print(f"DEBUG: _prepare_frontend_features finished. Prepared {len(features_for_frontend)} features for frontend.")
    return features_for_frontend


def feature_to_dict(feature):
    """
    Converts a Bio.SeqFeature object to a basic dictionary for frontend.
    This version is simplified as _prepare_frontend_features now handles aggregation.
    This is mainly used in old contexts or if not processed by _prepare_frontend_features.
    Ideally, _prepare_frontend_features replaces the direct use of this.
    """
    qualifiers = feature.qualifiers.copy()
    
    feature_strand = feature.location.strand if feature.location.strand is not None else 0

    feature_dict = {
        'id': feature.id, 
        'type': feature.type,
        'start': int(feature.location.start) + 1, 
        'end': int(feature.location.end),       
        'strand': int(feature_strand), 
        
        # Qualifiers that are directly displayed/edited
        'gene': qualifiers.get('gene', [''])[0],
        'product': qualifiers.get('product', [''])[0],
        'note': qualifiers.get('Note', [''])[0],
        'codon_start': qualifiers.get('codon_start', ['1'])[0], 
        'transcript_id': qualifiers.get('transcript_id', [''])[0],
        'protein_id': qualifiers.get('protein_id', [''])[0],
        'locus_tag': qualifiers.get('locus_tag', [''])[0],
        'db_xref': ','.join(qualifiers.get('Dbxref', [])) 
    }

    # Handle Parent attribute
    parent_ids = qualifiers.get('Parent', [])
    if parent_ids:
        feature_dict['parentId'] = parent_ids[0] 
    else:
        feature_dict['parentId'] = '' 

    # For this simplified version, just add the main location as an "exon"
    # The _prepare_frontend_features function will do the full aggregation.
    exons_list = [{
        'start': int(feature.location.start) + 1,
        'end': int(feature.location.end),
        'strand': int(feature_strand) 
    }]
    feature_dict['exons'] = exons_list
    return feature_dict


def dict_to_feature(feature_dict, sequence_length):
    """Converts a dictionary from frontend to a Bio.SeqFeature object."""
    
    exon_locations = []
    # Build CompoundLocation if multiple exons are provided from the frontend
    if 'exons' in feature_dict and feature_dict['exons']:
        for exon in feature_dict['exons']:
            start = max(0, int(exon['start']) - 1) 
            end = min(sequence_length, int(exon['end'])) 
            strand = int(exon['strand']) if exon['strand'] is not None else 0
            # Only add non-empty locations
            if start < end:
                exon_locations.append(FeatureLocation(start, end, strand=strand))
            else:
                print(f"WARNING: Skipping invalid (empty or negative length) exon from dict for feature {feature_dict.get('id')}: {exon}")
    else: # If no 'exons' key or empty, use the main feature's start/end
        start = max(0, int(feature_dict['start']) - 1)
        end = min(sequence_length, int(feature_dict['end']))
        strand = int(feature_dict['strand']) if feature_dict['strand'] is not None else 0
        if start < end:
            exon_locations.append(FeatureLocation(start, end, strand=strand))
        else:
             print(f"WARNING: Skipping invalid (empty or negative length) main location from dict for feature {feature_dict.get('id')}: {feature_dict['start']}-{feature_dict['end']}")


    exon_locations.sort(key=lambda loc: loc.start)

    location = None
    if len(exon_locations) > 1:
        location = CompoundLocation(exon_locations)
    elif exon_locations:
        location = exon_locations[0]
    else: 
        print(f"WARNING: No valid locations found for feature {feature_dict.get('id')}. Setting dummy location.")
        location = FeatureLocation(0, 0, strand=0) # Fallback dummy location
    
    qualifiers = {}
    
    # Store relevant attributes in qualifiers (note: BioPython's qualifiers are lists of values)
    if feature_dict.get('id'): qualifiers['ID'] = [feature_dict['id']]
    if feature_dict.get('gene'): qualifiers['gene'] = [feature_dict['gene']]
    if feature_dict.get('product'): qualifiers['product'] = [feature_dict['product']]
    if feature_dict.get('note'): qualifiers['Note'] = [feature_dict['note']]
    if feature_dict.get('locus_tag'): qualifiers['locus_tag'] = [feature_dict['locus_tag']]
    if feature_dict.get('transcript_id'): qualifiers['transcript_id'] = [feature_dict['transcript_id']]
    if feature_dict.get('protein_id'): qualifiers['protein_id'] = [feature_dict['protein_id']]
    if feature_dict.get('db_xref'): 
        # Split by comma and clean for multiple Db_xrefs
        qualifiers['Dbxref'] = [x.strip() for x in feature_dict['db_xref'].split(',') if x.strip()] 
    
    if feature_dict.get('parentId'):
        qualifiers['Parent'] = [feature_dict['parentId']] # Store Parent as list

    if feature_dict['type'].lower() == 'cds' and feature_dict.get('codon_start'):
        qualifiers['codon_start'] = [str(feature_dict['codon_start'])]

    feature = SeqFeature(
        location=location,
        type=feature_dict['type'],
        qualifiers=qualifiers
    )
    
    feature.id = feature_dict['id'] # Assign Biopython .id attribute
    return feature

def reverse_complement_features(seq_record):
    print("DEBUG: reverse_complement_features started.")
    if len(seq_record.seq) == 0:
        raise ValueError("Cannot reverse complement empty sequence")
    
    seq_record.seq = seq_record.seq.reverse_complement()
    seq_length = len(seq_record.seq)

    for feature in seq_record.features:
        locations_to_process = []
        if isinstance(feature.location, CompoundLocation):
            locations_to_process.extend(feature.location.parts)
        else:
            locations_to_process.append(feature.location)

        new_parts = []
        for loc in locations_to_process:
            new_start = seq_length - loc.end
            new_end = seq_length - loc.start
            new_strand = -1 * (loc.strand if loc.strand is not None else 0) 
            new_parts.append(FeatureLocation(new_start, new_end, strand=new_strand))
        
        new_parts.sort(key=lambda x: x.start)

        if len(new_parts) > 1:
            feature.location = CompoundLocation(new_parts)
        else:
            feature.location = new_parts[0]
        
        if feature.type.lower() == 'cds' and 'codon_start' in feature.qualifiers:
            feature.qualifiers['codon_start'] = ['1'] 

    return seq_record

def shift_origin(seq_record, new_start):
    print("DEBUG: shift_origin started.")
    if len(seq_record.seq) == 0:
        raise ValueError("Cannot shift origin of empty sequence")
    
    seq_length = len(seq_record.seq)
    new_start_0based = new_start - 1 
    
    if not (0 <= new_start_0based < seq_length):
        raise ValueError(f"Origin position must be between 1 and {seq_length}")
    
    shifted_seq = seq_record.seq[new_start_0based:] + seq_record.seq[:new_start_0based]
    seq_record.seq = shifted_seq
    
    for feature in seq_record.features:
        locations_to_process = []
        if isinstance(feature.location, CompoundLocation):
            locations_to_process.extend(feature.location.parts)
        else:
            locations_to_process.append(feature.location)

        new_parts = []
        for loc in locations_to_process:
            original_start = loc.start
            original_end = loc.end
            
            relative_start = original_start - new_start_0based
            relative_end = original_end - new_start_0based

            new_exon_start = (relative_start + seq_length) % seq_length
            new_exon_end = (relative_end + seq_length) % seq_length

            loc_strand = loc.strand if loc.strand is not None else 0
            new_parts.append(FeatureLocation(new_exon_start, new_exon_end, strand=loc_strand))
        
        new_parts.sort(key=lambda x: x.start)

        if len(new_parts) > 1:
            feature.location = CompoundLocation(new_parts)
        else:
            feature.location = new_parts[0]
        
        if feature.type.lower() == 'cds' and 'codon_start' in feature.qualifiers:
            feature.qualifiers['codon_start'] = ['1'] 
            
    print("DEBUG: shift_origin finished.")
    return seq_record

def validate_cds(feature, seq_record_dummy, codon_table_id): # Renamed seq_record to seq_record_dummy
    print(f"DEBUG: validate_cds started for feature {feature.id}, type {feature.type}.")
    if feature.type.lower() != 'cds':
        return {'valid': True, 'message': 'Not a CDS feature'}
    
    codon_table = next(
        (t for t in CODON_TABLES if t['id'] == codon_table_id),
        None
    )
    if not codon_table:
        print(f"ERROR: Invalid codon table ID '{codon_table_id}' selected for CDS validation.")
        return {'valid': False, 'message': 'Invalid codon table selected'}
    
    # Read sequence data from the file system, not session
    session_id = session.get('last_session_id_for_validation') # Stored during upload
    if not session_id:
        print("ERROR: Session ID not found for reading sequence data during validation.")
        return {'valid': False, 'message': 'Session expired or sequence data not found.'}
    
    fasta_filepath_fs = Path("static/igvjs/data") / session_id / "sequence.fasta"
    try:
        current_seq_record = SeqIO.read(fasta_filepath_fs, 'fasta')
        print(f"DEBUG: Sequence loaded from file for CDS validation.")
    except Exception as e:
        print(f"ERROR: Could not load sequence from {fasta_filepath_fs} for validation: {e}")
        return {'valid': False, 'message': f"Could not load sequence data: {e}"}

    cds_seq_obj = Seq("")
    all_locations = []
    if isinstance(feature.location, CompoundLocation):
        all_locations.extend(feature.location.parts)
    else:
        all_locations.append(feature.location)
    
    all_locations.sort(key=lambda x: x.start)
    print(f"DEBUG: CDS {feature.id} has {len(all_locations)} parts.")

    for exon_loc in all_locations:
        start_bound = max(0, exon_loc.start)
        end_bound = min(len(current_seq_record.seq), exon_loc.end) # Use current_seq_record length
        
        if start_bound >= end_bound:
            print(f"WARNING: Skipping zero or negative length exon part ({start_bound}-{end_bound}) for CDS {feature.id}.")
            continue 

        exon_seq = current_seq_record.seq[start_bound:end_bound] # Use current_seq_record
        
        if exon_loc.strand == -1: 
            exon_seq = exon_seq.reverse_complement()
        cds_seq_obj += exon_seq
    
    cds_seq = str(cds_seq_obj).upper()
    print(f"DEBUG: Translated CDS sequence length: {len(cds_seq)}")
    
    codon_start_val = 1
    try:
        if 'codon_start' in feature.qualifiers:
            codon_start_val = int(feature.qualifiers['codon_start'][0])
    except ValueError:
        print(f"WARNING: Could not parse codon_start for CDS {feature.id}. Defaulting to 1.")
        pass 

    if not (1 <= codon_start_val <= 3):
        print(f"ERROR: Invalid codon_start value for CDS {feature.id}: {codon_start_val}.")
        return {'valid': False, 'message': f"Invalid codon_start value: {codon_start_val}. Must be 1, 2, or 3."}

    cds_seq_framed = cds_seq[codon_start_val - 1:]

    errors = []
    
    if len(cds_seq_framed) < 3:
        errors.append(f"CDS length ({len(cds_seq_framed)}) too short for translation after applying codon_start. Must be at least 3bp.")
    elif len(cds_seq_framed) % 3 != 0:
        errors.append(f"CDS length ({len(cds_seq_framed)}) is not divisible by 3 after applying codon_start (remainder {len(cds_seq_framed) % 3}).")
    
    if len(cds_seq_framed) >= 3:
        start_codon = cds_seq_framed[:3]
        if start_codon not in codon_table['start_codons']:
            errors.append(f"Invalid start codon: {start_codon} (at relative position 1). Expected one of: {', '.join(codon_table['start_codons'])}")
        
        stop_codon = cds_seq_framed[-3:]
        if stop_codon not in codon_table['stop_codons']:
            errors.append(f"Invalid stop codon: {stop_codon} (at relative position {len(cds_seq_framed) - 2}). Expected one of: {', '.join(codon_table['stop_codons'])}")
        
        for i in range(3, len(cds_seq_framed) - 3, 3): 
            codon = cds_seq_framed[i:i+3]
            if codon in codon_table['stop_codons']:
                errors.append(f"Internal stop codon found at relative position {i+1} ({codon}).")
    
    print(f"DEBUG: CDS validation finished for {feature.id}. Errors: {errors if errors else 'None'}")    
    return {
        'valid': len(errors) == 0,
        'message': 'CDS is valid' if not errors else '; '.join(errors)
    }

@app.route('/')
def index():
    print("DEBUG: Serving index.html")
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    print("DEBUG: /upload endpoint hit.")
    if 'gff_file' not in request.files or 'fasta_file' not in request.files:
        print("ERROR: Missing gff_file or fasta_file in upload request.")
        return jsonify({'error': 'Missing files'}), 400
    
    try:
        gff_file = request.files['gff_file']
        fasta_file = request.files['fasta_file']
        print(f"DEBUG: Received files: GFF={gff_file.filename}, FASTA={fasta_file.filename}")
        
        seq_record = parse_files(gff_file, fasta_file)
        print(f"DEBUG: Files parsed into SeqRecord successfully.")
        
        # Call the new function to prepare features for frontend
        features_for_frontend = _prepare_frontend_features(seq_record)
        print(f"DEBUG: Features prepared for frontend: {len(features_for_frontend)} features.")

        gff_feature_lines_content = write_gff3(seq_record)
        print(f"DEBUG: GFF3 content regenerated for server storage. Length: {len(gff_feature_lines_content)}")
        print("--- GFF3 Preview Start (first 10 lines) ---")
        for i, line in enumerate(gff_feature_lines_content.splitlines()):
            if i >= 10: break
            print(line)
        print("--- GFF3 Preview End ---")

        session_id = str(uuid.uuid4())
        print(f"DEBUG: Generated new session ID: {session_id}")
        
        session[f'{session_id}_sequence_id'] = seq_record.id
        session[f'{session_id}_sequence_length'] = len(seq_record.seq)
        # REMOVED: session[f'{session_id}_sequence_data'] = str(seq_record.seq) -- AVOID LARGE COOKIES
        session['last_session_id_for_validation'] = session_id # Store for later sequence access
        print(f"DEBUG: Essential sequence metadata stored in session.")

        response_data = setup_igv_data_and_config(
            seq_record, 
            gff_feature_lines_content, 
            session_id=session_id
        )
        
        response_data['features'] = features_for_frontend
        print(f"DEBUG: Upload successful. Returning response.")
        
        return jsonify(response_data)
    except Exception as e:
        print(f"ERROR: Server error during upload: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/static/igvjs/data/<path:filename>')
def igv_data(filename):
    print(f"DEBUG: Serving IGV data file: {filename}")
    base_data_dir = os.path.join(app.static_folder, 'igvjs', 'data')
    
    try:
        response = send_from_directory(base_data_dir, filename)
        
        if USE_UNCOMPRESSED_GFF and filename.endswith('.gff3') and not filename.endswith('.gff3.gz'):
            response.headers['Content-Type'] = 'text/plain' 
            print(f"DEBUG: Served {filename} as text/plain.")
        elif filename.endswith('.gz'):
            response.headers['Content-Type'] = 'application/x-gzip'
            print(f"DEBUG: Served {filename} as application/x-gzip.")
        elif filename.endswith('.tbi'):
            response.headers['Content-Type'] = 'application/octet-stream'
            print(f"DEBUG: Served {filename} as application/octet-stream.")
        else:
            print(f"DEBUG: Served {filename} with default mimetype.")
        
        return response
    except Exception as e:
        print(f"ERROR: Error serving IGV data file {filename}: {e}")
        traceback.print_exc()
        return "Error serving file", 500

@app.route('/validate_cds', methods=['POST'])
def validate_cds_feature():
    print("DEBUG: /validate_cds endpoint hit.")
    try:
        data = request.get_json()
        if not data or 'feature' not in data or 'codon_table' not in data:
            print("ERROR: Invalid request data for /validate_cds.")
            return jsonify({'error': 'Invalid request data'}), 400
        
        feature_dict = data['feature']
        codon_table_id = data['codon_table']
        
        session_id_from_client = data.get('session_id') 
        if not session_id_from_client:
            print("ERROR: Session ID missing from client for CDS validation.")
            return jsonify({'error': 'Session ID missing for validation'}), 400

        # Set the session ID for internal use in validate_cds function
        session['last_session_id_for_validation'] = session_id_from_client

        # Sequence data will be read by validate_cds function directly from file
        # The seq_record parameter in validate_cds function is now a dummy.
        dummy_seq_record = SeqRecord(Seq(''), id=feature_dict.get('sequence_id', "validation_seq")) 

        feature_obj = dict_to_feature(feature_dict, data.get('sequence_length')) # Pass sequence_length from frontend if available
        result = validate_cds(feature_obj, dummy_seq_record, codon_table_id) 
        print(f"DEBUG: CDS validation result: {result}")
        return jsonify(result)
    except Exception as e:
        print(f"ERROR: Server error during CDS validation: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/flip', methods=['POST'])
def flip_sequence():
    print("DEBUG: /flip endpoint hit.")
    try:
        data = request.get_json()
        if not data or 'sequence_length' not in data or 'features' not in data:
            print("ERROR: Invalid request data for /flip.")
            return jsonify({'error': 'Invalid request data'}), 400
        
        session_id = data.get('session_id') 
        if not session_id:
            print("ERROR: Session ID missing for flip operation.")
            return jsonify({'error': 'Session ID missing for flip operation'}), 400

        # Read sequence from file for flip
        fasta_filepath_fs = Path("static/igvjs/data") / session_id / "sequence.fasta"
        try:
            seq_record = SeqIO.read(fasta_filepath_fs, 'fasta')
            print(f"DEBUG: Sequence loaded from file for flip operation.")
        except Exception as e:
            print(f"ERROR: Could not load sequence from {fasta_filepath_fs} for flip: {e}")
            return jsonify({'error': f"Could not load sequence data for flip: {e}"}), 500


        seq_record.features = [dict_to_feature(f, data['sequence_length']) for f in data['features']]
        print(f"DEBUG: Performing reverse complement on sequence and {len(seq_record.features)} features.")
        
        seq_record = reverse_complement_features(seq_record)
        
        # Write updated sequence back to file
        with open(fasta_filepath_fs, 'w') as f:
            SeqIO.write(seq_record, f, 'fasta')
        # Update .fai as well
        fai_filepath_fs = fasta_filepath_fs.with_suffix('.fasta.fai')
        with open(fai_filepath_fs, "w") as f:
            seq_len = len(seq_record.seq)
            line_length = 80 # Standard FASTA line length
            bytes_per_line = line_length + 1 # +1 for newline character
            offset = len(f">{seq_record.id}\n") 
            f.write(f"{seq_record.id}\t{seq_len}\t{offset}\t{line_length}\t{bytes_per_line}\n")
        
        # Update session for consistency (though actual seq content isn't there)
        session[f'{session_id}_sequence_id'] = seq_record.id
        session[f'{session_id}_sequence_length'] = len(seq_record.seq)


        features = _prepare_frontend_features(seq_record) # Use the new aggregation function
        print(f"DEBUG: Flip successful. Returning {len(features)} features (after re-aggregation).")
        
        return jsonify({
            'sequence_id': seq_record.id,
            'sequence_length': len(seq_record.seq),
            'features': features
        })
    except Exception as e:
        print(f"ERROR: Server error during flip: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/shift_origin', methods=['POST'])
def shift_origin_route():
    print("DEBUG: /shift_origin endpoint hit.")
    try:
        data = request.get_json()
        if not data or 'sequence_length' not in data or 'features' not in data or 'new_start' not in data:
            print("ERROR: Invalid request data for /shift_origin.")
            return jsonify({'error': 'Invalid request data'}), 400
        
        session_id = data.get('session_id') 
        if not session_id:
            print("ERROR: Session ID missing for shift origin operation.")
            return jsonify({'error': 'Session ID missing for shift origin operation'}), 400

        # Read sequence from file for shift origin
        fasta_filepath_fs = Path("static/igvjs/data") / session_id / "sequence.fasta"
        try:
            seq_record = SeqIO.read(fasta_filepath_fs, 'fasta')
            print(f"DEBUG: Sequence loaded from file for shift origin operation.")
        except Exception as e:
            print(f"ERROR: Could not load sequence from {fasta_filepath_fs} for shift origin: {e}")
            return jsonify({'error': f"Could not load sequence data for shift origin: {e}"}), 500


        seq_record.features = [dict_to_feature(f, data['sequence_length']) for f in data['features']]
        
        new_start = int(data['new_start'])
        print(f"DEBUG: Shifting origin to {new_start} for sequence and {len(seq_record.features)} features.")
        seq_record = shift_origin(seq_record, new_start)
        
        # Write updated sequence back to file
        with open(fasta_filepath_fs, 'w') as f:
            SeqIO.write(seq_record, f, 'fasta')
        # Update .fai as well
        fai_filepath_fs = fasta_filepath_fs.with_suffix('.fasta.fai')
        with open(fai_filepath_fs, "w") as f:
            seq_len = len(seq_record.seq)
            line_length = 80 # Standard FASTA line length
            bytes_per_line = line_length + 1 # +1 for newline character
            offset = len(f">{seq_record.id}\n") 
            f.write(f"{seq_record.id}\t{seq_len}\t{offset}\t{line_length}\t{bytes_per_line}\n")

        # Update session for consistency
        session[f'{session_id}_sequence_id'] = seq_record.id
        session[f'{session_id}_sequence_length'] = len(seq_record.seq)


        features = _prepare_frontend_features(seq_record) # Use the new aggregation function
        print(f"DEBUG: Shift origin successful. Returning {len(features)} features (after re-aggregation).")
        
        return jsonify({
            'sequence_id': seq_record.id,
            'sequence_length': len(seq_record.seq),
            'features': features
        })
    except Exception as e:
        print(f"ERROR: Server error during shift origin: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/export', methods=['POST'])
def export_files():
    print("DEBUG: /export endpoint hit.")
    try:
        data = request.get_json()
        if not data or 'sequence_length' not in data or 'features' not in data:
            print("ERROR: Invalid request data for /export.")
            return jsonify({'error': 'Invalid request data'}), 400
        
        sequence_id = data.get('sequence_id', 'unknown')
        sequence_length = data['sequence_length']

        session_id = data.get('session_id') 
        fasta_seq_str = ""
        if session_id:
            # Read actual sequence from file if session ID available
            fasta_filepath_fs = Path("static/igvjs/data") / session_id / "sequence.fasta"
            try:
                seq_record_from_file = SeqIO.read(fasta_filepath_fs, 'fasta')
                fasta_seq_str = str(seq_record_from_file.seq)
            except Exception as e:
                print(f"WARNING: Could not load sequence from {fasta_filepath_fs} for export. Using N's. Error: {e}")
                fasta_seq_str = "N" * sequence_length 
        else:
            fasta_seq_str = "N" * sequence_length 
            
        print(f"DEBUG: Exporting sequence ID: {sequence_id}, length: {sequence_length}")

        # The features sent from the frontend are already in the aggregated format,
        # so dict_to_feature will reconstruct CompoundLocation from its 'exons' list.
        seq_record = SeqRecord(
            Seq(fasta_seq_str),
            id=sequence_id,
            features=[dict_to_feature(f, sequence_length) for f in data['features']]
        )
        
        gff_output_features_only = write_gff3(seq_record)
        fasta_output = StringIO()
        SeqIO.write(seq_record, fasta_output, 'fasta')
        
        full_gff_export = "##gff-version 3\n"
        full_gff_export += f"##sequence-region {seq_record.id} 1 {len(seq_record.seq)}\n"
        full_gff_export += gff_output_features_only
        
        print("DEBUG: Export files generated successfully.")
        return jsonify({
            'gff': full_gff_export,
            'fasta': fasta_output.getvalue()
        })
    except Exception as e:
        print(f"ERROR: Server error during export: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/cleanup', methods=['POST'])
def cleanup_session():
    print("DEBUG: /cleanup endpoint hit.")
    try:
        data = request.get_json()
        if not data or 'session_id' not in data:
            print("ERROR: Session ID required for cleanup.")
            return jsonify({'error': 'Session ID required'}), 400
            
        cleanup_igv(data['session_id'])
        # Remove the 'last_session_id_for_validation' key if it matches the current session
        if session.get('last_session_id_for_validation') == data['session_id']:
            session.pop('last_session_id_for_validation', None)
        print(f"DEBUG: Cleanup successful for session {data['session_id']}.")
        return jsonify({'success': True})
    except Exception as e:
        print(f"ERROR: Server error during cleanup: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    codon_tables_path = Path('static/codon_tables.json')
    if not codon_tables_path.exists():
        codon_tables_path.parent.mkdir(parents=True, exist_ok=True)
        dummy_codon_data = [
            {"id": "1", "name": "Standard", "start": ["ATG"], "stop": ["TAA", "TAG", "TGA"], "table": {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "VTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}}
        ]
        with open(codon_tables_path, 'w', encoding='utf-8') as f:
            json.dump(dummy_codon_data, f, indent=4)
        print(f"Created dummy {codon_tables_path} for initial setup.")

    app.run(debug=True)
