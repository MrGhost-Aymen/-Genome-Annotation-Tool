import os
import json
from io import StringIO
from flask import Flask, render_template, request, jsonify, send_file, redirect, url_for, send_from_directory, Response
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

BGZIP_PATH = r"C:\msys64\mingw664\bin\bgzip.exe"
TABIX_PATH = r"C:\msys64\mingw664\bin\tabix.exe"
SORT_PATH = r"C:\msys64\usr\bin\sort.exe"

mimetypes.add_type('application/x-gzip', '.gz')
mimetypes.add_type('application/octet-stream', '.tbi')

app = Flask(__name__, static_folder='static', static_url_path='/static')
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

def load_codon_tables():
    try:
        with open('static/codon_tables.json', 'r', encoding='utf-8') as f:
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
        print(f"Error loading codon tables: {str(e)}")
        return []

CODON_TABLES = load_codon_tables()

@app.route('/api/codon_tables')
def get_codon_tables():
    return jsonify(CODON_TABLES)

def setup_jbrowse_data_with_external_tools(seq_record, gff_feature_lines_content, session_id):
    jbrowse_data_base_dir = Path("static/jbrowse2/data")
    session_data_dir = jbrowse_data_base_dir / session_id
    session_data_dir.mkdir(parents=True, exist_ok=True)
    session_data_dir_abs = session_data_dir.resolve()

    fasta_data_filename = "sequence.fasta"
    fasta_filepath_fs = session_data_dir_abs / fasta_data_filename
    with open(fasta_filepath_fs, 'w') as f:
        SeqIO.write(seq_record, f, 'fasta')
    
    fai_filepath_fs = fasta_filepath_fs.with_suffix('.fasta.fai')
    with open(fai_filepath_fs, "w") as f:
        seq_len = len(seq_record.seq)
        line_length = 80
        bytes_per_line = line_length + 1
        offset = len(f">{seq_record.id}\n")
        
        f.write(f"{seq_record.id}\t{seq_len}\t{offset}\t{line_length}\t{bytes_per_line}\n")

    gff_data_filename_base = "annotations.gff3"
    gff_gz_filepath_fs = session_data_dir_abs / (gff_data_filename_base + ".gz")
    gff_tbi_filepath_fs = session_data_dir_abs / (gff_data_filename_base + ".gz.tbi")

    header_lines = [
        "##gff-version 3",
        f"##sequence-region {seq_record.id} 1 {len(seq_record.seq)}"
    ]
    
    cleaned_feature_lines = [line.strip() for line in gff_feature_lines_content.splitlines() if line.strip()]
    
    header_file = session_data_dir_abs / "headers.gff3"
    feature_file = session_data_dir_abs / "features.gff3"
    temp_gff_uncompressed_filepath = session_data_dir_abs / "full.gff3"
    temp_gff_sorted_filepath = session_data_dir_abs / "sorted.gff3"
    
    with open(header_file, 'w', encoding='utf-8', newline='') as f:
        f.write('\n'.join(header_lines) + '\n')
    
    with open(feature_file, 'w', encoding='utf-8', newline='') as f:
        f.write('\n'.join(cleaned_feature_lines) + '\n')
    
    try:
        env = os.environ.copy()
        env['LC_ALL'] = 'C'

        sort_command = [
            SORT_PATH,
            '-k1,1',
            '-k4,4n',
            str(feature_file)
        ]
        
        print(f"Debug: Running sort command: {' '.join(sort_command)}")
        
        with open(temp_gff_sorted_filepath, 'wb') as outfile:
            process = subprocess.run(
                sort_command,
                check=True,
                stdout=outfile,
                stderr=subprocess.PIPE,
                env=env
            )
            if process.stderr:
                print(f"Sort stderr: {process.stderr.decode('utf-8', errors='replace')}")
        
        with open(temp_gff_uncompressed_filepath, 'w', encoding='utf-8', newline='') as outfile:
            with open(header_file, 'r', encoding='utf-8') as f_header:
                outfile.write(f_header.read())
            with open(temp_gff_sorted_filepath, 'r', encoding='utf-8', newline='') as f_sorted:
                outfile.write(f_sorted.read())
                
        print(f"Debug: Combined headers and sorted features into {temp_gff_uncompressed_filepath}")
        
        for temp_file in [header_file, feature_file, temp_gff_sorted_filepath]:
            try:
                if temp_file.exists():
                    temp_file.unlink()
                    print(f"Debug: Removed temporary file: {temp_file}")
            except Exception as e:
                print(f"Warning: Could not remove {temp_file}: {e}")

    except FileNotFoundError:
        raise RuntimeError(f"sort executable not found at {SORT_PATH}. Please check the path configuration.")
    except subprocess.CalledProcessError as e:
        stderr_output = e.stderr.decode('utf-8', errors='replace') if e.stderr else ""
        raise RuntimeError(f"sort failed with error (exit code {e.returncode}):\n{stderr_output}")
    except Exception as e:
        traceback.print_exc()
        raise RuntimeError(f"An unexpected error occurred during GFF3 sorting: {e}")

    try:
        abs_temp_gff_uncompressed_filepath = str(temp_gff_uncompressed_filepath)
        abs_gff_gz_filepath_fs = str(gff_gz_filepath_fs)

        print("Debug: Using external bgzip.exe for compression.")
        
        bgzip_command = [BGZIP_PATH, '-c', '-@', '2', abs_temp_gff_uncompressed_filepath]
        
        temp_bgzip_output = session_data_dir_abs / "temp_bgzip_output.gz"
        
        with open(temp_bgzip_output, 'wb') as outfile:
            print(f"Debug: Running bgzip command: {' '.join(bgzip_command)}")
            process = subprocess.run(
                bgzip_command,
                check=True,
                stdout=outfile,
                stderr=subprocess.PIPE
            )
            if process.stderr:
                print(f"bgzip stderr: {process.stderr.decode('utf-8', errors='replace')}")
        
        shutil.move(str(temp_bgzip_output), abs_gff_gz_filepath_fs)
        print(f"Debug: Moved compressed GFF to {abs_gff_gz_filepath_fs}")

        try:
            with gzip.open(abs_gff_gz_filepath_fs, 'rt', encoding='utf-8') as test_f:
                first_line = test_f.readline()
                if not first_line.startswith('##gff-version 3'):
                    raise RuntimeError("Decompressed content doesn't start with GFF header")
                
            print("Debug: BGZF validation passed (via Python's gzip)")
        except Exception as e:
            if os.path.exists(abs_gff_gz_filepath_fs):
                os.remove(abs_gff_gz_filepath_fs)
            raise RuntimeError(f"BGZF validation failed (via Python's gzip): {str(e)}")
            
        try:
            if temp_gff_uncompressed_filepath.exists():
                temp_gff_uncompressed_filepath.unlink()
                print(f"Debug: Removed uncompressed GFF: {temp_gff_uncompressed_filepath}")
        except Exception as e:
            print(f"Warning: Could not remove {temp_gff_uncompressed_filepath}: {e}")

    except FileNotFoundError:
        raise RuntimeError(f"bgzip executable not found at {BGZIP_PATH}. Please check the path configuration.")
    except subprocess.CalledProcessError as e:
        stderr_output = e.stderr.decode('utf-8', errors='replace') if e.stderr else ""
        raise RuntimeError(f"bgzip failed with error (exit code {e.returncode}):\n{stderr_output}")
    except Exception as e:
        traceback.print_exc()
        raise RuntimeError(f"An unexpected error occurred during bgzip compression: {e}")

    try:
        abs_gff_gz_filepath_fs = str(gff_gz_filepath_fs)
        abs_gff_tbi_filepath_fs = str(gff_tbi_filepath_fs)

        env = os.environ.copy()
        env['LC_ALL'] = 'C'

        tabix_command = [TABIX_PATH, '-f', '-p', 'gff', abs_gff_gz_filepath_fs]
        print(f"Debug: Running tabix command: {' '.join(tabix_command)}")
        process = subprocess.run(
            tabix_command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env
        )
        print(f"Debug: tabix stdout: {process.stdout.decode().strip()}")
        print(f"Debug: tabix stderr: {process.stderr.decode().strip()}")
        
        if not gff_tbi_filepath_fs.exists() or gff_tbi_filepath_fs.stat().st_size == 0:
            error_message = f"Tabix indexing failed: .tbi file was not created or is empty ({abs_gff_tbi_filepath_fs}). "
            error_message += f"Tabix stdout: {process.stdout.decode().strip()}\nTabix stderr: {process.stderr.decode().strip()}"
            print(f"Error: {error_message}")
            raise RuntimeError(error_message)

    except FileNotFoundError:
        raise RuntimeError(f"tabix executable not found at {TABIX_PATH}. Please check the path configuration.")
    except subprocess.CalledProcessError as e:
        stdout_output = e.stdout.decode('utf-8', errors='replace') if e.stdout else ""
        stderr_output = e.stderr.decode('utf-8', errors='replace') if e.stderr else ""
        raise RuntimeError(f"tabix failed with error (exit code {e.returncode}):\nSTDOUT: {stdout_output}\nSTDERR: {stderr_output}")
    except Exception as e:
        traceback.print_exc()
        raise RuntimeError(f"An unexpected error occurred during tabix indexing: {e}")

    time.sleep(1.0)

    base_web_url_for_session_data = f"/static/jbrowse2/data/{session_id}/"

    config_data = {
        "assemblies": [
            {
                "name": seq_record.id,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": "reference",
                    "adapter": {
                        "type": "IndexedFastaAdapter",
                        "fastaLocation": { "uri": f"{base_web_url_for_session_data}{fasta_data_filename}" },
                        "faiLocation": { "uri": f"{base_web_url_for_session_data}{fasta_data_filename}.fai" }
                    }
                }
            }
        ],
        "tracks": [
            {
                "type": "FeatureTrack",
                "trackId": "annotations_track",
                "name": "Genome Annotations",
                "assemblyNames": [seq_record.id],
                "adapter": {
                    "type": "Gff3TabixAdapter",
                    "gffGzLocation": { "uri": f"{base_web_url_for_session_data}{gff_data_filename_base}.gz" },
                    "index": {
                        "location": { "uri": f"{base_web_url_for_session_data}{gff_data_filename_base}.gz.tbi" }
                    }
                },
                "displays": [
                    {
                        "type": "LinearBasicDisplay",
                        "displayId": "annotations_track_display",
                        "renderer": {
                            "type": "SvgFeatureRenderer",
                            "color": "function(feature) { return feature.get('type') === 'gene' ? '#3366cc' : (feature.get('type') === 'CDS' ? '#ff9900' : '#888888'); }",
                            "showLabels": True,
                            "showDescriptions": True,
                            "label": "function(feature) { return feature.get('gene') || feature.get('Name') || feature.get('ID') || feature.get('type'); }",
                            "description": "function(feature) { return feature.get('product') || feature.get('Note') || feature.get('description'); }"
                        }
                    }
                ]
            }
        ],
        "defaultSession": {
            "view": {
                "type": "LinearGenomeView",
                "tracks": ["reference", "annotations_track"]
            }
        }
    }
    
    config_filepath_fs = session_data_dir_abs / 'config.json'
    with open(config_filepath_fs, "w") as f:
        json.dump(config_data, f, indent=2)
    
    return session_id

def cleanup_jbrowse(session_id):
    jbrowse_data_dir = Path(app.static_folder) / 'jbrowse2' / 'data' / session_id
    if jbrowse_data_dir.exists():
        shutil.rmtree(jbrowse_data_dir, ignore_errors=True)

def write_gff3(seq_record):
    gff_lines = []
    
    for feature in seq_record.features:
        if feature.type.lower() == 'region':
            continue

        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand_str = '+' if feature.location.strand == 1 else ('-' if feature.location.strand == -1 else '.')
        phase = feature.qualifiers.get('codon_start', ['0'])[0] if feature.type == 'CDS' else '.'

        feature_id = feature.id if feature.id and feature.id != '<unknown id>' else feature.qualifiers.get('ID', [''])[0]
        if not feature_id:
            feature_id = f"{feature.type}-{start}-{end}-{str(uuid.uuid4())[:8]}"
            feature.id = feature_id
        feature.qualifiers['ID'] = [feature_id]

        fields = [
            seq_record.id,
            'GenomeAnnotationTool',
            feature.type,
            str(start),
            str(end),
            '.',
            strand_str,
            phase
        ]
        
        attributes = []
        for key, values in feature.qualifiers.items():
            if key == 'codon_start' and feature.type == 'CDS':
                continue
            
            processed_values_for_key = []
            for val_item in values:
                cleaned_val = str(val_item).strip()
                escaped_val = (
                    cleaned_val.replace('%', '%25')
                               .replace(';', '%3B')
                               .replace('=', '%3D')
                               .replace(',', '%2C')
                               .replace('\n', '%0A')
                               .replace('\r', '%0D')
                )
                processed_values_for_key.append(escaped_val)

            attributes.append(f"{key}={','.join(processed_values_for_key)}")
        
        attributes.sort()
        
        gff_line = '\t'.join(fields) + '\t' + ';'.join(attributes)
        gff_lines.append(gff_line)
    
    return '\n'.join(gff_lines) + '\n'

def parse_files(gff_file_stream, fasta_file_stream):
    fasta_content = fasta_file_stream.read().decode('utf-8')
    try:
        seq_record = SeqIO.read(StringIO(fasta_content), 'fasta')
    except Exception as e:
        raise ValueError(f"FASTA parsing error: {str(e)}")
    
    seq_length = len(seq_record.seq)
    if seq_length == 0:
        raise ValueError("Empty sequence in FASTA file")
    
    gff_content = gff_file_stream.read().decode('utf-8')
    
    seq_record.features = []
    for line_num, line in enumerate(gff_content.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        parts = line.split('\t')
        if len(parts) < 9:
            print(f"Warning: Skipping malformed GFF3 line {line_num} (not enough columns): {line}")
            continue
            
        try:
            seqid, source, type_, start, end, score, strand_str, phase_str, attributes_str = parts
            
            start_pos = int(start) - 1
            end_pos = int(end)
            
            if start_pos < 0 or end_pos > seq_length or start_pos >= end_pos:
                raise ValueError(f"Invalid coordinates ({start}-{end}) or outside sequence length ({seq_length}) at line {line_num}")
                
            strand_val = 0
            if strand_str == '+':
                strand_val = 1
            elif strand_str == '-':
                strand_val = -1
            elif strand_str != '.':
                raise ValueError(f"Invalid strand ('{strand_str}') at line {line_num}")
            
            attr_dict = {}
            for attr_pair in attributes_str.split(';'):
                if '=' in attr_pair:
                    key, val = attr_pair.split('=', 1)
                    val = val.replace('%2C', ',').replace('%3D', '=').replace('%3B', ';').replace('%25', '%')\
                             .replace('%0A', '\n').replace('%0D', '\r')
                    if key.lower() == 'parent' or key.lower() == 'dbxref':
                        attr_dict.setdefault(key, []).extend(val.split(','))
                    else:
                        attr_dict.setdefault(key, []).append(val)
            
            codon_start = '1'
            if type_ == 'CDS' and phase_str != '.':
                codon_start = str(int(phase_str) + 1)
            if codon_start != '1':
                attr_dict.setdefault('codon_start', []).append(codon_start)

            feature = SeqFeature(
                FeatureLocation(start_pos, end_pos, strand=strand_val),
                type=type_,
                qualifiers=attr_dict
            )
            
            if 'ID' in attr_dict:
                feature.id = attr_dict['ID'][0]

            seq_record.features.append(feature)
            
        except ValueError as ve:
            print(f"Error parsing line {line_num} due to data validation: {ve}. Line: {line}")
        except Exception as e:
            print(f"General error parsing line {line_num}: {e}. Line: {line}")
    
    return seq_record


def feature_to_dict(feature):
    qualifiers = feature.qualifiers.copy()
    
    gene_name = qualifiers.get('gene', [''])[0] or \
                qualifiers.get('Name', [''])[0] or \
                qualifiers.get('gene_name', [''])[0] or \
                qualifiers.get('locus_tag', [''])[0]
    
    product_name = qualifiers.get('product', [''])[0] or \
                   qualifiers.get('description', [''])[0] or \
                   qualifiers.get('Note', [''])[0] or \
                   qualifiers.get('function', [''])[0]
    
    codon_start_val = qualifiers.get('codon_start', ['1'])[0]

    exons_list = []
    if hasattr(feature.location, 'parts') and feature.location.parts:
        for loc_part in sorted(feature.location.parts, key=lambda x: x.start):
            exons_list.append({
                'start': int(loc_part.start) + 1,
                'end': int(loc_part.end),
                'strand': int(loc_part.strand)
            })
    else:
        exons_list.append({
            'start': int(feature.location.start) + 1,
            'end': int(feature.location.end),
            'strand': int(feature.location.strand)
        })

    feature_id_for_dict = feature.id if feature.id and feature.id != '<unknown id>' else qualifiers.get('ID', [''])[0]
    if not feature_id_for_dict:
        feature_id_for_dict = str(uuid.uuid4())

    return {
        'id': feature_id_for_dict,
        'type': feature.type,
        'start': int(feature.location.start) + 1,
        'end': int(feature.location.end),
        'strand': int(feature.location.strand),
        'gene': gene_name,
        'product': product_name,
        'note': qualifiers.get('Note', [''])[0],
        'codon_start': codon_start_val,
        'exons': exons_list
    }

def dict_to_feature(feature_dict, seq_length):
    exon_locations = []
    for exon in feature_dict.get('exons', []):
        start = max(0, exon['start'] - 1)
        end = min(seq_length, exon['end'])
        strand = exon['strand']
        exon_locations.append(FeatureLocation(start, end, strand=strand))
    
    if len(exon_locations) > 1:
        location = CompoundLocation(sorted(exon_locations, key=lambda x: x.start))
    else:
        location = exon_locations[0] if exon_locations else FeatureLocation(0, 0, strand=1)
    
    qualifiers = {
        'gene': [feature_dict['gene']] if feature_dict.get('gene') else [],
        'product': [feature_dict['product']] if feature_dict.get('product') else [],
        'Note': [feature_dict['note']] if feature_dict.get('note') else []
    }
    
    if feature_dict['type'] == 'CDS' and feature_dict.get('codon_start'):
        qualifiers['codon_start'] = [str(feature_dict['codon_start'])]
    
    if feature_dict.get('id'):
        qualifiers.setdefault('ID', []).append(feature_dict['id'])

    feature = SeqFeature(
        location=location,
        type=feature_dict['type'],
        qualifiers=qualifiers
    )
    
    if feature_dict.get('id'):
        feature.id = feature_dict['id']
    
    return feature

def reverse_complement_features(seq_record):
    if len(seq_record.seq) == 0:
        raise ValueError("Cannot reverse complement empty sequence")
    
    seq_record.seq = seq_record.seq.reverse_complement()
    seq_length = len(seq_record.seq)

    new_features = []
    for feature in seq_record.features:
        locations_to_process = []
        if hasattr(feature.location, 'parts'):
            locations_to_process.extend(feature.location.parts)
        else:
            locations_to_process.append(feature.location)

        new_parts = []
        for loc in locations_to_process:
            new_start = seq_length - loc.end
            new_end = seq_length - loc.start
            new_strand = -1 * loc.strand
            new_parts.append(FeatureLocation(new_start, new_end, strand=new_strand))
        
        new_parts.sort(key=lambda x: x.start)

        if len(new_parts) > 1:
            feature.location = CompoundLocation(new_parts)
        else:
            feature.location = new_parts[0]
        
        if feature.type == 'CDS' and 'codon_start' in feature.qualifiers:
            feature.qualifiers['codon_start'] = ['1']

        new_features.append(feature)
    
    seq_record.features = new_features
    return seq_record

def shift_origin(seq_record, new_start):
    if len(seq_record.seq) == 0:
        raise ValueError("Cannot shift origin of empty sequence")
    
    seq_length = len(seq_record.seq)
    new_start_0based = new_start - 1
    
    if new_start_0based < 0 or new_start_0based >= seq_length:
        raise ValueError(f"Origin position must be between 1 and {seq_length}")
    
    shifted_seq = seq_record.seq[new_start_0based:] + seq_record.seq[:new_start_0based]
    seq_record.seq = shifted_seq
    
    new_features = []
    for feature in seq_record.features:
        locations_to_process = []
        if hasattr(feature.location, 'parts'):
            locations_to_process.extend(feature.location.parts)
        else:
            locations_to_process.append(feature.location)

        new_parts = []
        for loc in locations_to_process:
            exon_start = loc.start
            exon_end = loc.end
            
            new_exon_start = (exon_start - new_start_0based + seq_length) % seq_length
            new_exon_end = (exon_end - new_start_0based + seq_length) % seq_length
            
            new_parts.append(FeatureLocation(new_exon_start, new_exon_end, strand=loc.strand))
        
        new_parts.sort(key=lambda x: x.start)

        if len(new_parts) > 1:
            feature.location = CompoundLocation(new_parts)
        else:
            feature.location = new_parts[0]
        
        new_features.append(feature)

    seq_record.features = new_features
    return seq_record

def validate_cds(feature, seq_record, codon_table_id):
    if feature.type.lower() != 'cds':
        return {'valid': True, 'message': 'Not a CDS feature'}
    
    codon_table = next(
        (t for t in CODON_TABLES if t['id'] == codon_table_id),
        None
    )
    if not codon_table:
        return {'valid': False, 'message': 'Invalid codon table selected'}
    
    cds_seq_obj = Seq("")
    all_locations = []
    if hasattr(feature.location, 'parts'):
        all_locations.extend(feature.location.parts)
    else:
        all_locations.append(feature.location)
    
    all_locations.sort(key=lambda x: x.start)

    for exon_loc in all_locations:
        exon_seq = seq_record.seq[exon_loc.start:exon_loc.end]
        if feature.location.strand == -1:
            exon_seq = exon_seq.reverse_complement()
        cds_seq_obj += exon_seq
    
    cds_seq = str(cds_seq_obj).upper()
    
    codon_start_val = int(feature.qualifiers.get('codon_start', ['1'])[0])
    cds_seq_framed = cds_seq[codon_start_val - 1:]

    errors = []
    
    if len(cds_seq_framed) % 3 != 0:
        errors.append(f"CDS length ({len(cds_seq_framed)}) is not divisible by 3 after applying codon_start.")
    
    if len(cds_seq_framed) >= 3:
        start_codon = cds_seq_framed[:3]
        if start_codon not in codon_table['start_codons']:
            errors.append(f"Invalid start codon: {start_codon}")
        
        stop_codon = cds_seq_framed[-3:]
        if stop_codon not in codon_table['stop_codons']:
            errors.append(f"Invalid stop codon: {stop_codon}")
        
        for i in range(3, len(cds_seq_framed) - 3, 3):
            codon = cds_seq_framed[i:i+3]
            if codon in codon_table['stop_codons']:
                errors.append(f"Internal stop codon found at relative position {i+1} ({codon}).")
    elif len(cds_seq_framed) > 0:
        errors.append(f"CDS length ({len(cds_seq_framed)}) too short for codon analysis.")
        
    return {
        'valid': len(errors) == 0,
        'message': 'CDS is valid' if not errors else '; '.join(errors)
    }

@app.route('/')
def index():
    return send_from_directory('static', 'index.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    if 'gff_file' not in request.files or 'fasta_file' not in request.files:
        return jsonify({'error': 'Missing files'}), 400
    
    try:
        gff_file = request.files['gff_file']
        fasta_file = request.files['fasta_file']
        
        seq_record = parse_files(gff_file, fasta_file)
        gff_feature_lines_content = write_gff3(seq_record)
        
        session_id = str(uuid.uuid4())
        
        setup_jbrowse_data_with_external_tools(seq_record, gff_feature_lines_content, session_id)
        
        return jsonify({
            'sequence_id': seq_record.id,
            'sequence_length': len(seq_record.seq),
            'features': [feature_to_dict(f) for f in seq_record.features],
            'session_id': session_id
        })
    except Exception as e:
        print(f"Server error during upload: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/jbrowse/')
def jbrowse_viewer():
    session_id = request.args.get('session_id')
    if not session_id:
        return "Session ID required", 400
    
    return send_from_directory(os.path.join(app.static_folder, 'jbrowse2'), 'viewer.html')

@app.route('/static/jbrowse2/data/<path:filename>')
def jbrowse_data(filename):
    base_data_dir = os.path.join(app.static_folder, 'jbrowse2', 'data')
    
    response = send_from_directory(base_data_dir, filename)
    
    if filename.endswith('.gz'):
        response.headers['Content-Type'] = 'application/x-gzip'
    elif filename.endswith('.tbi'):
        response.headers['Content-Type'] = 'application/octet-stream'
        
    return response

@app.route('/validate_cds', methods=['POST'])
def validate_cds_feature():
    try:
        data = request.get_json()
        if not data or 'feature' not in data or 'codon_table' not in data:
            return jsonify({'error': 'Invalid request data'}), 400
        
        feature_dict = data['feature']
        codon_table_id = data['codon_table']
        sequence_length = data.get('sequence_length', 0)

        seq_record_for_validation = SeqRecord(Seq('N' * sequence_length), id='temp_validation_seq')

        feature_obj = dict_to_feature(feature_dict, sequence_length)
        result = validate_cds(feature_obj, seq_record_for_validation, codon_table_id)
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/flip', methods=['POST'])
def flip_sequence():
    try:
        data = request.get_json()
        if not data or 'sequence_length' not in data or 'features' not in data:
            return jsonify({'error': 'Invalid request data'}), 400
        
        seq_record = SeqRecord(Seq('N' * data['sequence_length']), id=data.get('sequence_id', 'unknown'))
        seq_record.features = [dict_to_feature(f, data['sequence_length']) for f in data['features']]
        
        seq_record = reverse_complement_features(seq_record)
        
        features = [feature_to_dict(f) for f in seq_record.features]
        
        return jsonify({
            'sequence_id': seq_record.id,
            'sequence_length': len(seq_record.seq),
            'features': features
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/shift_origin', methods=['POST'])
def shift_origin_route():
    try:
        data = request.get_json()
        if not data or 'sequence_length' not in data or 'features' not in data or 'new_start' not in data:
            return jsonify({'error': 'Invalid request data'}), 400
        
        seq_record = SeqRecord(Seq('N' * data['sequence_length']), id=data.get('sequence_id', 'unknown'))
        seq_record.features = [dict_to_feature(f, data['sequence_length']) for f in data['features']]
        
        new_start = int(data['new_start'])
        seq_record = shift_origin(seq_record, new_start)
        
        features = [feature_to_dict(f) for f in seq_record.features]
        
        return jsonify({
            'sequence_id': seq_record.id,
            'sequence_length': len(seq_record.seq),
            'features': features
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/export', methods=['POST'])
def export_files():
    try:
        data = request.get_json()
        if not data or 'sequence_length' not in data or 'features' not in data:
            return jsonify({'error': 'Invalid request data'}), 400
        
        seq_record = SeqRecord(
            Seq('N' * data['sequence_length']),
            id=data.get('sequence_id', 'unknown'),
            features=[dict_to_feature(f, data['sequence_length']) for f in data['features']]
        )
        
        gff_output_features_only = write_gff3(seq_record)
        fasta_output = StringIO()
        SeqIO.write(seq_record, fasta_output, 'fasta')
        
        full_gff_export = "##gff-version 3\n"
        full_gff_export += f"##sequence-region {seq_record.id} 1 {len(seq_record.seq)}\n"
        full_gff_export += gff_output_features_only
        
        return jsonify({
            'gff': full_gff_export,
            'fasta': fasta_output.getvalue()
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/cleanup', methods=['POST'])
def cleanup_session():
    try:
        data = request.get_json()
        if not data or 'session_id' not in data:
            return jsonify({'error': 'Session ID required'}), 400
            
        cleanup_jbrowse(data['session_id'])
        return jsonify({'success': True})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)

