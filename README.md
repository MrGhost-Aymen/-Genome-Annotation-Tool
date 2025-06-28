Genome Annotation Tool (Flask + JBrowse 2)

This Flask web application provides a user-friendly interface for uploading, processing, visualizing, and managing genome annotation files (FASTA and GFF3). It integrates external bioinformatics tools like sort, bgzip, and tabix to prepare data for interactive visualization using a bundled JBrowse 2 viewer.
Features

    File Upload: Securely upload FASTA sequence files and GFF3 annotation files.

    Automated Processing:

        Parses FASTA and GFF3 into Biopython SeqRecord objects.

        Combines GFF3 headers with sorted feature lines.

        Compresses GFF3 into BGZF format using bgzip.

        Indexes BGZF GFF3 with tabix for fast random access.

        Generates JBrowse 2-compatible config.json, FASTA, and FASTA index (.fai) files.

    Interactive Visualization: Seamlessly launches a JBrowse 2 instance to display your genome sequence and annotations.

    Annotation Management:

        View parsed annotations in a table.

        Validate CDS (Coding Sequence) features against selected codon tables for common issues (length divisibility by 3, start/stop codons, internal stop codons).

        Perform genome modifications:

            Reverse Complement the entire sequence and adjust feature coordinates.

            Shift the circular genome origin and update feature locations.

    Data Export: Download the modified FASTA and GFF3 files.

    Session Management: Creates unique sessions for each upload, allowing isolated data processing.

Prerequisites

Before running the application, ensure you have the following installed:

    Python 3.x: Download from python.org.

    Flask: A Python web framework.

    Biopython: For sequence and feature manipulation.

    External HTSlib Tools: bgzip and tabix (from HTSlib)

    External sort utility: A robust sorting tool (e.g., from GNU CoreUtils or MSYS2).

Acquiring HTSlib Tools and sort on Windows (Important!)

This application relies on external command-line tools (bgzip.exe, tabix.exe, sort.exe) to process genomic files. For Windows users, acquiring compatible versions that work well with JBrowse 2's client-side JavaScript decompressor can be tricky.

Recommended Approaches for Windows:

    1. Windows Subsystem for Linux (WSL2) - Highly Recommended for Reliability:
    This is the most robust solution for Windows users, providing a native Linux environment where HTSlib tools behave predictably.

        Install WSL2: Follow Microsoft's guide: Install WSL.

        Install a Linux Distribution: e.g., Ubuntu from the Microsoft Store.

        Inside WSL2, install samtools and sort:

        sudo apt update
        sudo apt install samtools # This includes bgzip and tabix
        sudo apt install coreutils # This includes sort

        Run Flask from WSL2: You can then clone this repository and run your Flask app entirely within your WSL2 Linux environment. Flask will be accessible from your Windows browser via localhost.

    2. MSYS2 (Recommended for Native Windows Binaries):
    MSYS2 provides a Unix-like environment for Windows, allowing you to install GNU tools like sort and bioinformatics tools like HTSlib (which includes bgzip and tabix).

        Download and install MSYS2: From msys2.org.

        Open the MSYS2 MinGW 64-bit terminal.

        Install necessary packages:

        pacman -Syu
        pacman -Su
        pacman -S mingw-w64-x86_64-htslib # For bgzip and tabix
        pacman -S msys/coreutils         # For sort

        Locate the executables: After installation, you will typically find them in C:\msys64\mingw64\bin\ (for bgzip.exe, tabix.exe) and C:\msys64\usr\bin\ (for sort.exe).

        Update paths in app.py: You must update the BGZIP_PATH, TABIX_PATH, and SORT_PATH variables in app.py to point to these exact locations on your system.

        BGZIP_PATH = r"C:\msys64\mingw64\bin\bgzip.exe"
        TABIX_PATH = r"C:\msys64\mingw64\bin\tabix.exe"
        SORT_PATH = r"C:\msys64\usr\bin\sort.exe"

    3. Pre-compiled Binaries (Less Reliable for JBrowse Compatibility):
    You can try downloading pre-compiled Windows binaries for bgzip, tabix, and sort from various bioinformatics tool websites. However, we have observed that certain Windows builds of bgzip.exe can produce BGZF files that JBrowse 2's JavaScript decompressor might find problematic, leading to "incorrect gzip header check" errors, even if bgzip.exe itself can decompress them. Using the MSYS2 or WSL2 builds is generally more reliable.

Setup Instructions

    Clone the repository:

    git clone https://github.com/MrGhost-Aymen/genome-annotation-tool.git
    cd genome-annotation-tool

    Create a virtual environment (recommended):

    python -m venv venv
    # On Windows:
    .\venv\Scripts\activate
    # On macOS/Linux:
    source venv/bin/activate

    Install Python dependencies:

    pip install -r requirements.txt

    (Create a requirements.txt file with Flask, biopython, werkzeug, pathlib, python-mimetypes).
    Example requirements.txt:

    Flask
    biopython
    werkzeug
    # pathlib is usually built-in
    # mimetypes is usually built-in

    Acquire and configure external tools:
    Follow the "Acquiring HTSlib Tools on Windows" section above to get bgzip.exe, tabix.exe, and sort.exe and update their paths in app.py.

    Download and setup JBrowse 2 static files:

        Go to the JBrowse 2 website.

        Download the latest "JBrowse 2 Web" distribution (look for the .zip or .tar.gz file).

        Extract the contents of the downloaded archive.

        Rename the extracted folder (e.g., jbrowse-web-2.x.x) to jbrowse2.

        Move the jbrowse2 folder into your Flask application's static directory.
        Your project structure should look something like this:

        genome-annotation-tool/
        ├── app.py
        ├── static/
        │   ├── jbrowse2/          <-- JBrowse 2 extracted here
        │   │   ├── index.html
        │   │   ├── config.json
        │   │   ├── data/
        │   │   └── ... (other JBrowse files)
        │   ├── codon_tables.json
        │   └── bgzip.exe          <-- Your bgzip.exe
        │   └── tabix.exe          <-- Your tabix.exe
        ├── requirements.txt
        └── venv/
        └── uploads/

        Ensure static/jbrowse2/viewer.html exists and is the main entry point for the JBrowse viewer (as referenced in app.py).

    Run the Flask application:

    python app.py

    The application will typically run on http://127.0.0.1:5000.

Usage

    Open your browser and navigate to http://127.0.0.1:5000.

    Upload your FASTA and GFF3 files using the provided forms.

    Upon successful upload and processing, a new JBrowse 2 viewer will open, displaying your genomic data.

    Use the application's interface to validate CDS features, flip the sequence, shift the origin, or export your processed files.

Troubleshooting

    "Error: problem decompressing block: incorrect gzip header check" in JBrowse 2:
    This is a known issue specific to certain Windows builds of bgzip.exe and JBrowse 2's JavaScript decompressor. Even if bgzip.exe and tabix.exe appear to work correctly on the server-side, JBrowse's client-side JavaScript can be very sensitive to the exact byte structure of the BGZF file.

        Solution: Follow the WSL2 setup instructions in the Prerequisites section. Running the Flask app and HTSlib tools within WSL2's Linux environment provides binaries that are highly compatible with JBrowse 2.

        Alternatively, try different builds of bgzip.exe (e.g., from newer MSYS2 installations, or other pre-compiled sources), but WSL2 is the most reliable fix.

    sort or bgzip or tabix FileNotFoundError or CalledProcessError:

        Check Paths: Ensure BGZIP_PATH, TABIX_PATH, and SORT_PATH in app.py are absolute and correct for your system.

        Permissions: Make sure the Flask app has execute permissions on the .exe files.

        Dependencies: Ensure the msys2 environment (if used) is correctly set up and all required DLLs are in place or on your system's PATH. Running from the MSYS2 MinGW 64-bit terminal often helps.

    "Empty sequence in FASTA file" or "FASTA parsing error":

        Ensure your FASTA file is correctly formatted and contains a sequence.

    "Invalid coordinates" or "Malformed GFF3 line":

        Review your GFF3 file for syntax errors, incorrect 1-based vs 0-based coordinates, or features extending beyond the sequence length.

Contributing

Feel free to fork this repository, open issues, and submit pull requests.
