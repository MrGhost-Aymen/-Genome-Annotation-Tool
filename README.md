# 🧬 Genome Annotation Tool

A web-based interactive tool for visualizing, editing, and exporting genomic annotations using FASTA and GFF3 files.

![Screenshot](![image](https://github.com/user-attachments/assets/13da2c84-f6a7-4d01-bd13-8e6d0a78d0d7)
) 

---

## 📌 Features

- **Visualize genomes** with IGV.js
- **Add/Edit/Delete annotations** (genes, CDS, exons, etc.)
- **Search & Filter features** by type, gene name, product, or ID
- **Export annotated sequences** in FASTA and GFF3 formats
- **Shift origin** or **flip sequence** to reorient the genome
- **Session management** for multi-user support
- **Responsive design** for desktop and mobile

---

## 🧰 Technologies Used

- **Frontend**: HTML5, Tailwind CSS, JavaScript, IGV.js
- **Backend**: Python Flask
- **Data Handling**: Biopython for parsing FASTA/GFF3
- **File Compression**: Optional bgzip/tabix support (configurable)

---

## 🚀 Getting Started

### Prerequisites

- Python 3.8+
- Node.js (for Tailwind CLI, optional)
- [Optional] `bgzip`, `tabix`, and `sort` binaries installed for compressed output

### Clone the Repository

```bash
git clone https://github.com/MrGhost-Aymen/Genome-Annotation-Tool.git
cd genome-annotation-tool
```

### Install Dependencies

#### Backend (Python Flask)

```bash
pip install -r requirements.txt
```

---

## ⚙️ Configuration

Edit the top of `app.py` to set paths for external tools:

```python
BGZIP_PATH = "/usr/bin/bgzip"  # Update this path
TABIX_PATH = "/usr/bin/tabix"
SORT_PATH = "/usr/bin/sort"

USE_UNCOMPRESSED_GFF = True  # Set to False if using bgzip/tabix
```

---

## 🏃‍♂️ Running the App

Start the Flask server:

```bash
python app.py
```

Open your browser and go to:  
👉 http://localhost:5000

---

## 📁 File Structure

```
.
├── app.py                  # Main Flask backend
├── static/
│   ├── igvjs/
├── templates/
│   └── index.html          # Main frontend page
└── uploads/                # Temporary upload storage
```

---

## 🧪 Testing

You can test the application using sample data:

- FASTA file: `test_data/sample.fasta`
- GFF3 file: `test_data/sample.gff3`

---

## 💾 Exporting Annotations

After making changes:
1. Click **Export GFF3 & FASTA**
2. Your browser will download both updated files

---

## 🧹 Session Cleanup

Each user gets a unique session. You can manually clean up sessions via:

```bash
POST /cleanup { "session_id": "your-session-id" }
```

Or let sessions expire naturally.

---

## 🛡️ Security Notes

- Input sanitization is important when accepting GFF3 files from users.
- Consider CSRF protection for production deployment.
- Validate all feature edits before writing back to GFF3.

---

## 🧩 Future Enhancements

- Add Undo/Redo functionality
- Support multiple tracks and layered annotations
- Implement user login and persistent storage
- Integrate BLAST search or functional annotation suggestions
- Add unit tests and CI/CD pipeline

---

## 🤝 Contributing

Contributions are welcome! Please open an issue or submit a pull request.


---

## 📬 Questions?

For questions or feedback, feel free to [open an issue](ouamoa@gmail.com).

---

Let me know if you'd like me to generate a ZIP-ready version or help with GitHub Actions CI setup as well!
