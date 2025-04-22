# Circos Input File Generator

This Python script generates input files (`karyotype.txt`, `labels.txt`, and `links.txt`) required to create Circos visualizations from GenBank and BLAST data. The script is designed to work with two genomes (e.g., plastome and mitochondrion) and their corresponding BLAST alignment report.

## Features

- **Customizable Genome Names**: Specify custom names for query and subject genomes via command-line arguments.
- **Strand Awareness**: Handles gene strand orientation (`+` or `-`) to ensure correct visualization in Circos.
- **Alignment Filtering**: Excludes short alignments based on a user-defined minimum alignment length.
- **Automated Workflow**: Generates all required input files for Circos in one step.

## Prerequisites

1. **Python 3.x**:
   - Install Python 3.x from [python.org](https://www.python.org/).

2. **Biopython**:
   - Install Biopython using pip:
     ```bash
     pip install biopython
     ```

3. **Circos**:
   - Download and install Circos from [http://circos.ca/](http://circos.ca/).
   - Ensure the `circos` executable is in your system's PATH.

4. **Input Files**:
   - Two GenBank files (one for each genome).
   - A BLAST report file in 12-column tabular format.

## Usage

### Step 1: Clone the Repository
Clone this repository to your local machine:

### Step 2: Prepare Input Files
Ensure you have the following files ready:
- Query GenBank file (e.g., `query.gbk`)
- Subject GenBank file (e.g., `subject.gbk`)
- BLAST report file (e.g., `blast_report.txt`)

The BLAST report should be in 12-column tabular format:
```
query_id    subject_id    identity    alignment_length    mismatches    gap_opens    q_start    q_end    s_start    s_end    evalue    bit_score
```

### Step 3: Run the Script
Run the script using the following command:
```bash
python circos_input_generator.py \
    --query_gbk query.gbk \
    --subject_gbk subject.gbk \
    --blast_report blast_report.txt \
    --query_name plastome \
    --subject_name mitochondrion \
    --min_alignment_length 50 \
    --output_dir ./output
```

### Command-Line Arguments
| Argument               | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `--query_gbk`          | Path to the query GenBank file (e.g., plastome).                           |
| `--subject_gbk`        | Path to the subject GenBank file (e.g., mitochondrion).                    |
| `--blast_report`       | Path to the BLAST report file (12-column tabular format).                  |
| `--query_name`         | Custom name for the query genome (e.g., `plastome`).                      |
| `--subject_name`       | Custom name for the subject genome (e.g., `mitochondrion`).               |
| `--min_alignment_length` | Minimum alignment length to include in links (default: 50).             |
| `--output_dir`         | Directory to save the output files (default: current directory).          |

### Example Output
The script generates the following files in the specified output directory:
- `karyotype.txt`: Defines the chromosomes for Circos.
- `labels.txt`: Contains gene labels for visualization.
- `links.txt`: Defines the links between query and subject genomes.

Example debugging output:
```
Parsed 120 genes from query GenBank file.
Parsed 150 genes from subject GenBank file.
Parsed 6 alignments from BLAST report.
Generated 4 links after filtering.
Circos input files generated successfully!
Karyotype file: ./output/karyotype.txt
Labels file: ./output/labels.txt
Links file: ./output/links.txt
```

## Visualizing with Circos

Once the input files are generated, you can use Circos to create the visualization:

1. Create a Circos configuration file (`circos.conf`) based on your preferences. Refer to the [Circos documentation](http://circos.ca/documentation/) for details.
2. Run Circos:
   ```bash
   circos -conf circos.conf
   ```
3. The output image (e.g., `circos.png`) will be generated in the working directory.

## Troubleshooting

- **Error: Circos executable not found**:
  - Ensure Circos is installed and added to your system's PATH.
- **Empty Links File**:
  - Check if the BLAST report contains high-quality alignments (e.g., identity ≥ 50%).
  - Verify that gene names in the BLAST report match those in the GenBank files.
- **Genes Displayed Incorrectly**:
  - Ensure strand orientation (`+` or `-`) is correctly handled in the GenBank files.

## Contributing

Contributions are welcome! If you encounter issues or have suggestions for improvement, please open an issue or submit a pull request.
