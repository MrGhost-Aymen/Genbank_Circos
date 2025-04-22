from Bio import SeqIO
import csv
import argparse
import os

def parse_genbank(gbk_file):
    """Parse GenBank file and extract gene names and positions."""
    genes = {}
    contig_length = 0
    for record in SeqIO.parse(gbk_file, "genbank"):
        contig_length = len(record.seq)  # Use the first contig's length as the genome length
        
        for feature in record.features:
            if feature.type == "CDS":
                # Use 'gene' qualifier if available, otherwise use 'locus_tag'
                gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                if gene_name == "unknown":
                    gene_name = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                
                start = int(feature.location.start)
                end = int(feature.location.end)
                
                # Use .location.strand instead of .strand and handle None values
                strand = "+" if feature.location.strand == 1 else "-"
                
                # Store gene information
                genes[gene_name] = {
                    "start": start,
                    "end": end,
                    "strand": strand
                }
    return genes, contig_length

def parse_blast_report(blast_file):
    """Parse BLAST report with 12 columns."""
    blast_data = []
    with open(blast_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            query_id = row[0]
            subject_id = row[1]
            identity = float(row[2])
            alignment_length = int(row[3])
            mismatches = int(row[4])
            gap_opens = int(row[5])
            q_start = int(row[6])
            q_end = int(row[7])
            s_start = int(row[8])
            s_end = int(row[9])
            evalue = float(row[10])
            bit_score = float(row[11])

            # Only keep high-quality matches (adjust threshold if needed)
            if identity >= 50:  # Adjusted threshold for debugging
                blast_data.append((query_id, subject_id, q_start, q_end, s_start, s_end))
    return blast_data

def generate_karyotype(query_length, subject_length, query_name, subject_name, output_file):
    """Generate karyotype file with two entries: one for query and one for subject."""
    with open(output_file, "w") as f:
        # Write query entry
        f.write(f"chr - {query_name} {query_name} 0 {query_length} green\n")

        # Write subject entry
        f.write(f"chr - {subject_name} {subject_name} 0 {subject_length} blue\n")

def generate_labels(query_genes, subject_genes, query_name, subject_name, output_file):
    """Generate labels file in the format: chr start end label."""
    with open(output_file, "w") as f:
        # Write query gene labels
        for gene_name, gene_info in query_genes.items():
            start = gene_info["start"]
            end = gene_info["end"]
            f.write(f"{query_name} {start} {end} {gene_name}\n")

        # Write subject gene labels
        for gene_name, gene_info in subject_genes.items():
            start = gene_info["start"]
            end = gene_info["end"]
            f.write(f"{subject_name} {start} {end} {gene_name}\n")

def generate_links(blast_data, query_genes, subject_genes, query_name, subject_name, output_file):
    """Generate links file."""
    with open(output_file, "w") as f:
        link_count = 0
        for query_id, subject_id, q_start, q_end, s_start, s_end in blast_data:
            # Find query gene in GenBank file
            if query_id in query_genes:
                query_gene = query_genes[query_id]
                q_genome_start = query_gene["start"] + q_start
                q_genome_end = query_gene["start"] + q_end
            else:
                print(f"Query gene {query_id} not found in query GenBank file.")
                continue

            # Find subject gene in GenBank file
            if subject_id in subject_genes:
                subject_gene = subject_genes[subject_id]
                s_genome_start = subject_gene["start"] + s_start
                s_genome_end = subject_gene["start"] + s_end
            else:
                print(f"Subject gene {subject_id} not found in subject GenBank file.")
                continue

            # Write link to file
            f.write(f"{query_name} {q_genome_start} {q_genome_end} "
                    f"{subject_name} {s_genome_start} {s_genome_end}\n")
            link_count += 1
        
        print(f"Generated {link_count} links.")

# Main script execution
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate Circos input files from BLAST report and GenBank files.")
    parser.add_argument("--query_gbk", required=True, help="Query GenBank file (e.g., plastome)")
    parser.add_argument("--subject_gbk", required=True, help="Subject GenBank file (e.g., mitochondrion)")
    parser.add_argument("--blast_report", required=True, help="BLAST report file (12-column tabular format)")
    parser.add_argument("--query_name", required=True, help="Custom name for the query genome (e.g., plastome)")
    parser.add_argument("--subject_name", required=True, help="Custom name for the subject genome (e.g., mitochondrion)")

    # Parse arguments
    args = parser.parse_args()

    # Automatically generate output filenames
    base_dir = os.path.dirname(args.query_gbk) or "."
    karyotype_file = os.path.join(base_dir, "karyotype.txt")
    labels_file = os.path.join(base_dir, "labels.txt")
    links_file = os.path.join(base_dir, "links.txt")

    # Parse GenBank files
    query_genes, query_length = parse_genbank(args.query_gbk)
    subject_genes, subject_length = parse_genbank(args.subject_gbk)

    # Print parsed genes for debugging
    print(f"Parsed {len(query_genes)} genes from query GenBank file.")
    print(f"Parsed {len(subject_genes)} genes from subject GenBank file.")

    # Parse BLAST report
    blast_data = parse_blast_report(args.blast_report)
    print(f"Parsed {len(blast_data)} alignments from BLAST report.")

    # Generate Circos input files
    generate_karyotype(query_length, subject_length, args.query_name, args.subject_name, karyotype_file)
    generate_labels(query_genes, subject_genes, args.query_name, args.subject_name, labels_file)
    generate_links(blast_data, query_genes, subject_genes, args.query_name, args.subject_name, links_file)

    print(f"Circos input files generated successfully!")
    print(f"Karyotype file: {karyotype_file}")
    print(f"Labels file: {labels_file}")
    print(f"Links file: {links_file}")