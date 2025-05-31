from Bio import SeqIO
import csv
import argparse
import os
import re

def parse_genbank(gbk_file):
    """Parse GenBank file and extract gene names and positions with all possible name variants."""
    genes = {}
    contig_length = 0
    for record in SeqIO.parse(gbk_file, "genbank"):
        contig_length = len(record.seq)
        
        for feature in record.features:
            if feature.type == "CDS" or feature.type == "tRNA" or feature.type == "rRNA":
                # Get all possible name variants
                name_variants = set()
                
                # Standard qualifiers to check
                for qualifier in ["gene", "locus_tag", "product", "note"]:
                    if qualifier in feature.qualifiers:
                        for value in feature.qualifiers[qualifier]:
                            # Clean up the value
                            cleaned = value.split(';')[0].split(',')[0].strip()
                            name_variants.add(cleaned)
                
                # If no names found, use feature ID
                if not name_variants:
                    name_variants.add(feature.id)
                
                # Get coordinates
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                
                # Store all name variants pointing to the same feature
                for name in name_variants:
                    # Simplify tRNA/rRNA names by removing parentheses and other characters
                    simple_name = re.sub(r'[\(\)\-\s]', '', name)
                    genes[simple_name] = {
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "original_name": name  # Keep original for reference
                    }
    return genes, contig_length

def parse_blast_report(blast_file):
    """Parse BLAST report with more flexible handling."""
    blast_data = []
    with open(blast_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue  # Skip incomplete rows
                
            query_id = re.sub(r'[\(\)\-\s]', '', row[0])  # Simplify query name
            subject_id = re.sub(r'[\(\)\-\s]', '', row[1])  # Simplify subject name
            identity = float(row[2])
            alignment_length = int(row[3])
            
            # Only keep high-quality matches
            if identity >= 70 and alignment_length >= 30:  # Adjusted thresholds
                blast_data.append((
                    query_id, 
                    subject_id, 
                    int(row[6]),  # q_start
                    int(row[7]),  # q_end
                    int(row[8]),  # s_start
                    int(row[9])   # s_end
                ))
    return blast_data

def generate_karyotype(query_length, subject_length, query_name, subject_name, output_file):
    """Generate karyotype file."""
    with open(output_file, "w") as f:
        f.write(f"chr - {query_name} {query_name} 0 {query_length} green\n")
        f.write(f"chr - {subject_name} {subject_name} 0 {subject_length} blue\n")

def generate_labels(genes, genome_name, output_file):
    """Generate labels file for one genome."""
    with open(output_file, "w") as f:
        for gene_name, gene_info in genes.items():
            start = gene_info["start"]
            end = gene_info["end"]
            f.write(f"{genome_name} {start} {end} {gene_info['original_name']}\n")

def generate_links(blast_data, query_genes, subject_genes, query_name, subject_name, output_file):
    """Generate links file with improved name matching."""
    with open(output_file, "w") as f:
        link_count = 0
        for query_id, subject_id, q_start, q_end, s_start, s_end in blast_data:
            # Try to find the query gene
            query_gene = query_genes.get(query_id)
            if not query_gene:
                # Try alternative name formats
                for name, gene in query_genes.items():
                    if query_id in name or name in query_id:
                        query_gene = gene
                        break
            
            # Try to find the subject gene
            subject_gene = subject_genes.get(subject_id)
            if not subject_gene:
                for name, gene in subject_genes.items():
                    if subject_id in name or name in subject_id:
                        subject_gene = gene
                        break
            
            if query_gene and subject_gene:
                # Calculate genomic coordinates
                q_genome_start = query_gene["start"] + q_start
                q_genome_end = query_gene["start"] + q_end
                s_genome_start = subject_gene["start"] + s_start
                s_genome_end = subject_gene["start"] + s_end
                
                # Write link
                f.write(f"{query_name} {q_genome_start} {q_genome_end} "
                       f"{subject_name} {s_genome_start} {s_genome_end}\n")
                link_count += 1
            else:
                if not query_gene:
                    print(f"Could not find query gene: {query_id}")
                if not subject_gene:
                    print(f"Could not find subject gene: {subject_id}")
        
        print(f"Generated {link_count} links.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Circos input files from BLAST and GenBank files.")
    parser.add_argument("--query_gbk", required=True, help="Query GenBank file")
    parser.add_argument("--subject_gbk", required=True, help="Subject GenBank file")
    parser.add_argument("--blast_report", required=True, help="BLAST report file")
    parser.add_argument("--query_name", required=True, help="Name for query genome")
    parser.add_argument("--subject_name", required=True, help="Name for subject genome")
    
    args = parser.parse_args()
    
    # Output files
    base_dir = os.path.dirname(args.query_gbk) or "."
    karyotype_file = os.path.join(base_dir, "karyotype.txt")
    query_labels_file = os.path.join(base_dir, "query_labels.txt")
    subject_labels_file = os.path.join(base_dir, "subject_labels.txt")
    links_file = os.path.join(base_dir, "links.txt")
    
    # Parse files
    query_genes, query_length = parse_genbank(args.query_gbk)
    subject_genes, subject_length = parse_genbank(args.subject_gbk)
    
    print(f"Parsed {len(query_genes)} genes from query GenBank file.")
    print(f"Parsed {len(subject_genes)} genes from subject GenBank file.")
    
    blast_data = parse_blast_report(args.blast_report)
    print(f"Parsed {len(blast_data)} alignments from BLAST report.")
    
    # Generate files
    generate_karyotype(query_length, subject_length, args.query_name, args.subject_name, karyotype_file)
    generate_labels(query_genes, args.query_name, query_labels_file)
    generate_labels(subject_genes, args.subject_name, subject_labels_file)
    generate_links(blast_data, query_genes, subject_genes, args.query_name, args.subject_name, links_file)
    
    print("\nCircos input files generated:")
    print(f" - Karyotype: {karyotype_file}")
    print(f" - Query labels: {query_labels_file}")
    print(f" - Subject labels: {subject_labels_file}")
    print(f" - Links: {links_file}")
