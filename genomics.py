"""
Sree Poojitha Paruchuri's Problen Set-1
"""
import matplotlib.pyplot as plt
import argparse

def read_reads_from_file(input_file):
    """Read sequences from a file and convert them into numerical representations."""
    with open(input_file, 'r') as file:
        reads = [line.strip() for line in file.readlines()]
    numerical_reads = [[ord(nucleotide) - ord('A') for nucleotide in read] for read in reads]
    return numerical_reads

def find_overlap(read1, read2, k):
    """Find the overlap between two reads."""
    for i in range(1, min(len(read1), len(read2), k) + 1):
        if read1[-i:] == read2[:i]:
            return i
    return 0

def merge_reads(read1, read2, overlap):
    """Merge two reads based on their overlap."""
    return read1 + read2[overlap:]

def assemble_reads(input_file, output_fasta, min_contigs, max_contigs, min_contig_length, max_contig_length, overlap_parameter):
    """Assemble reads into contigs based on overlap criteria."""
    reads = read_reads_from_file(input_file)
    contigs = [reads[0]]
    total_length = len(reads[0])

    for read in reads[1:]:
        overlap = find_overlap(contigs[-1], read, overlap_parameter)
        if overlap > 0:
            contigs[-1] = merge_reads(contigs[-1], read, overlap)
            total_length += len(read) - overlap
        else:
            contigs.append(read)
            total_length += len(read)

        if total_length >= sum(map(len, reads)):
            break

    contigs = [contig for contig in contigs if min_contig_length <= len(contig) <= max_contig_length]
    if len(contigs) < min_contigs or len(contigs) > max_contigs:
        contigs = contigs[:max(max_contigs, min_contigs)]

    write_fasta(output_fasta, contigs)
    return contigs

def write_fasta(output_fasta, contigs):
    """Write the assembled contigs to a FASTA file."""
    with open(output_fasta, 'w') as file:
        for i, contig in enumerate(contigs):
            contig_sequence = ''.join(chr(nucleotide + ord('A')) for nucleotide in contig)
            file.write(f">Contig_{i + 1}\n{contig_sequence}\n")

def visualize_coverage(contigs, reads):
    """Visualize the coverage of reads across the contigs."""
    contig_length = len(contigs[0])
    coverage = [[0] * contig_length for _ in range(len(contigs))]

    for i, contig in enumerate(contigs):
        for read in reads:
            for j in range(contig_length - len(read) + 1):
                if contig[j:j + len(read)] == read:
                    for k in range(j, j + len(read)):
                        coverage[i][k] += 1

    plt.figure(figsize=(10, 6))
    for i, contig_coverage in enumerate(coverage):
        plt.plot(contig_coverage, label=f"Contig_{i + 1}")
    plt.title("Read Density Across Contigs")
    plt.xlabel("Position")
    plt.ylabel("Read Density")
    plt.legend()
    plt.show()

def main():
    """Parse command-line arguments and assemble reads into contigs."""
    parser = argparse.ArgumentParser(description="Assemble reads into contigs and visualize coverage.")
    parser.add_argument('-i', '--input_file', required=True, help="Input file containing reads")
    parser.add_argument('-o', '--output_fasta', required=True, help="Output FASTA file for assembled contigs")
    parser.add_argument('-k', '--overlap_parameter', type=int, default=5, help="Overlap parameter (k)")
    args = parser.parse_args()

    contigs = assemble_reads(args.input_file, args.output_fasta, min_contigs=10, max_contigs=20,
                             min_contig_length=100, max_contig_length=500, overlap_parameter=args.overlap_parameter)
    visualize_coverage(contigs, read_reads_from_file(args.input_file))

if __name__ == "__main__":
    main()
