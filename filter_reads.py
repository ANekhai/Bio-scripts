import sys
from Bio import SeqIO


def parse_fasta(file):

    reads = list(SeqIO.parse(file, "fasta"))

    return reads


def extract_aln_data(file):
    file = open(file)
    data = []
    for line in file:
        line = line.strip("\n").split(" ")
        dataline = [i for i in line if i != '']
        data.append(dataline)

    return data


def find_all_unique_sequences(data, ref):
    sequences = set()
    seq_map = dict()

    for line in data:
        if len(line) == 7 and line[1] != "name" and line[1] != ref:
            sequences.add(line[1])
            if line[1] not in seq_map:
                seq_map[line[1]] = [[line[2], line[3], line[5]]]
            else:
                seq_map[line[1]].append([line[2], line[3], line[5]])

    return sequences, seq_map


def filter_mapped_reads(reads, sequences):
    filtered = []

    for read in reads:
        if read.id in sequences:
            filtered.append(read)

    return filtered


def filter_long_reads(reads):
    filtered = [read for read in reads if len(read.seq) >= 50000]
    return filtered


def write_fasta(reads, filename):
    SeqIO.write(reads, filename, "fasta")


def main(stream):
    reads = parse_fasta(stream[1])
    # aln_data = extract_aln_data(stream[2])

    # reference_name = "D38024.1"
    # unique_seq, seq_map = find_all_unique_sequences(aln_data, reference_name)

    # filtered_reads = filter_reads(reads, unique_seq)
    filtered_reads = filter_long_reads(reads)
    write_fasta(filtered_reads, filename="reads_50k.fasta")


if __name__ == "__main__":
    '''
    if len(sys.argv) != 3:
        print("ERROR: Proper usage is: filter_reads.py reads.fasta alignment.maf")
        exit(1)
    '''

    main(sys.argv)
