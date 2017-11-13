import sys


def extract_data(file):
    file = open(file)
    data = []
    for line in file:
        line = line.strip("\n").split(" ")
        dataline = [i for i in line if i != '']
        data.append(dataline)

    return data


def find_max_overlap(data):
    max = 0
    max_name = []
    for line in data:
        if len(line) == 7:
            if line[5] == 'seqSize':
                continue
            if int(line[5]) > max:
                max = int(line[5])
                max_name = line[1]
    return max, max_name


def find_all_unique_sequences(data, ref):
    sequences = set()
    seq_map = dict()

    for line in data:
        if len(line) == 7 and line[1] != "name" and line[1] != ref:
            if int(line[3]) < 500:
                continue
            sequences.add(line[1])
            if line[1] not in seq_map:
                seq_map[line[1]] = [[line[2], line[3], line[5]]]
            else:
                seq_map[line[1]].append([line[2], line[3], line[5]])

    return sequences, seq_map


def restrict_seq_map(seq_map, multiplicity):
    restricted_map = {}
    for name in seq_map.keys():
        if len(seq_map[name]) >= multiplicity:
            restricted_map[name] = seq_map[name]
    return restricted_map


def filter_short_reads(map, size):
    filtered = {}
    for key in map.keys():
        if int(map[key][0][2]) > size:
            filtered[key] = map[key]

    return filtered


def calculate_aln_coverage(map):
    coverage = []
    for key in map.keys():
        cov_sum = 0
        min, max = float('inf'), 0
        for hit in map[key]:
            cov_sum += int(hit[1])
            if int(hit[0]) < min:
                min = int(hit[0])
            if int(hit[0]) + int(hit[1]) > max:
                max = int(hit[0]) + int(hit[1])
        coverage.append((key, max - min, cov_sum))

    return coverage


def main(stream):
    data = extract_data(stream)
    # max, max_name = find_max_overlap(data)

    reference_name = "D38024.1"
    unique_seq, seq_map = find_all_unique_sequences(data, reference_name)
    filtered_map = filter_short_reads(seq_map, size=10000)

    coverage = calculate_aln_coverage(seq_map)

    print("Key, length between first and last alignment, total alignment length")
    print("alnStart, alnLength, readLength")
    print("")
    for vals in coverage:
        print(vals)
        print(seq_map[vals[0]])
        print("")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("ERROR: Proper usage is: " + sys.argv[0] + " file.maf")
        exit(1)

    main(sys.argv[1])
