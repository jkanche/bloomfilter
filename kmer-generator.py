import itertools

bases = ["A", "T", "G", "C"]

k = 13

kmer_list = [''.join(p) for p in itertools.product(bases, repeat=k)]

with open('kmer.txt', 'w') as f:
    for item in kmer_list:
        f.write("%s\n" % item)
