class Gff(object):
    def __init__(self, line):
        self.chrom, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attribute = line.split('\t')
        self.start = int(self.start)
        self.end = int(self.end)
    def print_gff(self):
        return f"{self.chrom}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.frame}\t{self.attribute}"
