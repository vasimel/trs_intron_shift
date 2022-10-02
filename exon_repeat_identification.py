import intervaltree
from intervaltree import Interval, IntervalTree
import re

#finding protein coding genes IDs and the number of their exons (key - ID, 
#value - number of exons)
coding_genes = {}
with open("/home/ubuntu/repeats_in_proteins/genes_exons.gff") as gff:
    for line in gff:
        if "protein_coding" in line.split("\t")[8] and line.split("\t")[2] == "gene":
            match = re.search("GeneID:\d+[;, ]+", line.split("\t")[8])
            if match:
                gene_id = match.group(0)
                coding_genes[gene_id] = 0
with open("/home/ubuntu/repeats_in_proteins/genes_exons.gff") as gff:
    for line in gff:
        if line.split("\t")[2] == "exon" and "pseudo=true" not in line:
            gene_id = re.search("GeneID:\d+[;, ]+", line.split("\t")[8]).group(0)
            #print(line)
            exon_number = re.search("ID=exon-\S+-\d+[;;, ]+", line.split("\t")[8]).group(0)
            #print(exon_number.split('-')[2][:-1])
            
            if gene_id in coding_genes.keys():
                try:
                    coding_genes[gene_id] = max(coding_genes[gene_id], int(exon_number.split('-')[2][:-1]))
                except:
                    continue 
               
multiexonic_list = list(filter(lambda id: coding_genes[id] > 1, coding_genes.keys())) #len = 18230

#gff with multiexonic genes
with open("/home/ubuntu/repeats_in_proteins/genes_exons.gff", "r") as genes_exons:
    with open("/home/ubuntu/repeats_in_proteins/mult_coding_genes_exons.gff", "w") as multiexonic:
        lines = genes_exons.readlines()
        for line in lines:
            match = re.search("GeneID:\d+[;, ]+", line)
            if match.group(0) in multiexonic_list:
                multiexonic.write(f"{line}")
      
         
#creating a list of genomic features intervals
interval_list=dict()
i = 0
with open("/home/ubuntu/repeats_in_proteins/mult_coding_genes_exons.gff") as gff:
    lines = gff.readlines()
    for line in lines:
        interval_list[i] = Gff(line.strip())
        i += 1

        
# remove exon entries with 0 interval length
for id, elem in interval_list.items():
    if elem.start == elem.end:
        print(f"{id} {elem.start} {elem.end}")
interval_list.pop(105867)


  
#interval tree of genomic features
exon_tree = intervalListToIntervalTree(interval_list)

#tandem repeats dictionary
trf_list=dict()
i = 0
with open("/home/ubuntu/repeats_in_proteins/trf_output.gff") as trf:
    lines = trf.readlines()
    for line in lines:
        trf_list[i] = Gff(line.strip())
        i += 1
        
        
# searching for intersections
with open("/home/ubuntu/repeats_in_proteins/intervaltree.out", "w") as output:
    for repeat in trf_list.values():
        hit = exon_tree[repeat.chrom][repeat.start:repeat.end]
        if hit:
            for entry in hit:
                f_id = entry.data
                hit_feature = interval_list[f_id]
                left_overlap = max(hit_feature.start, repeat.start)
                right_overlap = min(hit_feature.end, repeat.end)
                if hit_feature.start > repeat.start and hit_feature.end > repeat.end:
                    position = "left end"
                elif hit_feature.start < repeat.start and hit_feature.end < repeat.end:
                    position = "right end"
                elif hit_feature.start < repeat.start and hit_feature.end > repeat.end:
                    position = "repeat is inside"
                else:
                    position = "exon is inside"
                output.write(f"{hit_feature.print_gff()}\t{repeat.print_gff()}\t{position}\n")
 

#filter out plus strand exons with terminal repeats
!awk '$3 == "exon" && $7 == "+" {print $0}' /home/ubuntu/repeats_in_proteins/intervaltree.out | grep -v "repeat is inside" > /home/ubuntu/repeats_in_proteins/plus_strand_exons_with_terminal_repeats.tsv
#filter out minus strand exons with terminal repeats
!awk '$3 == "exon" && $7 == "-" {print $0}' /home/ubuntu/repeats_in_proteins/intervaltree.out | grep -v "repeat is inside" > /home/ubuntu/repeats_in_proteins/minus_strand_exons_with_terminal_repeats.tsv



#filtering plus strand exons
with open("/home/ubuntu/repeats_in_proteins/plus_strand_exons_with_terminal_repeats.tsv", "r") as exons:
    with open("/home/ubuntu/repeats_in_proteins/plus_strand_exons_with_terminal_repeats_filtered.tsv", "w") as filtered_exons:
        for line in exons:
            gene_id = re.search("GeneID:\d+[;, ]+", line.split("\t")[8]).group(0)
            exon_number = line.split('\t')[8].split(';')[0].split('-')[2]
            if exon_number == "1" and line.split('\t')[18] == "left end":
                continue
            elif exon_number == coding_genes[gene_id] and line.split('\t')[18] == "right end":
                continue
            else:
                filtered_exons.write(line)
              
            
            
            
#filtering minus strand exons
with open("/home/ubuntu/repeats_in_proteins/minus_strand_exons_with_terminal_repeats.tsv", "r") as exons:
    with open("/home/ubuntu/repeats_in_proteins/minus_strand_exons_with_terminal_repeats_filtered.tsv", "w") as filtered_exons:
        for line in exons:
            gene_id = re.search("GeneID:\d+[;, ]+", line.split("\t")[8]).group(0)
            exon_number = line.split('\t')[8].split(';')[0].split('-')[2]
            if exon_number == "1" and line.split('\t')[18] == "right end":
                continue
            elif exon_number == coding_genes[gene_id] and line.split('\t')[18] == "left end":
                continue
            else:
                filtered_exons.write(line)
                
