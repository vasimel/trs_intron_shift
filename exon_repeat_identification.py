import intervaltree
from intervaltree import Interval, IntervalTree
import Gff
import intervalListToIntervalTree
import re

#finding protein coding genes IDs and the number of their exons (key - ID, 
#value - number of exons)
def protein_coding_genes(gff_file_genes_and_exons):
    coding_genes = {}
    with open(gff_file_genes_and_exons) as gff:
        for line in gff:
            if "protein_coding" in line.split("\t")[8] and line.split("\t")[2] == "gene":
                match = re.search("GeneID:\d+[;, ]+", line.split("\t")[8])
                if match:
                    gene_id = match.group(0)
                    coding_genes[gene_id] = 0
                    
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
    return coding_genes


#create list of multiexonic genes
def multiexonic_list(coding_genes):
    return list(filter(lambda id: coding_genes[id] > 1, coding_genes.keys())) #len = 18230


#gff with multiexonic genes
def create_multiexonic_gff(gff_file_genes_and_exons, gff_file_multiexonic, multiexonic_list):
    with open(gff_file_genes_and_exons, "r") as genes_exons:
        with open(gff_file_multiexonic, "w") as multiexonic:
            lines = genes_exons.readlines()
            for line in lines:
                match = re.search("GeneID:\d+[;, ]+", line)
                if match.group(0) in multiexonic_list:
                    multiexonic.write(f"{line}")
    return 0
      
         
#creating a list of genomic features intervals
def features_interval_list(gff_file_multiexonic):
    interval_list=dict()
    i = 0
    with open(gff_file_multiexonic) as gff:
        for line in gff:
            interval_list[i] = Gff(line.strip())
            i += 1
    # remove exon entries with 0 interval length
    for id, elem in interval_list.items():
    if elem.start == elem.end:
        id_to_remove = id   
    interval_list.pop(id_to_remove)
    return interval_list


        
# searching for intersections
def find_intersections(exon_tree, output_file):
    with open(output_file, "w") as output:
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
    return 0
 

#filtering plus strand exons
def filter_plus_strand_exons(unfiltered_exons, filtered_exons, coding_genes):
    with open(unfiltered_exons, "r") as exons:
        with open(filtered_exons, "w") as filtered_exons:
            for line in exons:
                gene_id = re.search("GeneID:\d+[;, ]+", line.split("\t")[8]).group(0)
                exon_number = line.split('\t')[8].split(';')[0].split('-')[2]
                if exon_number == "1" and line.split('\t')[18] == "left end":
                    continue
                elif exon_number == coding_genes[gene_id] and line.split('\t')[18] == "right end":
                    continue
                else:
                    filtered_exons.write(line)
    return 0
              
            
            
            
#filtering minus strand exons
def filter_minus_strand_exons(unfiltered_exons, filtered_exons, coding_genes):
    with open(unfiltered_exons, "r") as exons:
        with open(filtered_exons, "w") as filtered_exons:
            for line in exons:
                gene_id = re.search("GeneID:\d+[;, ]+", line.split("\t")[8]).group(0)
                exon_number = line.split('\t')[8].split(';')[0].split('-')[2]
                if exon_number == "1" and line.split('\t')[18] == "right end":
                    continue
                elif exon_number == coding_genes[gene_id] and line.split('\t')[18] == "left end":
                    continue
                else:
                    filtered_exons.write(line)
    return 0
         
    
    
if __name__ == "__main__":
    wdir = "/home/ubuntu/repeats_in_proteins/"
    coding_genes = protein_coding_genes(f"{wdir}genes_exons.gff")
    create_multiexonic_gff(f"{wdir}genes_exons.gff", f"{wdir}mult_coding_genes_exons.gff", multiexonic_list(coding_genes))
    exon_tree = intervalListToIntervalTree(features_interval_list(f"{wdir}mult_coding_genes_exons.gff"))
    
    #tandem repeats dictionary
    trf_list=dict()
    i = 0
    with open(f"{wdir}trf_output.gff") as trf:
        for line in trf:
            trf_list[i] = Gff(line.strip())
            i += 1
    
    find_intersection(exon_tree, f"{wdir}intervaltree.out")
    #filter out plus strand exons with terminal repeats
    !awk '$3 == "exon" && $7 == "+" {print $0}' /home/ubuntu/repeats_in_proteins/intervaltree.out | grep -v "repeat is inside" > /home/ubuntu/repeats_in_proteins/plus_strand_exons_with_terminal_repeats.tsv
    #filter out minus strand exons with terminal repeats
    !awk '$3 == "exon" && $7 == "-" {print $0}' /home/ubuntu/repeats_in_proteins/intervaltree.out | grep -v "repeat is inside" > /home/ubuntu/repeats_in_proteins/minus_strand_exons_with_terminal_repeats.tsv

    
    filter_plus_strand_exons(f"{wdir}plus_strand_exons_with_terminal_repeats.tsv", 
                             f"{wdir}plus_strand_exons_with_terminal_repeats_filtered.tsv", coding_genes)
    
    filter_minus_strand_exons(f"{wdir}minus_strand_exons_with_terminal_repeats.tsv", 
                             f"{wdir}minus_strand_exons_with_terminal_repeats_filtered.tsv", coding_genes)
    
