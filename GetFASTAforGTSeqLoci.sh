#---
# Title: Get fasta for GTSeq loci
# Author: Peiwen Li, Queen's University
# Date: Apr 7, 2021
# Description: This script extracts the flanking regions of the target SNP from a reference 
# genome using a list of SNP positions, and compiles the sequences into a fasta format for 
# GTSeq primer design
#---

##########################################################################################
#                                    Step 1. Prepare inputs                              #
##########################################################################################
# Input1: a whitelist of the targeted SNPs. Two columns separated by a tab, in which the 
# first column is the CHROM, and the second column is the POS. e.g.,
NC_036867.1	22906519
# The file is /Users/Peiwen/Desktop/Dropbox/TSFN/Fish Data_DNA/GTSeq/GTseq_SNPwhitelist.txt
# Now add the REF and ALT info to this file
nano 3074.txt # This file has all 3074 SNPs information
join -j 1 -o 1.2,1.3,1.4,1.5 <(<3074.txt awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<whitelist314.txt awk '{print $1"-"$2" "$0}' | sort -k1,1) > fullwhitelist314.txt
# Join can only join by comparing one column, so I need to first combine the chrom and pos columns into one in both files


# Input2: the genome fasta file which can be found in graham /home/peiwen/projects/rrg-lough/peiwen/RNA_seq/GCF_002910315.2_ASM291031v2_genomic.fna


##########################################################################################
#                              Step 2. Extract sequences                                 #
##########################################################################################
Nathan: As for the input format for primer design.  My script takes sequences in fasta format 
as input data.  The target SNP should be identified in square brackets and any non-target variation 
in the sequences should be masked with either an "N" or with a base ambiguity code (K,R,S,W....).

>ExampleLocus123
AGACCGAGCTTACCGATTTTCCGAACCAGATCNNACGATCCAGAT[C/T]CCTAGACCCAGATGAGGGGATCTACCAGAACCC
 
Ideally, we have about 100 bases on either side of the target SNP to design primers but I understand that is not possible with RAD data. 

cd ~/projects/rrg-lough/peiwen/
mkdir GTSeq
cd GTSeq/

# Identify non-target variations nearby
# First find the non-target variations, using the mac3 outputted VCF file as the reference (the VCF file before the last thinning step)
cp /home/peiwen/projects/def-lough/peiwen/TSFN6789/bwa_aln_samtools/filter_set_3/second_round_filtering/mac3.recode.vcf .
# Exclude the targeted variations (314 SNPs) from the VCF file
module load vcftools
vcftools --vcf mac3.recode.vcf --exclude-positions whitelist314.txt --recode --out rm314 #4302 variations left
rm mac3.recode.vcf

# Remove the headers, and keep only column 1,2,4,5 then sort
grep -v -P "#" rm314.recode.vcf | awk '{print $1,$2,$4,$5}' | sort > rm314.txt
# Sort mac3.txt by the first column, then by the second column numerically. 
# Find those variations out, using R
module load r
R
x<-read.csv("rm314.txt", header=F, sep=' ')
x<-x[order(x$V1,as.numeric(as.character(x$V2))),] # sort in R in numeric order (1,2,3,10 not 1,10,2,3)
x$V3<-nchar(x$V3) # replace the 3rd column by the # of characters of that colum
x<-x[,-4] # and remove the last ref column
write.table(x, file = "rm314-sorted.txt", sep = " ", row.names = F, col.names = F, quote = F)
quit()

# Technically, we just need to mark any variations that are <250 bps from our target variation as N. But figuring
# out which variations is <250bps from the target is hard. Instead, I will mark all these 4302 variation as N in the reference sequences

# Now extract the full scaffolds for all 314 targeted + 4302 nontargeted variations
awk '{print $1}' fullwhitelist314.txt | sort | uniq > chrom.txt
awk '{print $1}' rm314-sorted.txt | sort | uniq >> chrom.txt
awk '{print $1}' chrom.txt | sort | uniq > uniquechrom.txt
wc -l uniquechrom.txt #673

while read -r chrom; do
    fullseq="$(sed -n -e "/$chrom/,/>/ p" ../RNA_seq/GCF_002910315.2_ASM291031v2_genomic.fna | sed 1d | sed '$d'| tr -d '\n')" # Get the full scaffold sequence, by first extract the lines between $chrom and the next line starting with a >, then remove the first and last lines, and remove the newline characters. TESTED: if $chrom=NC_00086.1, which is the last scaffold in the GCF..fna file, this command STILL work despite there is no line starting with a '>' following that scaffold, just dont need to remove the last line (no sed '$d' is needed)
    echo \>$chrom >> fullseq.txt
    echo $fullseq >> fullseq.txt
done < uniquechrom.txt
wc -l fullseq.txt #should be 673*2=1346 lines

####
#!/bin/bash
#SBATCH -J mask 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --cpus-per-task 1 
#SBATCH --mem-per-cpu=4g 
#SBATCH --account=rrg-lough
#SBATCH -t 2-0:0:0 
#SBATCH --mail-type=BEGIN,END,FAIL 
#SBATCH --mail-user=15pl16@queensu.ca 
#SBATCH -o mask_%j.o
# Mask the non-target variations. Store the new $fullseq to fullseq.txt
while read -r chrom pos num; do
    fullseq="$(grep -A 1 "$chrom" fullseq.txt | sed 1d)"
    # here mask the non-target variations in the fullseq
    for i in $(seq 1 $num); do 
    if [ $num -gt 1 ] # if the reference allele has more than 1 bases, e.g. REF=ATT and ALT=ATTTT, need to replace all ATT to NNN
    then
    fullseq="$(echo $fullseq | sed "s/./N/$((pos+i-1))")"
    fi
    if [ $num -eq 1 ] # if the reference allele has 1 bases, e.g. REF=A and ALT=T, need to replace just the A to N
    then
    fullseq="$(echo $fullseq | sed "s/./N/$pos")"
    fi
    done
    sed -i "/$chrom/,+1 d" fullseq.txt # deleted the original sequence and header
    echo \>$chrom >> fullseq.txt # add the header back
    echo $fullseq >> fullseq.txt # add the masked sequence back, so that in the next iteration, it works on the new masked sequence from last iteration (this is designed for nontarget variations taht are on the same scaffold)
done < rm314-sorted.txt
#### put line 80-108 in a job

#### Check:
wc -l fullseq.txt # should still be 1346 lines
# see if all positions in one same chromosome were changed to N
# e.g.1:
NC_036851.1 13008257
chrom=NC_036851.1
fullseq="$(sed -n -e "/$chrom/,/>/ p" ./fullseq.txt | sed 1d | sed '$d'| tr -d '\n')" 
fullseqorigin="$(sed -n -e "/$chrom/,/>/ p" ..//RNA_seq/GCF_002910315.2_ASM291031v2_genomic.fna | sed 1d | sed '$d'| tr -d '\n')" 
echo "${fullseq:(13008257-101):250}" # the 101th position should be an N
grep 'NC_036851.1' rm314-sorted.txt | awk '{s+=$3}END{print s}' #89
# and since there are 89 non-target variations in this chromosome,
echo $fullseq | grep -o 'N' | wc -l
echo $fullseqorigin | grep -o 'N' | wc -l # there should be 89 more N in $fullseq than in $fullseqorigin
# e.g.2
NC_036853.1 37109232 
chrom=NC_036853.1
fullseq="$(sed -n -e "/$chrom/,/>/ p" ./fullseq.txt | sed 1d | sed '$d'| tr -d '\n')" 
fullseqorigin="$(sed -n -e "/$chrom/,/>/ p" ..//RNA_seq/GCF_002910315.2_ASM291031v2_genomic.fna | sed 1d | sed '$d'| tr -d '\n')" 
echo "${fullseq:(37109232-101):250}" # the 101th position should be an N
grep 'NC_036853.1' rm314-sorted.txt | awk '{s+=$3}END{print s}' #174
# and since there are 174 non-target variations in this chromosome,
echo $fullseq | grep -o 'N' | wc -l
echo $fullseqorigin | grep -o 'N' | wc -l # there should be 174 more Ns in $fullseq than in $fullseqorigin, but only 173 more. why? maybe some positions that needs to be replaced is already N: in pos 43706811, 7 bases were replaced by N, including 43706815. and 43706815 appears again--> counted twice
####

echo '>' >> fullseq.txt # add a '>' to the end of fullseq.txt so that line 144 wont get an error
# Now extract the 250 bps sequences in which the 101th position is the targeted SNP. Store the sequences in snps.fa. This then can be sent to Nathan
while read -r chrom pos ref alt; do
    echo \>$chrom\_$pos\_$ref\_$alt >> snps.fa
    fullseq="$(sed -n -e "/$chrom/,/>/ p" ./fullseq.txt | sed 1d | sed '$d'| tr -d '\n')" # Get the full masked scaffold sequence, by first extract the lines between $chrom and the next line starting with a >, then remove the first and last lines, and remove the newline characters. TESTED: if $chrom=NW_019957564.1, which is the last scaffold in the fullseq.txt file, this command STILL work because i added an extra '>' to the end of fullseq.txt in line 155
    newseq=$(echo "${fullseq:(pos-101):250}") # Extract: start at position $pos-101 (indexed at 0 here, but in a VCF file the pos is indexed at 1) and print 250 bases, the 101th character is the targeted SNP. Some SNPs are indels, so I get 250 bases instead of 201
    newseq1=$(echo $newseq | tr [:lower:] [:upper:]) # Change all  lower cases to upper cases
    newseq2="${newseq1:0:100}"["${ref}"/"${alt}"]"${newseq1:100}" 
    echo $newseq2 | sed "s/\]$ref/\]/g" >> snps.fa # line 161-162: highlight the target SNP, and store the sequence to snps.fa
done < fullwhitelist314.txt
wc -l snps.fa #should be 314*2=628 lines


####
# Now we have 314+51=365 positions. I want to add another 132 snps from the 3074 snps dataset, to make it 500 (=50 from genes + 446 from DNA-seq + 4 spaces for SDY gene)
# The additional 132 loci are in file add132.txt
join -j 1 -o 1.2,1.3,1.4,1.5 <(<3074.txt awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<add132.txt awk '{print $1"-"$2" "$0}' | sort -k1,1) > full132.txt
# Now extract these positions from the masked fullseq.txt
# Since these positions were already masked as N, we need to change the code a bit
while read -r chrom pos ref alt; do
    echo \>$chrom\_$pos\_$ref\_$alt >> snps.fa
    fullseq="$(sed -n -e "/$chrom/,/>/ p" ./fullseq.txt | sed 1d | sed '$d'| tr -d '\n')" # Get the full masked scaffold sequence, by first extract the lines between $chrom and the next line starting with a >, then remove the first and last lines, and remove the newline characters. TESTED: if $chrom=NW_019957564.1, which is the last scaffold in the fullseq.txt file, this command STILL work because i added an extra '>' to the end of fullseq.txt in line 155
    newseq=$(echo "${fullseq:(pos-101):250}") # Extract: start at position $pos-101 (indexed at 0 here, but in a VCF file the pos is indexed at 1) and print 250 bases, the 101th character is the targeted SNP. Some SNPs are indels, so I get 250 bases instead of 201
    newseq1=$(echo $newseq | tr [:lower:] [:upper:]) # Change all  lower cases to upper cases
    newseq2="${newseq1:0:100}"["${ref}"/"${alt}"]"${newseq1:100}" 
    echo $newseq2 | sed "s/\]N/\]/g" >> snps.fa # line 161-162: highlight the target SNP, and store the sequence to snps.fa
done < full132.txt
wc -l snps.fa # should be (314+132)*2=892 lines

