# RECONcile
Take clustered element fragment coordinates from RECON and write sequences to fasta for alignment and consensus calling.

# Workflow  

## Get total count and names of candidate TEs for RECON  
awk '/^>/ {print $1;} ' Unclustered_candidate_TEs.fa | wc -l | sed 's/^ *//g' > denovo_names.txt  
awk '/^>/ {print $1;} ' Unclustered_candidate_TEs.fa | sed 's/>//g' | sed 's/^ *//g' >> temp_names.txt  

## Sort in lexical order  
cat temp_names.txt | sort >> denovo_names.txt  
rm temp_names.txt  

## Run All vs All BLAST  
blastn -subject Unclustered_candidate_TEs.fa -query Unclustered_candidate_TEs.fa \  
-out denovo_all_v_all.tab -evalue 0.001 -outfmt 6 -task blastn -max_target_seqs 1000  

## Keep all non-self hits with bitscore > 100  
awk '{if ($1 != $2 && $12 >= 100) print;}' denovo_all_v_all.tab > filtered_hits.tab  

## MSPCollect fails, use [blast2MSP.pl](https://gist.github.com/sestaton/ea770a26032983e49189#file-blast2msp-pl) to convert to blast result to MSP format instead:  
./blast2MSP.pl filtered_hits.tab > MSP_denovo.out  

## Run RECON  
./RECON-1.08/scripts/recon.pl denovo_names.txt MSP_denovo.out 1  

## Clean up eles output file from RECON  
awk '!/^#/ {print;}' eles | sed 's/\( \)*/\1/g' | cut -d' ' -f 2- > Clean_Eles.txt  

## Fetch fragments for each family  
./RECON_cluster.py -i Unclustered_candidate_TEs.fa -e Clean_Eles.txt -d output_clusters  
  
## Align clustered fragments with MAFFT  
mkdir mafftaligns  
for i in $(find output_clusters -name '*.fa' | sort);do  
base=$(basename $i)  
name="${base%.*}"  
mafft $i > 'mafftaligns/'$name'.aln'  
done  

## Build concensus sequences from alignments with [cons](http://www.bioinformatics.nl/cgi-bin/emboss/help/cons)  
mkdir consensus_fasta  
for i in $(find mafftaligns -name '*.aln' | sort);do  
base=$(basename $i)  
name="${base%.*}"  
cons -sequence $i -identity 0 -snucleotide1 -supper1 -name $name'_consensus' \  
-outseq 'consensus_fasta2/'$name'_consensus.fa' -auto  
done  

## Concatenate all consensus sequences and convert to uppercase  
rm RECON_consensus.fa  
touch RECON_consensus.fa  
for x in $(find consensus_fasta -name '*.fa' | sort);do  
tr [:lower:] [:upper:] < $x >> RECON_consensus.fa  
done  

**Note:** need script to trim terminal Ns from consensus sequences.  
**Note:** Add '#Unknown' to all consensus names to make RM compatible.  
