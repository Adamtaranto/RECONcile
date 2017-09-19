# RECONcile

# Table of contents

* [About RECONcile](#about-reconcile)
* [Options and usage](#options-and-usage)
    * [Installing RECONcile](#installing-reconcile)
    * [Example Workflow](#example-workflow)
    * [Standard options](#standard-options)
* [License](#license)

# About RECONcile

Take clustered element fragment coordinates from RECON and write sequences to FASTA for alignment and consensus calling.

# Options and usage

## Installing RECONcile

Install from this repository.

```bash
git clone https://github.com/Adamtaranto/RECONcile.git && cd RECONcile
```

## Example Workflow

**Get total count and names of candidate TEs for RECON**  

```bash
awk '/^>/ {print $1;} ' Unclustered_candidate_TEs.fa | wc -l | sed 's/^ *//g' > denovo_names.txt  
awk '/^>/ {print $1;} ' Unclustered_candidate_TEs.fa | sed 's/>//g' | sed 's/^ *//g' >> temp_names.txt  
```

**Sort in lexical order**

```bash
cat temp_names.txt | sort >> denovo_names.txt  
rm temp_names.txt  
```

**Run All vs All BLAST**  

```bash
blastn -subject Unclustered_candidate_TEs.fa -query Unclustered_candidate_TEs.fa \  
-out denovo_all_v_all.tab -evalue 0.001 -outfmt 6 -task blastn -max_target_seqs 1000  
```

**Keep all non-self hits with bitscore > 100**  

```bash
awk '{if ($1 != $2 && $12 >= 100) print;}' denovo_all_v_all.tab > filtered_hits.tab  
```

**MSPCollect fails, use [blast2MSP.pl](https://gist.github.com/sestaton/ea770a26032983e49189#file-blast2msp-pl) to convert to blast result to MSP format instead:**  

```bash
./blast2MSP.pl filtered_hits.tab > MSP_denovo.out  
```

**Run RECON**  

```bash
./RECON-1.08/scripts/recon.pl denovo_names.txt MSP_denovo.out 1  
```

**Clean up eles output file from RECON**  

```bash
awk '!/^#/ {print;}' eles | sed 's/\( \)*/\1/g' | cut -d' ' -f 2- > Clean_Eles.txt  
```

**Fetch fragments for each family**  

```bash
./RECON_cluster.py -i Unclustered_candidate_TEs.fa -e Clean_Eles.txt -d output_clusters  
```  

**Align clustered fragments with MAFFT**  

```bash
mkdir mafftaligns  
for i in $(find output_clusters -name '*.fa' | sort);do  
base=$(basename $i)  
name="${base%.*}"  
mafft $i > 'mafftaligns/'$name'.aln'  
done  
```

**Build concensus sequences from alignments with [cons](http://www.bioinformatics.nl/cgi-bin/emboss/help/cons)**  

```bash
mkdir consensus_fasta  
for i in $(find mafftaligns -name '*.aln' | sort);do  
base=$(basename $i)  
name="${base%.*}"  
cons -sequence $i -identity 0 -snucleotide1 -supper1 -name $name'_consensus' \  
-outseq 'consensus_fasta2/'$name'_consensus.fa' -auto  
done  
```

**Concatenate all consensus sequences and convert to uppercase**  

```bash
rm RECON_consensus.fa  
touch RECON_consensus.fa  
for x in $(find consensus_fasta -name '*.fa' | sort);do  
tr [:lower:] [:upper:] < $x >> RECON_consensus.fa  
done  
``` 

## Standard options

```bash
Usage: RECONcile [-h] -i INFASTA [-e ELES] [-d OUTDIR]

Extracts and orients TE frags belonging to RECON clusters.

Optional arguments:
  -h, --help        Show this help message and exit
  -i, --inFasta     Multi fasta containing all TE sequences.
  -e, --eles        Space delimited 'ele' file from RECON.
  -d, --outDir      Directory for new cluster files to be written to.
```

# License

Software provided under MIT license.