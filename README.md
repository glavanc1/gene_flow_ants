# Gene flow in ants

This is the repository for the scripts used in this study:
>Widespread but cryptic introgression shapes genetic diversity in natural populations
>by Guillaume Lavanchy, Ludovic Ruedi, Olivier Broennimann, Kristine Jecha, Marianna Tzivanopoulou, Jérôme Goudet, Tanja Schwander

## 1.Demultiplexed read acquisition
The SRA accession numbers for the demultiplexed reads used in this study will be available in `somewhere`.
They were demultiplexed using the `process_radtags` command of [stacks](https://catchenlab.life.illinois.edu/stacks/) v2.3e with the `-c -q -r -t 143 --filter_illumina`. options.


## General workflow
Each genus was analyzed separately but all genera were analysed the same way. One directory was created per genus (hereafter `.`).

## 2.Mapping to the reference genomes
We mapped the reads onto published assemblies (one reference per genus). Accession numbers are given below:

| Genus          | Species used as reference | NCBI accession number or link |
|:---------------|:--------------------------| ---------------------------:|
| *Camponotus*   | *fallax*                  | GCF_003227725               |
| *Formica*      | *selysi*                  | GCA_009859135               |
| *Lasius*       | *niger*                   | GCA_041902855.1             |
| *Myrmica*      | *rubra*                   | GCA_048181765.1             |
| *Tapinoma*     | *erraticum*               | https://zenodo.org/records/5705739      |
| *Temnothorax*  | *unifasciatus*            | GCA_048541725.1             |
| *Tetramorium*  | *immigrans*               | GCA_011636585.1             |


We used the script `map_and_clean.sh` to map the reads with the `mem` algorithm of `bwa` version 0.7.17 and to process them using using `samtools` version 1.4. 
We checked the proportion of reads mapping onto the reference genome using the script `mapping_stats.sh` and visualized the results in `R` using `Mapping.R`.


## 3.SNP calling
We created a `popmap.txt` file with all individuals belonging to the same population.

```bash
ls clean/ | grep '.fq.gz' | sed 's/.fq.gz//' > individuals.tmp
NIND=$(wc -l individuals.tmp | cut -f 1 -d ' ')
yes 1 | head -n $NIND > ones.tmp
paste individuals.tmp ones.tmp > popmap.txt
rm *.tmp
```
We then called SNPs using `stacks` v2.3e. We ran `gstacks` and `populations` using the following commands:

```bash
gstacks -t 16 --phasing-dont-prune-hets --ignore-pe-reads -M ./popmap.txt -I ./mapped -O ./stacks
populations -t 16 -M ./popmap.txt -P ./stacks -O ./out --write-single-snp --vcf
```

## 4.Removing contaminations and filtering
We visualized the effect of different filtering criteria on the number of retained SNPs using the script `filtering_tests.sh` and visualised the results in `R` using the script `Filtering.R`.
We selected *a posteriori* a minimum depth cutoff of 8 reads (which gives below 1% of allelic dropout), a minor allele count of 3, and a maximum 75% of individuals with missing genotypes per locus. This was performed in `vcftools` v0.1.14:
```
vcftools --vcf FILENAME.vcf --minDP 8 --max-meanDP 200 --mac 2 --max-missing 0.75 --recode --out OUTPUTNAME.vcf
```

After this first filtering, we discarded individuals with more than 75% of missing data. We then estimated the presence of contaminations and filtered them based on Allelic Depth Ratio (ADR) using the script `ADR_filtering.R` with the default cutoff (0.01). Both these things were done in `R` with the script `ADR_filtering.R`.

We then ran the same `vcftools` command as above to filter again after having discarded the individuals with high missing data and corrected depth at some loci. The output was considered the final dataset for the first step (Taxonomy).

## 5.Taxonomy
We ran a series of analyses:
- phylogeny of a fragment of COI. We aligned our COI sequences, along with reference individuals from published studies when available, with `mafft` v7.475: `mafft --localpair  --maxiterate 16 --phylipout --inputorder all_COI.fasta > all_clean_trimmed_alignment_COI.phylip`. We then reconstructed a ML phylogeny with `iqtree` v2.0.6: `iqtree2 -T 1 -s clean_trimmed_alignment_COI.phylip -m MFP -B 1000`.
- `admixture` with a range of K values (depending on the number of expected species in each genus, based on the number of morphological species). This was done using the script `runadmixture.sh`.
- Multi-Dimensional Scaling. This, along with the visualization of the results of the other analyses and the selection of individuals for further steps, was done in `R` with the scripts `GENUSNAME_Analyses.R`.
- 
## 6.Hybridization rates

## 7.Correlates of hybridization


