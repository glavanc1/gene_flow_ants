# Widespread but cryptic introgression shapes genetic diversity in natural populations

This is the repository for the scripts used in the study by Guillaume Lavanchy, Ludovic Ruedi, Olivier Broennimann, Kristine Jecha, Marianna Tzivanopoulou, Jérôme Goudet, Tanja Schwander
`the link will be added when available`

Each genus was analysed separately, but all genera were analysed the same way.

## 1.Demultiplexed read acquisition
The SRA accession numbers for the demultiplexed reads are available in the table below.
They were demultiplexed using the `process_radtags` command of [stacks](https://catchenlab.life.illinois.edu/stacks/) v2.3e with the `-c -q -r -t 143 --filter_illumina`. options.

| Genus          |NCBI accession number      |
|:---------------|:--------------------------| 
| *Camponotus*   | PRJNA1458683              |
| *Formica*      | PRJNA1466771              |
| *Lasius*       | PRJNA1176511              |
| *Myrmica*      | PRJNA1188312              |
| *Tapinoma*     | PRJNA1476096              |
| *Temnothorax*  | PRJNA1476131              | 
| *Tetramorium*  | PRJNA1476143              |


## 2. Competitive mapping
To filter out potential cross-contaminations during wetlab, we performed competitive mapping, mapping the reads simultaneously to the concatenated assemblies of all genera and retaining only reads that map to their expected assembly. Accession numbers are given below:

| Genus          | Species used as reference | NCBI accession number or link |
|:---------------|:--------------------------| ---------------------------:|
| *Camponotus*   | *fallax*                  | GCF_003227725               |
| *Formica*      | *selysi*                  | GCA_009859135               |
| *Lasius*       | *niger*                   | GCA_041902855.1             |
| *Myrmica*      | *rubra*                   | GCA_048181765.1             |
| *Tapinoma*     | *erraticum*               | https://zenodo.org/records/5705739      |
| *Temnothorax*  | *unifasciatus*            | GCA_048541725.1             |
| *Tetramorium*  | *immigrans*               | GCA_011636585.1             |

We concatenated all assemblies
```
cat *.fasta > cat_gen.fasta
bwa index cat_gen.fasta
``` 

and mapped each sample to this concatenated file

` bwa mem -t 20 cat_gen.fasta $READ.fq.gz| samtools sort --output-fmt=SAM --threads 20 | samtools view --no-header -o ./sam/catgen_$READ.sam - `  

We then identified the contaminated reads using the script `FilterContaminates.R` and filtered out the reads that did not map to the 

We then filtered out reads that did not map to the genus of interest using the `filterbyname.sh` script from the BBmap suite (https://github.com/BioInfoTools/BBMap/blob/master/sh/filterbyname.sh)
```
module load bbmap 
filterbyname.sh in=${SAMPS}/${READ}.fq.gz names=${IDDIR}/${READ}_ID.txt out=${OUT}/${READ}.filtered.fq.gz include=t
```

## 3.Mapping to the reference genomes
We then mapped the retained reads to the assembly of their own genome only (same as above).

We used the script `map_and_clean.sh` to map the reads with the `mem` algorithm of `bwa` version 0.7.17 and to process them using using `samtools` version 1.4. 

We checked the proportion of reads mapping onto the reference genome using the script `mapping_stats.sh` and visualized the results in `R` using `Mapping.R`.

## 4.SNP calling
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

## 5.Removing contaminations and filtering
We visualized the effect of different filtering criteria on the number of retained SNPs using the script `filtering_tests.sh` and visualised the results in `R` using the script `Filtering.R`.
We selected *a posteriori* a minimum depth cutoff of 8 reads (which gives below 1% of allelic dropout), a minor allele count of 3, and a maximum 75% of individuals with missing genotypes per locus. This was performed in `vcftools` v0.1.14:
```
vcftools --vcf FILENAME.vcf --minDP 8 --max-meanDP 200 --mac 2 --max-missing 0.75 --recode --out OUTPUTNAME.vcf
```

After this first filtering, we discarded individuals with more than 75% of missing data, estimated the presence of contaminations and filtered them based on Allelic Depth Ratio (ADR) using the script `filter_vcf_ADR.R` with the default cutoff (0.01). See Jecha et al. 2024 (https://www.biorxiv.org/content/10.1101/2024.11.27.625433v1.full) for the rationale behind the approach.

We then ran the same `vcftools` command as above to filter again after having discarded the individuals with high missing data and corrected depth at some loci. The output was considered the final dataset for the first step (Taxonomy).

## 6.Taxonomy
We ran a series of analyses:
### phylogeny of a fragment of COI.
We aligned our COI sequences, along with reference individuals from published studies when available, with `mafft` v7.475:
`mafft --localpair  --maxiterate 16 --phylipout --inputorder all_COI.fasta > all_clean_trimmed_alignment_COI.phylip`.
We then reconstructed a ML phylogeny with `iqtree` v2.0.6: `iqtree2 -T 1 -s clean_trimmed_alignment_COI.phylip -m MFP -B 1000`.

### Clustering
We ran `admixture` with a range of K values (depending on the number of expected species in each genus, based on the number of morphological species). This was done using the script `runadmixture.sh`.

### Multi-Dimensional Scaling.
This, along with the visualization of the results of the other analyses and the selection of individuals for further steps, was done in `R` with the scripts `GENUSNAME_Analyses.R`.

## 7.Introgression rates
We ran `admixture` again, this time in supervised mode using the population maps generated in the scripts `GENUSNAME_Analyses.R`, using the script 
## 8.Contribution to genetic diversity

## 9.Correlates of hybridization and visualization


