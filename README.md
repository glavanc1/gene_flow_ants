# Gene flow in ants

This is the repository for the scripts used in this study:
>Preprint ref one day maybe who knows?

## Demultiplexed read acquisition
The SRA accession numbers for the demultiplexed reads used in this study will be available in `somewhere`.
They were demultiplexed using the `process_radtags` command of [stacks](https://catchenlab.life.illinois.edu/stacks/) v2.3e with the `-c -q -r -t 143 --filter_illumina`. options.

## Reference genomes
We mapped the reads onto published assemblies (one reference per genus). Accession numbers are given below (*italics = needs to be completed once they are published*):

| Genus        | Species used as reference | Bioproject accession number |
|:-------------|:--------------------------| ---------------------------:|
| Camponotus   | floridanus                | PRJNA445978                 |
| Formica      | selysi                    | PRJNA557079                 |
| Lasius       | niger                     | *GAGA*                      |
| Leptothorax  | acervorum                 | *GAGA*                      |
| Myrmica      | rubra                     | *GAGA*                      |
| Tapinoma     | erraticum                 | *got it from Jonathan*      |
| Temnothorax  | unifasciatus              | *GAGA*                      |
| Tetramorium  | immigrans                 | PRJNA533534                 |

## General workflow
All genera were analysed the same way. One directory was created per genus (hereafter `.`). Several subdirectories were created:
```bash
                      .
                      |
 --------------------------------------------
 |      |       |        |       |          |
raw   clean   stacks    out   results   admixture
```

## Mapping to the reference genomes
We used the `mem` algorithm of `bwa` version 0.7.17 and converted to bam format using `samtools` version 1.4. The script `bwa_mem_SLURM_script_maker.sh` was used to produce one SLURM script for each individual.

## SNP calling
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

## Filtering
We tested the effect of different filtering criteria on the number of retained SNPs using the script `filtering_tests.sh` We then visualised the results in `R` using the script `Filtering.R`.



