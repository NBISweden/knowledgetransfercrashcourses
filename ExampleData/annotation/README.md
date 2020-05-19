# Annotation gtf

The mouse genome annotation gtf-file is downloaded from the same source as the Reference fasta file and pruned to contain only entries for chromosome 19, as follows:

```
wget -O - ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz | zcat |awk '{if($1=="19" || match($1, "#")==1) print $0}' | gzip > Mus_musculus.GRCm38.99.chromosome.19.gtf.gz
```
