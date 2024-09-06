Genome pre-processing pipeline
==============================

This is a local pipeline to pre-process downloaded genomes before feeding them
to <https://github.com/nf-core/pairgenomealign>.

## Introduction

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

## Usage

First make a sample sheet with usual _nf-core_ pipelines.

`samplesheet.csv`:

```csv
sample,fasta
query_1,path-to-query-genome-file-one.fasta
query_2,path-to-query-genome-file-two.fasta
```

Then run the pipeline as usual:

```bash
nextflow run oist/Luscombe_pairgenomealign-preprocess \
   -profile oist \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

Test the pipeline (adapt the `-w` option for your own case!

```bash
nextflow run oist/Luscombe_pairgenomealign-preprocess \
   -profile oist,test \
   -w /flash/LuscombeU/`whoami`/cache/deletemeTest \
   --outdir results_test
```

## Pipeline output

TBD

## Credits

nf-core/pairgenomealignmask was originally written by [Mahdi](https://github.com/U13bs1125).
