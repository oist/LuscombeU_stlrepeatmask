# Genome repeat masking pipeline

This is a local pipeline to repeat-mask genomes before feeding them to
<https://github.com/nf-core/pairgenomealign>.  The rationale is that it is
better when all soft masks have been produced by the same pipeline…

## What it does:

This pipeline takes genomes as inputs and soft-masks their repeats with the following software:

- [tantan](https://gitlab.com/mcfrith/tantan), version 51.  Tantan is our default choice from a long time because TRF used to be non-free.
- [WindowMasker](https://doi.org/10.1093/bioinformatics/bti774), version 1.0.0 distributed with [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK569845/#ckbk_Createmaskedb.Create_masking_inform_1) 2.17.0.
- [RepeatMasker](https://www.repeatmasker.org/) version 4.1.9.

The input of repeatmasker can be any of:

- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) version 2.0.7 (default)
- [Dfam](https://www.dfam.org/home) (optional)
- A custom repeat library (optional)

RepeatMasker and RepeatModeler are run from the same bioconda package as the standard _nf-core_ module.

The pipeline then merges the soft masks of the RepeatMasker runs, and then merges that with the tantan and WindowMasker runs.

Finally, the pipeline prepares a MultiQC report that shows the extent of masking for each tool.

## Disclaimer

This is not an official pipeline. This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

## Usage

**This pipeline does not work with conda.**

First make a sample sheet with usual _nf-core_ pipelines.

`samplesheet.csv`:

```csv
sample,fasta
query_1,path-to-query-genome-file-one.fasta.gz
query_2,path-to-query-genome-file-two.fasta.gz
```

If the input is not compressed then pass the `--gzipped_input=false` parameter.
Note that mixing compressed and uncompressed input is not supported, partly
because [WindowMasker does not handle `stdin` input](https://github.com/ncbi/ncbi-cxx-toolkit-public/issues/21).

Then run the pipeline as usual:

```bash
nextflow run oist/LuscombeU_stlrepeatmask \
   -profile oist \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

Test the dev branch of the pipeline (adapt the `-w` option for your own case!

```bash
nextflow run oist/LuscombeU_stlrepeatmask \
   -profile oist,test \
   -w /flash/LuscombeU/`whoami`/cache/deletemeTest \
   --outdir results_test
```

Test the local checkout

```bash
nextflow run ./main.nf \
   -profile oist,test \
   -w /flash/LuscombeU/`whoami`/cache/deletemeTest \
   --outdir results_test
```

## Options

- Point `--repeatlib` to a FASTA file to have an extra RepeatMasker run using it as a library.
- Set `--taxon` to a taxon name to have an extra RepeatMasker run using the `-species` option set to that taxon.
- Set the `--gzipped_input=false` parameter when the input is not compressed..
- Point `--dfam` to a directory containing a `famdb` subdirectory with the FamDB files in HDF5 format (not gzipped).

### Dfam

Most containers that provide RepeatMasker do not contain a full copy of Dfam,
which is huge.  However they sometimes have a stub, for instance under
`/usr/local/share/RepeatMasker/Libraries` or `/opt/RepeatMasker/Libraries`.
Interstingly, when using the `--libdir` option, Dfam has to be in a subfolder
named `famdb`, although it is named `FamDB` in the download website.

## Pipeline output

### `tantan`, `repeatmodeler`, `windowmasker`, `dfam` (optional), `extlib` (optional), `mergedmasks`

- Masked genome file (compressed with `bgzip`).
- BED file representing the masked regions.
- Summary statistics of the softmasked genome.

### Only in `repeatmodeler`

- De novo detected repeats (`.fa`, `.log`, `.stk` and BLAST database files `.n*`)

## Resource usage

On a test run on haplotype-merged and diploid assemblies of _Oikopleura dioica_ (2n = 60 Mbp):

- CPU usage was ~50 % for most processes. RepeatModeller was allocated 24 cores and used ~10 on average.
- Memory usage was less than 1 GB for all processes except RepeatModeller (~6 GB, max 8 GB).
- All processes needed only 10 % of the allocated time, except for RepeatModeller, which took between 100 and 500 minutes.
- On a couple of primate genomes, RepeatModeller managed to keep its 24 cores 60% busy for ~30 hours using 40 GB memory.

## Future directions

- It may be interesting to add TRF and ULTRA, and compare and combine their results to the ones of tantan.

## Credits

This pipeline was originally written by [Mahdi](https://github.com/U13bs1125) and then
taken over by @charles-plessy.
