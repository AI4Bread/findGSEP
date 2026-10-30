# findGSEP Web Server (GenomeGMM backend)

Shiny application for the GenomeGMM backend: k-mer based genome size
estimation for polyploid species.

## Files

| File | Purpose |
|---|---|
| `ui.R` | Shiny UI |
| `server.R` | Shiny server logic |
| `genome_gmm_script.R` | Command-line wrapper called by `server.R` |
| `genome_gmm.R` | GenomeGMM algorithm (single-file backend) |

All four files must stay in the same directory.

## Running

```r
shiny::runApp(".", port = 3838, host = "0.0.0.0", launch.browser = FALSE)
```

`server.R` resolves `Rscript` in this order: the `GSEGMM_RSCRIPT`
environment variable, `/home/shiny/miniconda3/envs/GSEGMM/bin/Rscript`,
then `Rscript` on `PATH`.

## Ploidy input

- **Estimate ploidy n** — ploidy is inferred from the k-mer spectrum.
- **Provide ploidy n** — supply the ploidy number plus a ploidy state:
  - `Autoploidy`: the fitted curve uses `n` equal copies.
  - `Alloploidy`: the subgenome number is additionally resolved from `n`.

The findGSEP paper itself does not use the autopolyploid/allopolyploid
distinction; it models "1..p homologous chromosomes". The ploidy state is an
input prior, not an inference from the k-mer histogram.

Example data are pre-filled as follows, based on published literature:

| Example | Ploidy | State | Reference |
|---|---|---|---|
| wheat | 6 | allopolyploid | AABBDD allohexaploid |
| potato | 4 | autopolyploid | autotetraploid |
| strawberry | 8 | allopolyploid | *Fragaria x ananassa* |
| javanica | 4 | allopolyploid | Nat Commun 2023, doi:10.1038/s41467-023-42700-w |
| Chinese sturgeon | 8 | autopolyploid | Genomics Proteomics Bioinformatics 2024, doi:10.1093/gpbjnl/qzad002 |

## Known issue

In `Alloploidy` mode, when the resolved base pattern is smaller than the
supplied ploidy, the reported genome size is inflated by
`ploidy / base_pattern`.

Cause: `genome_size = subgenome_unit * actual_ploidy`, while
`subgenome_unit` is derived from the base-pattern harmonic. The multiplier
should be `base_pattern`, not `actual_ploidy`.

Wheat (n = 6, base pattern resolved as 3, multiplier 2):

| Mode | Reported | Bread wheat 1C |
|---|---|---|
| Autoploidy n=6 | 15,237 Mb | about 15,400 Mb |
| Alloploidy n=6 | 29,372 Mb | about 15,400 Mb |

Not triggered when `base_pattern == ploidy` (e.g. javanica: 305.97 Mb,
matching flow cytometry 297 +/- 27 Mb).

## Known limitation

`www/` static assets are not included here, so the page header logo and the
custom CSS do not render.
