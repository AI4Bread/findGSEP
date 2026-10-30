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

## Example data

Pre-filled ploidy and ploidy state, based on published literature:

| Example | Ploidy | State | Reference |
|---|---|---|---|
| wheat | 6 | allopolyploid | AABBDD allohexaploid |
| potato | 4 | autopolyploid | autotetraploid |
| strawberry | 8 | allopolyploid | *Fragaria x ananassa* |
| javanica | 4 | allopolyploid | Nat Commun 2023, doi:10.1038/s41467-023-42700-w |
| Chinese sturgeon | 8 | autopolyploid | Genomics Proteomics Bioinformatics 2024, doi:10.1093/gpbjnl/qzad002 |

## Note

`www/` static assets are not included in this directory, so the page header
logo and the custom CSS are not rendered.
