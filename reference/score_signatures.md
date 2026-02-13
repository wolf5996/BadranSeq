# Score Gene Set Signatures

Scores cells against MSigDB gene sets using
[`AddModuleScore`](https://satijalab.org/seurat/reference/AddModuleScore.html).
Fetches gene sets via
[`msigdbr`](https://igordot.github.io/msigdbr/reference/msigdbr.html),
computes per-cell module scores, and stores them as metadata columns
with clean pathway names.

## Usage

``` r
score_signatures(
  object,
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL,
  prefix = NULL,
  report_only = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object. Gene names (rownames) must be HGNC symbols. Ensembl
  IDs and Entrez IDs are detected and rejected with an informative
  error.

- species:

  character. Species name for msigdbr (default: `"Homo sapiens"`).
  Available species:

  - `"Homo sapiens"` – human

  - `"Mus musculus"` – mouse

  - `"Rattus norvegicus"` – rat

  - `"Danio rerio"` – zebrafish

  - `"Drosophila melanogaster"` – fruit fly

  - `"Saccharomyces cerevisiae"` – yeast

  - `"Caenorhabditis elegans"`

  - `"Sus scrofa"` – pig

  - `"Bos taurus"` – cow

  - `"Canis lupus familiaris"` – dog

  - `"Felis catus"` – cat

  - `"Gallus gallus"` – chicken

  - `"Equus caballus"` – horse

  - `"Pan troglodytes"` – chimpanzee

  - `"Macaca mulatta"` – rhesus macaque

  - `"Anolis carolinensis"` – green anole

  - `"Monodelphis domestica"` – opossum

  - `"Ornithorhynchus anatinus"` – platypus

  - `"Xenopus tropicalis"` – tropical clawed frog

  - `"Schizosaccharomyces pombe 972h-"`

- collection:

  character. MSigDB collection (default: `"H"`). Available collections:

  - `"H"` – Hallmark gene sets (50 sets). Well-defined biological states
    and processes with reduced redundancy.

  - `"C1"` – Positional gene sets (302 sets). Genes organised by
    chromosomal position.

  - `"C2"` – Curated gene sets. From publications, pathway databases,
    and expert knowledge. Subcollections:

    - `"CGP"` – Chemical and Genetic Perturbations

    - `"CP"` – Canonical Pathways

    - `"CP:BIOCARTA"` – BioCarta Pathways

    - `"CP:KEGG_LEGACY"` – KEGG Legacy Pathways

    - `"CP:KEGG_MEDICUS"` – KEGG Medicus Pathways

    - `"CP:PID"` – PID Pathways

    - `"CP:REACTOME"` – Reactome Pathways

    - `"CP:WIKIPATHWAYS"` – WikiPathways

  - `"C3"` – Regulatory target gene sets. Genes sharing a regulatory
    motif. Subcollections:

    - `"MIR:MIRDB"` – miRDB microRNA targets

    - `"MIR:MIR_LEGACY"` – Legacy microRNA targets

    - `"TFT:GTRD"` – GTRD transcription factor targets

    - `"TFT:TFT_LEGACY"` – Legacy transcription factor targets

  - `"C4"` – Computational gene sets. Defined by mining large
    collections of cancer-oriented microarray data. Subcollections:

    - `"3CA"` – Curated Cancer Cell Atlas

    - `"CGN"` – Cancer Gene Neighborhoods

    - `"CM"` – Cancer Modules

  - `"C5"` – Ontology gene sets. Derived from controlled vocabularies.
    Subcollections:

    - `"GO:BP"` – GO Biological Process

    - `"GO:CC"` – GO Cellular Component

    - `"GO:MF"` – GO Molecular Function

    - `"HPO"` – Human Phenotype Ontology

  - `"C6"` – Oncogenic signature gene sets (189 sets). Defined from
    microarray data of cancer gene perturbations.

  - `"C7"` – Immunologic signature gene sets. From immunology studies.
    Subcollections:

    - `"IMMUNESIGDB"` – ImmuneSigDB

    - `"VAX"` – HIPC Vaccine Response

  - `"C8"` – Cell type signature gene sets (866 sets). Curated from
    single-cell sequencing studies.

- subcollection:

  character. MSigDB sub-collection (default: `NULL`). Only applicable
  for collections with sub-collections (C2, C3, C4, C5, C7). See
  `collection` for all valid subcollection values.

- prefix:

  character. Prefix for metadata column names (default: `NULL` =
  automatically derived from `collection`, e.g. `"Hallmark_"` for
  collection `"H"`).

- report_only:

  logical. If `TRUE`, return a tibble reporting which genes from each
  gene set are found in the object instead of scoring. The tibble has
  columns `gene_set`, `gene`, and `found`. All original genes from the
  msigdbr fetch are included (not just those present in the object).
  Default: `FALSE`.

- ...:

  Additional arguments passed to
  [`AddModuleScore`](https://satijalab.org/seurat/reference/AddModuleScore.html).

## Value

If `report_only = FALSE` (default), the Seurat object with new metadata
columns (one per gene set). Column names are cleaned pathway names with
the specified prefix. If `report_only = TRUE`, a tibble with columns
`gene_set` (character), `gene` (character), and `found` (logical).

## Examples

``` r
if (FALSE) { # \dontrun{
# Score Hallmark gene sets (default)
obj <- score_signatures(obj)

# Score KEGG pathways
obj <- score_signatures(obj, collection = "C2",
                        subcollection = "CP:KEGG_MEDICUS")

# Score mouse cell type signatures
obj <- score_signatures(obj, species = "Mus musculus", collection = "C8")

# Custom prefix
obj <- score_signatures(obj, collection = "C6", prefix = "Onco_")

# Report-only mode: see which genes are found
report <- score_signatures(obj, collection = "H", report_only = TRUE)
} # }
```
