CEL-Seq2 counts (both raw and filtered) and R codes for differential gene expression analysis

CEL-Seq2\_counting\_report.txt

Contains transcript raw counts for each Aqu2.1 gene for each larval cell type sample

E - epithelial cells

F - flask cells

G - globular cells

I - cells of the inner cell mass

P - cells constituting the posterior pigment ring

S - spicule cells

CEL-Seq2\_counting\_report\_highest75%var.txt

Raw counts with least variable genes (25%) removed to improve DESeq2 sensitivity (as per Sha et al. 2016: 15% to 30% recommended)

vst-preDESeq2filter.R

Use:

Transform raw counts with DESeq2 package function (variance-stabilising transformation). Using transformed counts, remove 25% genes with least variable counts.

Input:

Raw count table provided [CEL-Seq2\_counting\_report.txt]

Output:

vst transformed counts (and of top 75% most variable genes)

DESeq2.R

Use:

Perform groupwise comparison (e.g. e cells vs not-e cells) to contrast cell types. Compute genes differentially expressed for selected cell type.

Input:

Count table provided [CEL-Seq2\_counting\_report\_highest75%var.txt]

Output:

Differentially expressed genes for selected cell types (with fold change and p-value).