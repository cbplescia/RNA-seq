# RNA-seq methods (draft for the paper)

Items in [square brackets] are not recorded in the scripts and need filling in or
checking. Everything else is taken from scripts 00, 02-07 and config.sh.

## RNA sequencing and read processing

[Cell type / treatments / doses / time points; RNA extraction; library preparation kit
and strandedness chemistry.] Libraries were sequenced on an Element Biosciences AVITI
system (2 x 75 bp paired-end) and demultiplexed with bases2fastq. Adapters were removed
and low-quality bases were trimmed with cutadapt [version], using an Illumina-compatible
adapter on both mates, a Phred quality threshold of 20 at both read ends (-q 20,20) and
a minimum read length of 50 nt. The dataset comprised 33 libraries: 11 conditions
(interferon: none, IFN-beta, IFN-gamma; crossed with virus or virus-like-particle
conditions: none, EBOV, MARV, EVLP, MVLP) with three biological replicates each.
[Confirm condition layout.]

## Reference construction

Transcript sequences were extracted from the Ensembl GRCm39 release 109 annotation
(GTF) and primary-assembly genome using gffread. A decoy-aware Salmon index was built
(Salmon v1.10.3; k = 31) from a combined "gentrome" of (i) the Ensembl transcriptome,
(ii) the sequences of the recombinant construct [name; describe], as additional
quantification targets, followed by (iii) the full GRCm39 primary assembly, whose
sequences were supplied as decoys. Decoy sequences absorb reads that originate from
unannotated, intronic or intergenic sequence so that they are not assigned to
similar-looking transcripts. Before indexing, sequence names across the three inputs
were checked for duplicates. Each construct feature was mapped to itself in the
transcript-to-gene table so it could be reported individually.

## Quantification

Trimmed read pairs were quantified with salmon quant v1.10.3 against this index using
selective alignment (the default mapping mode in this version), with the library type
set explicitly to ISR (inward, stranded, read 1 reverse; established unambiguously on a
test library with a compatible-fragment ratio of 1.0) for all samples, and with
sequence-specific (--seqBias) and fragment-GC (--gcBias) bias correction enabled.
Per-sample mapping rates, library types and decoy-assigned fragment counts were taken
from Salmon's meta_info.json for quality control [report range of mapping rates].

## Gene-level summarisation and normalisation

Transcript-level estimates were imported into R [version] with tximport [version]
(type "salmon", transcript versions ignored) and summarised to the gene level using the
transcript-to-gene table above, with average transcript length offsets carried into
DESeq2 [version] via DESeqDataSetFromTximport. Library-size/composition normalisation
was estimated from host (mouse) genes only: construct-derived features remained in the
count matrix, but were excluded from normalisation using the controlGenes argument of
estimateSizeFactors, so that construct RNA, which can be a large and uneven share of
reads in infected samples, would not distort the normalisation of host genes. Genes were
retained if they had at least 10 normalised counts in at least three samples (construct
features were always retained). Variance-stabilised expression values (VST,
blind = FALSE) were used for visualisation.

## Differential expression

Differential expression was tested with DESeq2 using a single-factor design in which
each interferon x virus condition was one level (design ~ group) and Mock (untreated)
was the reference. IFN-beta versus Mock and IFN-gamma versus Mock were extracted as
contrasts; [other contrasts: EBOV, MARV, EVLP, MVLP and the IFN x virus combinations
versus Mock]. P values were from the DESeq2 Wald test and adjusted for multiple testing
by the Benjamini-Hochberg procedure (padj). Reported log2 fold-changes are shrunk with
the adaptive shrinkage estimator (ashr) [confirm ashr was installed; otherwise state
"normal" shrinkage] and, where indicated, unshrunk values from the Wald test were used.
Interferon-stimulated genes (ISGs) were defined operationally as genes with padj < 0.05
and log2 fold-change > 1 (at least 2-fold increase) versus Mock, and interferon-
downregulated genes (IDGs) as padj < 0.05 and log2 fold-change < -1.

## Gene annotation and filtering

Gene symbols and biotypes were taken from the gene_name and gene_biotype attributes of
the same Ensembl release 109 GTF used for quantification. Analyses restricted to
protein-coding genes used the biotype "protein_coding".

## Assessment of the re-quantification (optional; for supplement)

To determine whether decoy-aware, bias-corrected re-quantification changed conclusions,
differential expression results were compared with those from an earlier quantification
of the same trimmed reads using an index without decoys, automatic library-type
detection and no bias correction. Log2 fold-changes were compared by Pearson correlation
and significant calls were classified as retained, lost or gained (same thresholds as
above). [Report: IFN-beta r = 0.79 shrunk / 0.77 unshrunk; IFN-gamma r = 0.85 / 0.87;
of significant genes in the original analysis, 1,084/1,200 (IFN-beta) and 1,844/2,004
(IFN-gamma) were retained when unshrunk fold-changes were compared.] Size factors
estimated on host genes did not differ systematically between conditions.

## Software

Salmon v1.10.3; cutadapt [version]; gffread [version]; R [version, e.g. 4.5.x];
tximport [version]; DESeq2 [version]; ashr [version]; readr, dplyr, ggplot2, pheatmap,
ggrepel [versions]. Analysis scripts are available at [repository / DOI]. Raw and
processed data: [GEO accession].
