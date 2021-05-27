library(data.table)
library(ggplot2)
library(cowplot)
library(ggsci)
library(viridisLite)
library(gghalves)
library(ggridges)
library(pbapply)
library(extrafont)
library(Seurat)
library(ComplexUpset)
library(stringr)

library(RhpcBLASctl)

blas_p = blas_get_num_procs()
blas_set_num_threads(1)
omp_p = omp_get_num_procs()
omp_set_num_threads(1)


darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) = names(color)
  col
}

cellline_colors = c("HepG2"="#e6194B",    
                    "HT29"="#ffe119",
                    "A549"="#4363d8",
                    "Jurkat"="#42d4f4",
                    "REH"="#469990",
                    "KG1"="#dcbeff",
                    "BJ"="#800000",
                    "THP-1"="#000075")

cellline_dark_colors = darken(cellline_colors)

method_colors = c("Sandberg"="#1b9e77",
                  "CleanTag"="#d95f02",
                  "Giraldez"="#7570b3",
                  "Jensen"="#e7298a",
                  "CATS"="#e6ab02",
                  "Faridani"="#66a61e",
                  "Shore"="#4daf4a")

adapter_3p_colors = c("CleanTag"="#fdc086",
                      "Sandberg"="#7fc97f",
                      "4N"="#beaed4",
                      "Rand"="#d6604d",
                      "CATS"="#ffff99")

adapter_5p_colors = c("C3"="#d9d9d9",
                      "Sandberg"="#7fc97f",
                      "4N"="#beaed4",
                      "CleanTag"="#fdc086",
                      "CleanTag UMI6"="#b2182b",
                      "Block"="#fb8072",
                      "CATS"="#ffff99")

top8_protocols = names(cellline_colors)

metadata = fread("data/figure3_input/samples.csv")
metadata[, `Sample name`:=factor(`Sample name`, levels=metadata[order(Group, `Sample name`)]$`Sample name`)]
cutadapt = fread("data/figure3_input/cutadapt_stats.csv")
star = fread("data/figure3_input/mapping_overview.csv")
read_compo = fread("data/figure3_input/overall_composition.detailed.raw.csv")
cutadapt_len_stats = fread("data/figure3_input/overall.too_short.len_stats.csv")
adapter_dimers = colSums(cutadapt_len_stats[Length <= 3][,2:ncol(cutadapt_len_stats)])

sample_stats = merge(
    merge(
      merge(cutadapt, metadata, by.y = "Sample name", by.x="Sample"),
      star, by="Sample"),
  read_compo, by="Sample")
sample_stats[, adapter_dimers_ratio:=(adapter_dimers[match(Sample, names(adapter_dimers))])/r_processed]
sample_stats[, read_lost_in_qc:=(r_processed-r_written)/r_processed-adapter_dimers_ratio]
sample_stats[, unmapped_reads:=(r_written-mapped_total)/r_processed]
sample_stats[, mapped_ratio:=(mapped_total - miRNA)/r_processed]
sample_stats[, mirna_ratio:=miRNA/r_processed]

######### sample stats table #############
supp_table8 = sample_stats[, c("Sample", "Group",
                               "r_processed", "r_written",
                               "uniquely_mapped", "multimapped", "mapped_total",
                               "miRNA", "protein_coding_gene",
                               "intergenic", "lncRNA",
                               "snoRNA_ensembl",
                               "coding_exons", "rRNA_ensembl", "GtRNAdb", "piRBase"
)]

supp_table8[, `miRNA (% of total)`:=miRNA/r_processed * 100]
supp_table8[, `protein coding gene (% of total)`:=protein_coding_gene/r_processed * 100]
supp_table8[, `intergenic (% of total)`:=intergenic/r_processed * 100]
supp_table8[, `lncRNA (% of total)`:=lncRNA/r_processed * 100]
supp_table8[, `snoRNA (% of total)`:=snoRNA_ensembl/r_processed * 100]
supp_table8[, `other mapped reads (% of total)`:=(sample_stats$mapped_total - miRNA - protein_coding_gene - intergenic - lncRNA - snoRNA_ensembl - GtRNAdb - rRNA_ensembl - coding_exons - piRBase) / r_processed * 100]
supp_table8[, `Uniquely mapped to hg38 (% of total)`:=uniquely_mapped / r_processed * 100]
supp_table8[, `Multimapped to hg38 (% of total)`:=multimapped / r_processed * 100]
supp_table8[, `Unmapped to hg38 (% of total)`:=(r_written - mapped_total) / r_processed * 100]
supp_table8[, `coding exons (% of total)`:=coding_exons/r_processed * 100]
supp_table8[, `rRNA (% of total)`:=rRNA_ensembl/r_processed * 100]
supp_table8[, `GtRNAdb (% of total)`:=GtRNAdb/r_processed * 100]
supp_table8[, `piRBase (% of total)`:=piRBase/r_processed * 100]

setnames(supp_table8,
         c("Group", "r_processed", "r_written"),
         c("Protocol", "Total reads", "Cleaned reads"))
supp_table8 = supp_table8[, c("Sample", "Protocol",
                              "Total reads", "Cleaned reads", 
                              "Uniquely mapped to hg38 (% of total)", "Multimapped to hg38 (% of total)", "Unmapped to hg38 (% of total)",
                              "miRNA (% of total)", "protein coding gene (% of total)", "intergenic (% of total)",
                              "lncRNA (% of total)", "snoRNA (% of total)", "GtRNAdb (% of total)", "rRNA (% of total)", "piRBase (% of total)", "coding exons (% of total)", "other mapped reads (% of total)")]

fwrite(supp_table8, "results/tables/supp_data5.csv", sep='\t')

sample_stats_plot_df = melt(sample_stats[, c("Sample", "adapter_dimers_ratio", "read_lost_in_qc", "unmapped_reads", "mapped_ratio", "mirna_ratio")], id.vars="Sample")
sample_stats_plot_df[, Group:=metadata$Group[match(Sample, metadata$`Sample name`)]]

sample_stats_mean_plot_df = sample_stats_plot_df[, list(value=mean(value), sd=sd(value)), by=c("variable", "Group")]
sample_stats_mean_plot_df = sample_stats_mean_plot_df[order(variable, decreasing = TRUE)]
sample_stats_mean_plot_df[, ypos:=cumsum(value), by="Group"]

sample_stats_per_group_with_errorbars = ggplot(sample_stats_mean_plot_df, aes(x=factor(Group, levels=sample_stats_mean_plot_df[variable == "mirna_ratio"][order(value)]$Group),
                                                                              y=value, fill=variable)) +
  geom_col() +
  geom_errorbar(aes(ymax=ypos-sd, ymin=ypos-sd, color=variable), position="identity", width=0.5, show.legend = F) +
  geom_linerange(aes(ymin = ypos-sd, ymax = ypos, color=variable), show.legend=F) + 
  xlab("") + 
  scale_y_continuous(labels=scales::percent_format(), name="Total sequenced reads", expand = expansion(mult = c(0,0))) +
  coord_flip() +
  scale_fill_manual(values=c(read_lost_in_qc="#d3d3d3", adapter_dimers_ratio="#adadad", unmapped_reads="#3B4992FF", mapped_ratio="#008B45FF", mirna_ratio="#c51b8a"),
                    labels=c(read_lost_in_qc="Lost in QC", adapter_dimers_ratio="Adapter dimers", unmapped_reads="Unmapped", mapped_ratio="Mapped", mirna_ratio="miRNAs"),
                    name="",
                    guide=guide_legend(nrow=2)) +
  scale_color_manual(values=c(read_lost_in_qc=darken("#d3d3d3", 1.5), adapter_dimers_ratio=darken("#adadad", 1.5), unmapped_reads=darken("#3B4992FF"), mapped_ratio=darken("#008B45FF"), mirna_ratio=darken("#c51b8a")),
                     labels=c(read_lost_in_qc="Lost in QC", adapter_dimers_ratio="Adapter dimers", unmapped_reads="Unmapped", mapped_ratio="Mapped", mirna_ratio="miRNAs"),
                     name="",
                     guide=guide_legend(nrow=2)) +
  theme_cowplot(11) + theme(legend.position="bottom")


mir_expr = fread("data/figure3_input/mirna_quantification_rpmmm_norm.csv")
mir_expr[, miRNA:=paste(precursor, miRNA, sep='|')]
mir_expr[, precursor:=NULL]

detected_mirs = data.frame(num_mirnas=colSums(mir_expr[, 2:ncol(mir_expr)] > 0))
detected_mirs$Group = metadata$Group[match(rownames(detected_mirs), metadata$`Sample name`)]
detected_mirs = setDT(detected_mirs)

detected_mirs_p = ggplot(detected_mirs, aes(x=factor(Group, levels=detected_mirs[, mean(num_mirnas), by="Group"][order(V1)]$Group), y=num_mirnas, color=Group, fill=Group)) + 
  geom_half_point() +
  geom_half_boxplot(outlier.colour = NA) +
  theme_cowplot(11) +
  scale_y_continuous(breaks=scales::pretty_breaks(10)) +
  scale_fill_manual(values=cellline_colors) + scale_color_manual(values=cellline_dark_colors) +
  coord_flip() +
  xlab("") + ylab("Detected miRNAs") + theme(legend.position = "none")


expr = fread("data/figure3_input/read_counts.csv")

expr_mat = as.matrix(expr[, 2:ncol(expr)])
rownames(expr_mat) = expr$Geneid
seu = CreateSeuratObject(counts=expr_mat,
                         meta.data=data.frame(metadata[match(colnames(expr_mat), `Sample name`)], row.names=colnames(expr_mat)))

seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes = rownames(seu)
seu = ScaleData(seu, features = all.genes)

seu = RunPCA(seu, features = VariableFeatures(object = seu), npcs = 47)
seu = RunUMAP(seu, dims = 1:20, n.neighbors = 12, min.dist = 0.1,  seed.use=42)

seu_umap = data.table(seu@reductions$umap@cell.embeddings, keep.rownames = TRUE)
seu_umap[, Protocol:=metadata$Group[match(rn, metadata$Sample)]]

seu_umap_plot = ggplot(seu_umap, aes(x=UMAP_1, y=UMAP_2, color=Protocol, label=rn)) +
  geom_point(size=1) +
  scale_color_manual(values=cellline_colors, guide=guide_legend(ncol=2), name="Cell line") +
  theme_cowplot(11) + 
  theme(legend.position = c(0.01,0.75), axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")


expressed_mirs_per_protocol = lapply(unique(metadata$Group), function(x) {
  samples = as.character(metadata[Group == x]$`Sample name`)
  mir_expr$miRNA[rowSums(mir_expr[, ..samples] > 0) > 0]
})

names(expressed_mirs_per_protocol) = unique(metadata$Group)


qlist_intersect = lapply(names(expressed_mirs_per_protocol), function(n){
  upset_query(intersect = n, fill="#fd8d3c", color="#fd8d3c", only_components=c('Common\nmiRNAs', 'intersections_matrix'))
})

upset_p = upset(UpSetR::fromList(expressed_mirs_per_protocol), intersect=names(expressed_mirs_per_protocol),
      base_annotations = list(`Common\nmiRNAs`=intersection_size(counts=TRUE, text=list(angle=90, vjust=0.5, hjust=1, size=2), bar_number_threshold = 1)),
      queries=c(qlist_intersect),
      matrix=intersection_matrix(
        geom=geom_point(
          size=2
        )
      ),
      themes=upset_default_themes(text=element_text(size=9)),
      name="",
      sort_sets="ascending",
      n_intersections=20,
      height_ratio = 1.8,
      width_ratio = 0.3
)

save_plot("results/figures/fig3f.pdf", upset_p, base_width = 100, base_height = 50, unit="mm")

### gsea according to decreasing avg. expression
expr.dedup = fread("data/figure3_input/mirna_quantification_dedup_raw_mirna.csv")

expr_perc = apply(expr.dedup[, 2:ncol(expr.dedup)], 2, function(e) e / sum(e) * 100)
rownames(expr_perc) = expr.dedup$miRNA

expr_perc_df = as.data.table(reshape2::melt(expr_perc))
expr_perc_df[, Group:=metadata$Group[match(Var2, metadata$`Sample name`)]]
mean_expr = expr_perc_df[, mean(value), by=c("Var1", "Group")][order(V1, decreasing=TRUE)]

dir.create("results/mieaa_input", showWarnings = F, recursive = T)
for(c in unique(mean_expr$Group)){
  print(c)
  writeLines(as.character(mean_expr[Group == c & V1 > 0]$Var1), sprintf("results/mieaa_input/%s_gsea.txt", c))
}

fig3ab = plot_grid(sample_stats_per_group_with_errorbars, detected_mirs_p, labels = c("a", "b"))

fig3ctof = plot_grid(
  seu_umap_plot,
  NULL,
  ncol=2,
  labels=c("c", "d")
)
fig3atof = plot_grid(fig3ab, fig3ctof, ncol=1, rel_heights = c(1,1))
save_plot("results/figures/fig3.pdf", fig3atof, base_height = 100, base_width = 210, unit="mm")

setnames(sample_stats_mean_plot_df, c("Group", "value"), c("Protocol", "mean"))
fwrite(sample_stats_mean_plot_df[, c("variable", "Protocol", "mean", "sd"), with=F], "results/figures/fig3a.txt", sep='\t')

setnames(detected_mirs, "Group", "Protocol")
fwrite(detected_mirs, "results/figures/fig3b.txt", sep='\t')

setnames(seu_umap, c("rn", "UMAP_1", "UMAP_2"), c("Sample", "UMAP 1", "UMAP 2"))
fwrite(seu_umap, "results/figures/fig3c.txt", sep='\t')

fwrite(UpSetR::fromList(expressed_mirs_per_protocol), "results/figures/fig3d.txt", sep='\t')
