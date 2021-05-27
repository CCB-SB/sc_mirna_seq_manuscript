library(data.table)
library(ggplot2)
library(ggrepel)
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
library(RColorBrewer)


darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) = names(color)
  col
}

uni_color = "#7570b3"

metadata = fread("data/figure4_input/samples.csv")

###### Cells per patient ######
metadata[, cells:=.N, by=Patient]
metadata[, Patient:=factor(Patient, levels = c(paste("P", 1:7, sep = ""), "N"))]

patient_colors = c(brewer.pal(7, "Set1"), "#d3d3d3")
names(patient_colors) = levels(metadata$Patient)

cell_numbers_p = ggplot(metadata, aes(x=Patient)) +
  geom_bar(fill=uni_color) +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size=2.5) +
  scale_y_continuous(limits=c(0,NA), expand = expansion(mult=c(0,0.15))) + 
  ylab("#Cells") + xlab("Patient") + theme_cowplot(12)

##### Sequencing stats #######
cutadapt = fread("data/figure4_input/cutadapt_stats.csv")
star = fread("data/figure4_input/mapping_overview.csv")
read_compo = fread("data/figure4_input/overall_composition.detailed.raw.csv")

sample_stats = merge(
  merge(
    merge(cutadapt, metadata, by="Sample"),
    star, by="Sample"),
    read_compo, by="Sample")
sample_stats[, read_lost_in_qc:=(r_processed-r_written)/r_processed]
sample_stats[, unmapped_reads:=(r_written-mapped_total)/r_processed]
sample_stats[, mapped_ratio:=(mapped_total)/r_processed]
sample_stats[, mirna_ratio:=miRNA/r_processed]
sample_stats[, mapped_ratio_remaining:=mapped_total/r_written]
sample_stats[, mirna_ratio_remaining:=miRNA/r_written]

###### filter cells miRNAs #####
mir_mol_count_expr = fread("data/figure4_input/mirna_quantification_dedup_raw_mirna.csv")
mir_mol_counts = colSums(mir_mol_count_expr[, 2:ncol(mir_mol_count_expr)])

mir_mol_counts_df = data.frame(count=mir_mol_counts)
mir_mol_counts_df$Sample = rownames(mir_mol_counts_df)
mir_mol_counts_df$Patient = metadata$Patient[match(rownames(mir_mol_counts_df), metadata$Sample)]
mir_mol_counts_df$Negative = ifelse(mir_mol_counts_df$Patient == "N", "Negative", "CTC")
mir_mol_counts_df$Negative2 = ifelse(mir_mol_counts_df$Patient == "N", "Negative control", "CTC")
mir_mol_counts_df$discard = ifelse(mir_mol_counts_df$count < 50, "Discarded", "Kept")
mir_mol_counts_df = setDT(mir_mol_counts_df)

mol_counts_plot = ggplot(mir_mol_counts_df, 
                         aes(x=factor(Sample, levels=mir_mol_counts_df[order(count)]$Sample),
                             y=count, fill=discard, label=Negative2)) + 
              geom_col() + 
  geom_text_repel(
    data          = mir_mol_counts_df[Patient == "N"],
    nudge_y       = 1500,
    arrow = arrow(length = unit(0.02, "npc")),
    direction     = "y",
    hjust         = 1
  ) +
  geom_vline(xintercept = sum(mir_mol_counts_df$count < 50) + 0.5, linetype="dashed", color="grey50") +
  coord_flip() +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult=c(0, 0.05)), breaks=scales::pretty_breaks(5)) +
  scale_fill_manual(values=c("Kept"=uni_color, "Discarded"="#e66101"), name="") +
  xlab("Cells") +
  ylab("Molecules") + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "none")

save_plot("results/figures/supp_fig11.pdf", mol_counts_plot, base_height = 100, base_width = 100, unit="mm")

mir_mol_source_df = copy(mir_mol_counts_df[, c("count", "Sample", "Patient", "Negative2", "discard")])
setnames(mir_mol_source_df, c("count", "Negative2", "discard"), c("Cells", "Negative", "Discard"))
fwrite(mir_mol_source_df, "results/figures/supp_fig11.txt", sep='\t')

valid_cells = names(mir_mol_counts[mir_mol_counts>=50])

###### RNA composition #######
rna_class_colors_df = fread("data/rna_colors.csv")
rna_class_colors = rna_class_colors_df$color
names(rna_class_colors) = rna_class_colors_df$RNA

rna_composition = fread("data/figure4_input/overall_composition.detailed.mapped_perc.csv")

rna_comp_df = melt(rna_composition)
rna_comp_df[, Patient:=metadata$Patient[match(Sample, metadata$Sample)]]
rna_comp_df[, variable:=factor(variable, levels=rna_comp_df[, mean(value), by="variable"][order(V1, decreasing = T)]$variable)]

# sum all types < 0.5 mean expression in category "Other"
other_rna_classes = rna_comp_df[,mean(value), by="variable"][V1 < 0.5]$variable
rna_comp_df[variable %in% other_rna_classes, variable:="other"]
rna_comp_df = rna_comp_df[, list(value=sum(value)), by=c("Sample", "variable", "Patient")]

compo_grouped = ggplot(rna_comp_df, aes(x=factor(Sample, levels=rna_comp_df[variable == "miRNA"][order(value)]$Sample),
                                             y=value, fill=variable)) + 
  facet_grid(rows=vars(Patient), scales = "free_y", space = "free", switch = "y") +
  geom_col() +
  theme_cowplot(12) + scale_y_continuous(labels = scales::percent_format(accuracy = 1, scale=1), expand = expansion(mult=c(0, 0.05))) +
  coord_flip() + xlab("Patient") +
  ylab("Read composition") +
  scale_fill_manual(values=rna_class_colors, 
                    labels=function(x) str_trim(gsub("_", " ", gsub("ensembl", "", x))), 
                    name="Class",
                    guide=guide_legend(ncol=3, title.position = "top")) +
  theme(legend.position = "bottom", axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size=9),
        panel.spacing = unit(0, "pt"),
        legend.margin = margin(-25, 0, 0, 0, "pt"))


supp_table9 = sample_stats[, c("Sample", "Patient",
                               "r_processed", "r_written",
                               "uniquely_mapped", "multimapped", "mapped_total",
                               "miRNA", "protein_coding_gene",
                               "intergenic", "lncRNA",
                               "snoRNA_ensembl",
                               "coding_exons", "rRNA_ensembl", "GtRNAdb", "piRBase",
                               "Mt_tRNA_ensembl"
)]

supp_table9[, `miRNA (% of total)`:=miRNA/r_processed * 100]
supp_table9[, `protein coding gene (% of total)`:=protein_coding_gene/r_processed * 100]
supp_table9[, `intergenic (% of total)`:=intergenic/r_processed * 100]
supp_table9[, `lncRNA (% of total)`:=lncRNA/r_processed * 100]
supp_table9[, `snoRNA (% of total)`:=snoRNA_ensembl/r_processed * 100]
supp_table9[, `other mapped reads (% of total)`:=(sample_stats$mapped_total - miRNA - protein_coding_gene - intergenic - lncRNA - snoRNA_ensembl - GtRNAdb - rRNA_ensembl - coding_exons - piRBase - Mt_tRNA_ensembl) / r_processed * 100]
supp_table9[, `Uniquely mapped to hg38 (% of total)`:=uniquely_mapped / r_processed * 100]
supp_table9[, `Multimapped to hg38 (% of total)`:=multimapped / r_processed * 100]
supp_table9[, `Unmapped to hg38 (% of total)`:=(r_written - mapped_total) / r_processed * 100]
supp_table9[, `coding exons (% of total)`:=coding_exons/r_processed * 100]
supp_table9[, `rRNA (% of total)`:=rRNA_ensembl/r_processed * 100]
supp_table9[, `GtRNAdb (% of total)`:=GtRNAdb/r_processed * 100]
supp_table9[, `piRBase (% of total)`:=piRBase/r_processed * 100]
supp_table9[, `Mt tRNA (% of total)`:=Mt_tRNA_ensembl/r_processed * 100]

setnames(supp_table9,
         c("r_processed", "r_written"),
         c("Total reads", "Cleaned reads"))
supp_table9 = supp_table9[, c("Sample", "Patient",
                              "Total reads", "Cleaned reads", 
                              "Uniquely mapped to hg38 (% of total)", "Multimapped to hg38 (% of total)", "Unmapped to hg38 (% of total)",
                              "miRNA (% of total)", "protein coding gene (% of total)", "intergenic (% of total)",
                              "lncRNA (% of total)", "snoRNA (% of total)", "GtRNAdb (% of total)", "Mt tRNA (% of total)", "rRNA (% of total)", "piRBase (% of total)", "coding exons (% of total)", "other mapped reads (% of total)")]

fwrite(supp_table9, "results/tables/supp_data6.csv", sep='\t')

###### Reads per UMI #######
expr = fread("data/figure4_input/mirna_quantification_raw_mirna.csv")
expr.dedup = fread("data/figure4_input/mirna_quantification_dedup_raw_mirna.csv")

#UMI to count ratio
dedup_expr = expr.dedup
raw_expr_umi_samples = expr[, colnames(dedup_expr), with=F]

counts_per_umi = raw_expr_umi_samples[,2:ncol(raw_expr_umi_samples)] / dedup_expr[,2:ncol(dedup_expr)]
counts_per_umi[is.na(counts_per_umi) | counts_per_umi == Inf] = 0
mirs_w_expr = rowSums(counts_per_umi) > 0
counts_per_umi = counts_per_umi[mirs_w_expr]
counts_per_umi_plotdf = melt(counts_per_umi)
counts_per_umi_plotdf = counts_per_umi_plotdf[value > 0]
counts_per_umi_plotdf[, Patient:=metadata$Patient[match(variable, metadata$Sample)]]

counts_per_umi_p = ggplot(counts_per_umi_plotdf, aes(x=factor(variable, levels=counts_per_umi_plotdf[,median(value), by="variable"][order(V1)]$variable), y=value, color=Patient, fill=Patient)) +
  facet_grid(rows=vars(Patient), scales = "free_y", space = "free", switch = "y") +
  geom_boxplot(outlier.alpha = 1) +
  scale_fill_manual(values=patient_colors) +
  scale_color_manual(values=darken(patient_colors)) +
  theme_cowplot(12) + 
  scale_y_continuous(breaks=scales::pretty_breaks(10)) +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size=9),
        panel.spacing = unit(0, "pt"),
        legend.position = "none") +
  xlab("Patient") + ylab("Reads per UMI")

######### Most expressed miRNAs #######
expr_perc = apply(expr.dedup[,..valid_cells], 2, function(e) e / sum(e) * 100)
rownames(expr_perc) = expr.dedup$miRNA

ranks = apply(expr_perc, 2, function(e) rank(-e))

mean_ranks = rowMeans(ranks)
mean_expr = rowMeans(expr_perc)

# get top 20 most expressed mirnas on average
top_mirnas = names(mean_expr[order(-mean_expr)[1:20]])

top_mirna_df = setDT(reshape2::melt(expr_perc[top_mirnas,2:ncol(expr_perc)], value.name="Expression"))

top_mirna_log_p = ggplot(top_mirna_df, aes(x=factor(Var1, levels=rev(levels(Var1))), y=Expression)) + 
  geom_half_boxplot(outlier.alpha = 0, color=darken(uni_color), fill=uni_color) +
  geom_half_dotplot(binwidth = 0.05, method="histodot", dotsize = 1, color=darken(uni_color)) +
  ylab("Expression (%)") + xlab("") +
  scale_y_sqrt(expand=expansion(mult=c(0,0)), limits=c(0, NA),
               breaks=c(0, 1, 2, 5, 10, 20, 30, 50)) +
  coord_flip() +
  theme_cowplot(12)


########## miRNA expression stats ##########
sd_expr = apply(expr_perc, 1, sd)[names(mean_expr)]
median_expr = apply(expr_perc, 1, median)[names(mean_expr)]
expressed_cells = apply(expr_perc, 1, function(e) sum(e > 0))[names(mean_expr)]

df = data.table(miRNA=names(mean_expr),
                "Mean expression (%)"=mean_expr,
                "Standard deviation (%)"=sd_expr,
                "Median expression (%)"=median_expr,
                "#Expressed cells"=expressed_cells)
fwrite(df[order(`Mean expression (%)`, decreasing = TRUE)], "results/tables/supp_data7.csv", sep='\t')


####### Most found miRNAs #######
detected_mirs_in_cells_df = data.table(miRNA=rownames(expr_perc),
                                   cells=rowSums(expr_perc > 0))[order(cells, decreasing = TRUE)][, cells_p:=cells/ncol(expr_perc)]


detected_mirs_in_cells = ggplot(detected_mirs_in_cells_df[cells_p >= 0.9],
                                aes(x=factor(miRNA, levels=rev(detected_mirs_in_cells_df$miRNA)), y=cells_p)) + 
  geom_col(fill="#7570b3", width = 0.75) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), labels = scales::percent_format(), limits = c(0, 1)) +
  coord_flip() +
  ylab("Detected in cells") +
  xlab("") +
  theme_cowplot(12) +
  theme(plot.margin = margin(5, 15, 1, 1))

save_plot("results/figures/supp_fig12.pdf", detected_mirs_in_cells, base_height = 100, base_width = 100, unit="mm")

detected_mirs_source_df = copy(detected_mirs_in_cells_df[cells_p >= 0.9])
detected_mirs_source_df[, cells_p:=cells_p*100]
setnames(detected_mirs_source_df, c("cells_p"), "Detected in cells (%)")

fwrite(detected_mirs_source_df, "results/figures/supp_fig12.txt", sep='\t')

####### Most variable mirnas ########
expr_mat = as.matrix(expr.dedup[, ..valid_cells])
rownames(expr_mat) = expr.dedup$miRNA

expr_mat_norm = log2(scale(expr_mat, center=F, scale=colSums(expr_mat)) * 1e4 + 1)

expr_mat.filtered = expr_mat[, c(TRUE, !metadata$Patient[match(colnames(expr_mat)[2:ncol(expr_mat)], metadata$Sample)] %in% c("N", "P7", "P6"))]
expr_mat_norm.filtered = expr_mat_norm[, c(TRUE, !metadata$Patient[match(colnames(expr_mat_norm)[2:ncol(expr_mat_norm)], metadata$Sample)] %in% c("N", "P7", "P6"))]

seu.filtered = CreateSeuratObject(counts=expr_mat.filtered,
                                  meta.data=data.frame(metadata[match(colnames(expr_mat.filtered), Sample)], row.names=colnames(expr_mat.filtered)))
seu.filtered = NormalizeData(seu.filtered, normalization.method = "LogNormalize", scale.factor=10000)
seu.filtered = FindVariableFeatures(seu.filtered, selection.method = "vst", nfeatures = 30)

res = data.table(reshape2::melt(expr_mat_norm.filtered, id.vars = "miRNA", variable.name = "Sample", value.name = "Count"))
colnames(res) = c("miRNA", "Sample", "Count")
res[, Patient:=metadata$Patient[match(Sample, metadata$Sample)]]

most_variable_vst_fts = head(VariableFeatures(seu.filtered), 6)
most_var_mirs_p = ggplot(res[miRNA %in% most_variable_vst_fts], aes(x=Count, y=Patient, fill=Patient, color=Patient)) +
  facet_wrap(~miRNA) +
  stat_density_ridges(jittered_points = TRUE, quantile_lines = TRUE, position="points_sina", point_size=0.1, size=0.1, alpha=0.5, point_alpha=0.5) +
  scale_fill_manual(values=patient_colors) + scale_color_manual(values=darken(patient_colors)) +
  theme_cowplot(12) +
  xlab(bquote(Normalized~log[2]~expression)) + 
  theme(legend.position="none",
        strip.text.x = element_text(size=9))

########### Seurat embedding #########

seu = CreateSeuratObject(counts=expr_mat,
                         meta.data=data.frame(metadata[match(colnames(expr_mat), Sample)], row.names=colnames(expr_mat)))

seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 200)


all.genes = rownames(seu)
seu = ScaleData(seu, features = all.genes)

seu = RunPCA(seu, features = VariableFeatures(object = seu))
seu = RunUMAP(seu, dims = 1:15, n.neighbors = 15, min.dist = 0.1)
seu = FindNeighbors(seu, dims = 1:15)
seu = FindClusters(seu, resolution = 1.4)

seu_umap = data.table(seu@reductions$umap@cell.embeddings, keep.rownames = TRUE)
seu_umap[, Sex:=metadata$Sex[match(rn, metadata$Sample)]]
seu_umap[, Patient:=metadata$Patient[match(rn, metadata$Sample)]]
seu_umap[, Cluster:=seu@meta.data$RNA_snn_res.1.4]

seu_umap_plot_patient = ggplot(seu_umap, aes(x=UMAP_1, y=UMAP_2, color=Patient)) +
  geom_point(size=1) +
  scale_color_manual(values=patient_colors) +
  theme_cowplot(12) + 
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")

seu_umap_plot_cluster = ggplot(seu_umap, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point(size=1) +
  scale_color_brewer(type="qual", palette = 7) +
  theme_cowplot(12) + 
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")


seu13 = CreateSeuratObject(counts=expr_mat,
                           meta.data=data.frame(metadata[match(colnames(expr_mat), Sample)], row.names=colnames(expr_mat)))

seu13 = seu13[, seu13$orig.ident %in% c("18-043", "19-032", "19-231")]
seu13 = NormalizeData(seu13, normalization.method = "LogNormalize", scale.factor = 10000)

seu13 = FindVariableFeatures(seu13, selection.method = "vst", nfeatures = 200)


all.genes = rownames(seu13)
seu13 = ScaleData(seu13, features = all.genes)

seu13 = RunPCA(seu13, features = VariableFeatures(object = seu13), npcs = 30)
seu13 = RunUMAP(seu13, dims = 1:15, n.neighbors = 15, min.dist = 0.1)
seu13 = FindNeighbors(seu13, dims = 1:15)
seu13 = FindClusters(seu13, resolution = 1.05)

seu13_umap = data.table(seu13@reductions$umap@cell.embeddings, keep.rownames = TRUE)
seu13_umap[, Sex:=metadata$Sex[match(rn, metadata$Sample)]]
seu13_umap[, Patient:=metadata$Patient[match(rn, metadata$Sample)]]
seu13_umap[, Cluster:=seu13@meta.data$RNA_snn_res.1.05]

seu13_umap_plot_patient = ggplot(seu13_umap, aes(x=UMAP_1, y=UMAP_2, color=Patient)) +
  geom_point(size=1) +
  scale_color_manual(values=patient_colors) +
  theme_cowplot(12) + 
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")

seu13_umap_plot_cluster = ggplot(seu13_umap, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point(size=1) +
  scale_color_brewer(type="qual", palette = 7) +
  theme_cowplot(12) + 
  theme(legend.position = "right", axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")


###### figure 3 ####
fig4ab = plot_grid(NULL, cell_numbers_p, labels=c("a", "b"))
fig4cd = plot_grid(compo_grouped, counts_per_umi_p, labels=c("c", "d"))
fig4ef = plot_grid(top_mirna_log_p, most_var_mirs_p, labels=c("e", "f"))
fig4g = plot_grid(plot_grid(seu_umap_plot_patient, seu_umap_plot_cluster, nrow=1), labels="g")
fig4h = plot_grid(plot_grid(seu13_umap_plot_patient, seu13_umap_plot_cluster, nrow=1), labels="h")

fig4 = plot_grid(fig4ab, fig4cd, fig4ef, plot_grid(fig4g, fig4h, nrow=1), ncol=1, rel_heights = c(1, 2, 1.5, 1))
save_plot("results/figures/fig4.pdf", fig4, base_height = 297, base_width = 210, unit="mm")

fwrite(metadata[, list("#Cells"=.N), by="Patient"][order(`Patient`)], "results/figures/fig4b.txt", sep='\t')

setnames(rna_comp_df, c("value", "variable"), c("Reads (%)", "Class"))
fwrite(rna_comp_df, "results/figures/fig4c.txt", sep='\t')

setnames(counts_per_umi_plotdf, c("variable", "value"), c("Sample", "Reads per UMI"))
fwrite(counts_per_umi_plotdf, "results/figures/fig4d.txt", sep='\t')

setnames(top_mirna_df, c("Var1", "Var2"), c("miRNA", "Sample"))
fwrite(top_mirna_df, "results/figures/fig4e.txt", sep='\t')

most_var_mir_expr = copy(res[miRNA %in% most_variable_vst_fts])
setnames(most_var_mir_expr, "Count", "Normalized log2 expression")
fwrite(most_var_mir_expr, "results/figures/fig4f.txt", sep='\t')

setnames(seu_umap, c("rn", "UMAP_1", "UMAP_2"), c("Sample", "UMAP 1", "UMAP 2"))
fwrite(seu_umap[, c("Sample", "UMAP 1", "UMAP 2", "Patient", "Cluster")], "results/figures/fig4g.txt", sep='\t')

setnames(seu13_umap, c("rn", "UMAP_1", "UMAP_2"), c("Sample", "UMAP 1", "UMAP 2"))
fwrite(seu13_umap[, c("Sample", "UMAP 1", "UMAP 2", "Patient", "Cluster")], "results/figures/fig4h.txt", sep='\t')

###### mieaa input #####
mirna_exp_ordered = mean_expr[order(mean_expr, decreasing = T)]
mirna_exp_ordered = mirna_exp_ordered[mirna_exp_ordered > 0]

dir.create("results/mieaa_input", showWarnings = F, recursive = T)
writeLines(names(mirna_exp_ordered), "results/mieaa_input/CTC_SCLC_gsea.txt")
