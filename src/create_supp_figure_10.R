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

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) = names(color)
  col
}

protocol_colors = c("SB"="#2b98d9",    # top
                    "SB_C3"="#d0d1e6",
                    "SB_CL"="#0570b0",
                    "SB_4N"="#74a9cf",  # top
                    "CL"="#cc4c02",              # top
                    "CL_SB"="#662506",
                    "CL_4N"="#993404",
                    "CL_Rand"="#ec7014",
                    "CL_C3"="#fe9929",
                    "CL_UMI6"="#fec44f",
                    "CL_16C"="#fee391",
                    "CL_Block"="#f2c813",
                    "4N"="#19b576",             # top
                    "4N_C3"="#238443",           # top
                    "4N_CL"="#59de5a",           #top
                    "SBN"="#b02cad",             #top
                    "SBN_CL"="#810f7c",         #top
                    "SBN_4N"="#8c6bb1",
                    "CATS"="#35978f")

protocol_dark_colors = darken(protocol_colors)

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

top8_protocols = c("CL", "SB_4N", "4N", "4N_C3", "SB", "SBN", "4N_CL", "SBN_CL")

metadata = fread("data/supp_fig10_input/samples_sce.csv")
metadata[, `Sample name`:=factor(`Sample name`, levels=metadata[order(Group, `Sample name`)]$`Sample name`)]
cutadapt = fread("data/supp_fig10_input/cutadapt_stats.csv")
star = fread("data/supp_fig10_input/mapping_overview.csv")
read_compo = fread("data/supp_fig10_input/overall_composition.detailed.raw.csv")
cutadapt_len_stats = fread("data/supp_fig10_input/overall.too_short.len_stats.csv")
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

metadata_me = fread("data/supp_fig7_input/samples_me.csv")
metadata_me[, `Sample name`:=factor(`Sample name`, levels=metadata_me[order(Group, `Sample name`)]$`Sample name`)]
cutadapt_me = fread("data/supp_fig7_input/cutadapt_stats.csv")
star_me = fread("data/supp_fig7_input/mapping_overview.csv")
read_compo_me = fread("data/supp_fig7_input/overall_composition.detailed.raw.csv")
cutadapt_len_stats_me = fread("data/supp_fig7_input/overall.too_short.len_stats.csv")
adapter_dimers_me = colSums(cutadapt_len_stats_me[Length <= 3][,2:ncol(cutadapt_len_stats_me)])


sample_stats_me = merge(
  merge(
    merge(cutadapt_me, metadata_me, by.y = "Sample name", by.x="Sample"),
    star_me, by="Sample"),
  read_compo_me, by="Sample")
sample_stats_me[, adapter_dimers_ratio:=(adapter_dimers_me[match(Sample, names(adapter_dimers_me))])/r_processed]
sample_stats_me[, read_lost_in_qc:=(r_processed-r_written)/r_processed-adapter_dimers_ratio]
sample_stats_me[, unmapped_reads:=(r_written-mapped_total)/r_processed]
sample_stats_me[, mapped_ratio:=(mapped_total - miRNA)/r_processed]
sample_stats_me[, mirna_ratio:=miRNA/r_processed]

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

sample_stats_per_group_legend = get_legend(sample_stats_per_group_with_errorbars + theme(legend.box.margin = margin(0, 0, 0, 12), legend.justification="center", legend.box.just = "bottom"))

adapter_dimers_stage12 = data.table(stage1=na.omit(sample_stats_me[match(sample_stats$Sample, Sample)]$adapter_dimers_ratio),
                               stage2=sample_stats[sample_stats$Sample %in% sample_stats_me$Sample]$adapter_dimers_ratio,
                               Sample=sample_stats[sample_stats$Sample %in% sample_stats_me$Sample]$Sample)
adapter_dimers_stage12[, Group:=metadata$Group[match(Sample, metadata$`Sample name`)]]

adapter_dimers_stage12_box = data.table(Stage=c(rep("1", sum(sample_stats_me$Group %in% top8_protocols)),
                                                rep("2", sum(sample_stats$Group %in% top8_protocols))),
                                        Reads=c(sample_stats_me[sample_stats_me$Group %in% top8_protocols]$adapter_dimers_ratio,
                                                sample_stats[sample_stats$Group %in% top8_protocols]$adapter_dimers_ratio
                                        ),
                                        Sample=c(sample_stats_me[sample_stats_me$Group %in% top8_protocols]$Sample,
                                                 sample_stats[sample_stats$Group %in% top8_protocols]$Sample),
                                        Protocol=c(sample_stats_me[sample_stats_me$Group %in% top8_protocols]$Group,
                                                   sample_stats[sample_stats$Group %in% top8_protocols]$Group
                                        )
)


adapter_dimers_stage12_boxp = ggplot(adapter_dimers_stage12_box, aes(x=factor(Protocol, adapter_dimers_stage12_box[, median(Reads), by="Protocol"]$Protocol), y=Reads, color=Stage, fill=Stage)) + 
  geom_half_boxplot() + geom_half_point() +
  scale_color_manual(values=darken(ggsci::pal_npg()(2))) +
  scale_fill_npg() +
  coord_flip() +
  scale_y_continuous(labels=scales::label_percent(accuracy = 0.1), breaks=scales::pretty_breaks(6), expand = expansion(mult=c(NA,0)), limits = c(0, 1)) +
  ylab("Adapter dimers") +
  xlab("") +
  theme_cowplot(11) +
  theme(legend.position = c(0.8, 0.9), axis.text.x = element_text(size=8))

mir_expr = fread("data/supp_fig10_input/mirna_quantification_rpmmm_norm.csv")
mir_expr[, miRNA:=paste(precursor, miRNA, sep='|')]
mir_expr[, precursor:=NULL]

detected_mirs = data.frame(num_mirnas=colSums(mir_expr[, 2:ncol(mir_expr)] > 0))
detected_mirs$Group = metadata$Group[match(rownames(detected_mirs), metadata$`Sample name`)]
detected_mirs = setDT(detected_mirs)

detected_mirs_p = ggplot(detected_mirs, aes(x=factor(Group, levels=detected_mirs[, mean(num_mirnas), by="Group"][order(V1)]$Group), y=num_mirnas, color=Group, fill=Group)) + 
  geom_half_point() +
  geom_half_boxplot(outlier.colour = NA) +
  theme_cowplot(11) +
  scale_y_continuous(breaks=scales::pretty_breaks(5)) +
  scale_fill_manual(values=protocol_colors, name="Protocol") + scale_color_manual(values=protocol_dark_colors, name="Protocol") +
  coord_flip() +
  xlab("") + ylab("Detected miRNAs")

detected_mirs_legend = get_legend(detected_mirs_p + theme(legend.box.margin = margin(0, 0, 0, 12), legend.justification="center", legend.box.just = "bottom", legend.direction = "horizontal"))

expr_norm_full = fread("data/supp_fig10_input/read_counts.csv")
expr_norm_full[, 2:ncol(expr_norm_full)] = as.data.table(as.matrix(expr_norm_full[,2:ncol(expr_norm_full)]) %*% diag(1/sample_stats[match(colnames(expr_norm_full)[2:ncol(expr_norm_full)], Sample)]$mapped_total) * 1e6)

# only keep genes that are present in at least 12 samples (~2 protocols)
expr = expr_norm_full[rowSums(expr_norm_full[,2:ncol(expr_norm_full)] > 0) >= 12]

### for every protocol - intra protocol distance and inter protocol distance  ####
dmat_top8 = as.matrix(dist(log2(t(as.matrix(expr[, 2:ncol(expr)]) + 1)), method="euclidean"))

rep_dist = lapply(top8_protocols, function(x) {
  samples = as.character(metadata[Group == x]$`Sample name`)
  rep_dist = dmat_top8[samples,samples]
  rep_dist[upper.tri(rep_dist)]
})
names(rep_dist) = top8_protocols

protocol_dist = as.data.table(melt(rep_dist))
protocol_dist[, Protocol:="Same"]

other_dist = lapply(top8_protocols, function(x) {
  samples = as.character(metadata[Group == x]$`Sample name`)
  other_dist = dmat_top8[samples,!colnames(dmat_top8) %in% samples]
  other_dist[upper.tri(other_dist)]
})
names(other_dist) = top8_protocols

dist_same_diff_protocols = as.data.table(rbind(protocol_dist, as.data.table(melt(other_dist))[, Protocol:="Different"]))
dist_same_diff_protocols[, Protocol:=factor(Protocol, levels=c("Same", "Different"))]
dist_same_diff_protocols[, mean_dist:=mean(value), by=c("L1", "Protocol")]
dist_same_diff_protocols[, L1:=factor(L1, levels=unique(dist_same_diff_protocols[Protocol == "Same"][order(mean_dist)]$L1))]

distance_plot = ggplot(dist_same_diff_protocols, aes(x=factor(L1, levels = dist_same_diff_protocols[Protocol == "Same", median(mean_dist), by="L1"][order(V1, decreasing = TRUE)]$L1), y=value, fill=Protocol, color=Protocol)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_point(size=0.5) +
  scale_fill_manual(values=c(Same="#41ab5d", Different="#d95f02")) +
  scale_color_manual(values=c(Same=darken("#41ab5d"), Different=darken("#d95f02"))) +
  coord_flip() +
  theme_cowplot(11) + theme(legend.position = "right") +
  scale_y_continuous(breaks=scales::pretty_breaks(5)) +
  ylab("Euclidean distance") + xlab("")


expr = fread("data/supp_fig10_input/read_counts.csv")

expr_mat = as.matrix(expr[, 2:ncol(expr)])
rownames(expr_mat) = expr$Geneid
seu = CreateSeuratObject(counts=expr_mat,
                         meta.data=data.frame(metadata[match(colnames(expr_mat), `Sample name`)], row.names=colnames(expr_mat)))

seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

all.genes = rownames(seu)
seu = ScaleData(seu, features = all.genes)

seu = RunPCA(seu, features = VariableFeatures(object = seu), npcs = 23)
seu = RunUMAP(seu, dims = 1:15, n.neighbors = 9)

seu_umap = data.table(seu@reductions$umap@cell.embeddings, keep.rownames = TRUE)
seu_umap[, Protocol:=metadata$Group[match(rn, metadata$Sample)]]

seu_umap_plot = ggplot(seu_umap, aes(x=UMAP_1, y=UMAP_2, color=Protocol)) +
  geom_point(size=1) +
  scale_color_manual(values=protocol_colors) +
  theme_cowplot(11) + 
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")


mir_detected_in_exps = data.table(miRNA=mir_expr$miRNA, experiments=rowMeans(mir_expr[, 2:ncol(mir_expr)] > 0))[order(experiments, decreasing = TRUE)]

# only keep one miRNA/precursor combination for this plot to remove redundancy
detected_mirs_in_expr_df = mir_detected_in_exps[, miRNA2:=gsub("^.*\\|", "", miRNA)][!duplicated(miRNA2)]
detected_mirs_in_exp = ggplot(detected_mirs_in_expr_df[1:10], aes(x=factor(miRNA2, levels=rev(detected_mirs_in_expr_df$miRNA2)), y=experiments)) + 
  geom_col(fill="#7570b3", width = 0.75) +
  scale_y_continuous(expand=expansion(mult=c(0,0)), labels=scales::label_percent()) +
  coord_flip() +
  ylab("Detected in experiments") +
  xlab("") +
  theme_cowplot(11)


num_mirs_detected_at_least1 = nrow(mir_detected_in_exps[experiments > 0])

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
                queries=c(list(
                  upset_query(set="SBN", fill=unname(protocol_colors["SBN_CL"])),  # SBN_CL
                  upset_query(set="SBN_CL", fill=unname(protocol_colors["SB"])),         # SB
                  upset_query(set="SB", fill=protocol_colors["4N_CL"]),  # 4N_CL
                  upset_query(set="SB_4N", fill=protocol_colors["SBN"]),             # 
                  upset_query(set="4N_CL", fill=protocol_colors["4N"]), # 4N
                  upset_query(set="4N", fill=unname(protocol_colors["SB_4N"])),   # 
                  upset_query(set="4N_C3", fill=unname(protocol_colors["CL"])),    #CL
                  upset_query(set="CL", fill=protocol_colors["4N_C3"])     # 4N_C3
                ), qlist_intersect),
                matrix=intersection_matrix(
                  geom=geom_point(
                    size=2
                  )
                ),
                themes=upset_default_themes(text=element_text(size=9)),
                name="",
                sort_sets="ascending",
                n_intersections=20,
                height_ratio = 1.5,
                width_ratio = 0.3
)

save_plot("results/figures/supp_fig10g.pdf", upset_p, base_width = 100, base_height = 70, unit="mm")

supp_fig10abc = plot_grid(plot_grid(plot_grid(sample_stats_per_group_with_errorbars + theme(legend.position = "none"),
                                       sample_stats_per_group_legend, ncol=1, rel_heights = c(5, 1)), labels="a"),
                   plot_grid(
                     plot_grid(adapter_dimers_stage12_boxp, detected_mirs_p + theme(legend.position="none"), labels=c("b", "c")),
                               detected_mirs_legend, ncol=1, rel_heights = c(5,1)),
                   nrow=1, rel_widths = c(1,2))

supp_fig10dtog = plot_grid(
  seu_umap_plot,
  distance_plot,
  detected_mirs_in_exp,
  NULL,
  nrow=2,
  labels=c("d", "e", "f", "g")
)

supp_fig10atog = plot_grid(supp_fig10abc, supp_fig10dtog, ncol=1, rel_heights = c(1,2))
save_plot("results/figures/supp_fig10.pdf", supp_fig10atog, base_height = 270, base_width = 210, unit="mm")

setnames(sample_stats_mean_plot_df, c("Group", "value"), c("Protocol", "mean"))
fwrite(sample_stats_mean_plot_df[, c("variable", "Protocol", "mean", "sd"), with=F], "results/figures/supp_fig10a.txt", sep='\t')

fwrite(adapter_dimers_stage12_box, "results/figures/supp_fig10b.txt", sep='\t')

setnames(detected_mirs, "Group", "Protocol")
fwrite(detected_mirs, "results/figures/supp_fig10c.txt", sep='\t')

setnames(seu_umap, c("rn", "UMAP_1", "UMAP_2"), c("Sample", "UMAP 1", "UMAP 2"))
fwrite(seu_umap, "results/figures/supp_fig10d.txt", sep='\t')

setnames(dist_same_diff_protocols, c("value", "L1", "Protocol"), c("distance", "Protocol", "Comparison"))
fwrite(dist_same_diff_protocols[, c("distance", "Protocol", "Comparison"), with=F], "results/figures/supp_fig10e.txt", sep='\t')

fwrite(detected_mirs_in_expr_df[1:10], "results/figures/supp_fig10f.txt", sep='\t')

fwrite(UpSetR::fromList(expressed_mirs_per_protocol), "results/figures/supp_fig10g.txt", sep='\t')
