library(data.table)
library(ggplot2)
library(cowplot)
library(viridisLite)
library(gghalves)
library(ggridges)
library(pbapply)
library(extrafont)
library(uwot)

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
                  "Faridani"="#66a61e")

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

metadata = fread("data/supp_fig7_input/samples_me.csv")
metadata[, `Sample name`:=factor(`Sample name`, levels=metadata[order(Group, `Sample name`)]$`Sample name`)]
cutadapt = fread("data/supp_fig7_input/cutadapt_stats.csv")
star = fread("data/supp_fig7_input/mapping_overview.csv")
spike_in = fread("data/supp_fig7_input/map_spike_stats.csv")
setnames(spike_in, "Mapped reads", "Spike-in reads")
read_compo = fread("data/supp_fig7_input/overall_composition.detailed.raw.csv")
cutadapt_len_stats = fread("data/supp_fig7_input/overall.too_short.len_stats.csv")
adapter_dimers = colSums(cutadapt_len_stats[Length <= 3][,2:ncol(cutadapt_len_stats)])

############ overall sample stats ##############

sample_stats = merge(
  merge(
  merge(
    merge(cutadapt, metadata, by.y = "Sample name", by.x="Sample"),
    star, by="Sample"),
  spike_in, by="Sample"),
  read_compo, by="Sample")
sample_stats[, adapter_dimers_ratio:=(adapter_dimers[match(Sample, names(adapter_dimers))])/r_processed]
sample_stats[, mapped_to_spike_ratio:=`Spike-in reads`/r_processed]
sample_stats[, read_lost_in_qc:=(r_processed-r_written)/r_processed-adapter_dimers_ratio]
sample_stats[, unmapped_reads:=(r_written-mapped_total)/r_processed]
sample_stats[, unmapped_spikes:=(r_written-`Spike-in reads`)/r_processed]
sample_stats[, mapped_ratio:=(mapped_total - miRNA)/r_processed]
sample_stats[, mirna_ratio:=miRNA/r_processed]
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
                    name="") +
  scale_color_manual(values=c(read_lost_in_qc=darken("#d3d3d3", 1.5), adapter_dimers_ratio=darken("#adadad", 1.5), unmapped_reads=darken("#3B4992FF"), mapped_ratio=darken("#008B45FF"), mirna_ratio=darken("#c51b8a")),
                     labels=c(read_lost_in_qc="Lost in QC", adapter_dimers_ratio="Adapter dimers", unmapped_reads="Unmapped", mapped_ratio="Mapped", mirna_ratio="miRNAs"),
                     name="") +
  theme_cowplot(11) + theme(legend.position="right")

sample_stats_per_group_legend = get_legend(sample_stats_per_group_with_errorbars + guides(color=guide_legend(nrow=2), fill=guide_legend(nrow=2)) + theme(legend.box.margin = margin(0, 0, 0, 12)))


########### mapping to miRExplore ############

sample_stats_spike_plot_df = melt(sample_stats[, c("Sample", "read_lost_in_qc", "adapter_dimers_ratio", "unmapped_spikes", "mapped_to_spike_ratio")], id.vars="Sample")
sample_stats_spike_plot_df[, Group:=metadata$Group[match(Sample, metadata$`Sample name`)]]
sample_stats_spike_mean_plot_df = sample_stats_spike_plot_df[, list(value=mean(value), sd=sd(value)), by=c("variable", "Group")]
sample_stats_spike_mean_plot_df = sample_stats_spike_mean_plot_df[order(variable, decreasing = TRUE)]
sample_stats_spike_mean_plot_df[, ypos:=cumsum(value), by="Group"]

sample_stats_spike_per_group_with_errorbars = ggplot(sample_stats_spike_mean_plot_df, aes(x=factor(Group, levels=sample_stats_spike_mean_plot_df[variable == "mapped_to_spike_ratio"][order(value)]$Group), y=value, fill=variable)) +
  geom_col() +
  geom_errorbar(aes(ymax=ypos-sd, ymin=ypos-sd, color=variable), position="identity", width=0.5, show.legend = F) +
  geom_linerange(aes(ymin = ypos-sd, ymax = ypos, color=variable), show.legend=F) + 
  xlab("") + 
  scale_y_continuous(labels=scales::percent_format(), name="Total sequenced reads", expand = expansion(mult = c(0,0))) +
  coord_flip() +
  scale_fill_manual(values=c(read_lost_in_qc="#d3d3d3", adapter_dimers_ratio="#adadad", unmapped_spikes="#3B4992FF", mapped_to_spike_ratio="#ff7f00"),
                    labels=c(read_lost_in_qc="Lost in QC", adapter_dimers_ratio="Adapter dimers", unmapped_spikes="Unmapped", mapped_to_spike_ratio="Mapped to miRExplore"),
                    name="",
                    guide=guide_legend(nrow = 2)) +
  scale_color_manual(values=c(read_lost_in_qc=darken("#d3d3d3", 1.5), adapter_dimers_ratio=darken("#adadad", 1.5), unmapped_spikes=darken("#3B4992FF", 1.5), mapped_to_spike_ratio=darken("#ff7f00", 1.5)),
                    labels=c(read_lost_in_qc="Lost in QC", adapter_dimers_ratio="Adapter dimers", unmapped_spikes="Unmapped", mapped_to_spike_ratio="Mapped to miRExplore"),
                    name="",
                    guide=guide_legend(nrow = 2)) +
  theme_cowplot(11) + theme(legend.position="bottom")

########## detected spike in sequences ##########

spike_expr = fread("data/supp_fig7_input/counts.map_norm.csv")
detected_spikes = data.frame(detected=colSums(spike_expr[, 2:ncol(spike_expr)] > 0))
detected_spikes$Sample = rownames(detected_spikes)
detected_spikes = setDT(detected_spikes)
detected_spikes[, Group:=metadata$Group[match(Sample, metadata$`Sample name`)]]

detected_spikes_p = ggplot(detected_spikes, aes(x=factor(Group, levels=detected_spikes[, mean(detected), by="Group"][order(V1)]$Group), y=detected, color=Group, fill=Group)) + 
  geom_half_boxplot() +
  geom_half_point() +
  theme_cowplot(12) +
  scale_y_continuous(breaks=scales::pretty_breaks(5)) +
  scale_fill_manual(values=protocol_colors) + scale_color_manual(values=protocol_dark_colors) +
  coord_flip() +
  xlab("") + ylab("Detected miRExplore sequences") + theme(legend.position = "none")

######### distance matrix ###########
dmat_all = as.matrix(dist(log2(t(as.matrix(spike_expr[, 2:ncol(spike_expr)]) + 1)), method="euclidean"))

top8_protocols = c("CL", "SB_4N", "4N", "4N_C3", "SB", "SBN", "4N_CL", "SBN_CL")
dmat_top8 = dmat_all[metadata$Group[match(rownames(dmat_all), metadata$`Sample name`)] %in% top8_protocols, metadata$Group[match(colnames(dmat_all), metadata$`Sample name`)] %in% top8_protocols]

### for every protocol - intra protocol distance and inter protocol distance
rep_dist = lapply(top8_protocols, function(x) {
  samples = as.character(metadata[Group == x]$`Sample name`)
  rep_dist = dmat_top8[samples,samples]
  rep_dist[upper.tri(rep_dist)]
})
names(rep_dist) = top8_protocols

protocol_dist = as.data.table(melt(rep_dist))
protocol_dist[, Protocol:="Same"]

# average spearman cor between diff protocols
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
  scale_y_continuous(limits = c(0,NA), expand=expansion(mult=c(0,0.05)), breaks=scales::pretty_breaks(5)) +
  #guides(fill=guide_legend(direction="horizontal"), color=guide_legend(direction="horizontal")) +
  ylab("Euclidean distance") + xlab("")


############## UMAP embedding ###############

spike_expr_map_norm = fread("data/supp_fig7_input/counts.map_norm.csv")
set.seed(142)
all_umap = as.data.table(umap(transpose(spike_expr_map_norm[, 2:ncol(spike_expr_map_norm)]), n_neighbors = 9, scale=FALSE, min_dist=0.05))
all_umap[, Sample:=colnames(spike_expr_map_norm)[2:ncol(spike_expr_map_norm)]]
all_umap[, Group:=as.character(metadata$Group[match(Sample, metadata$`Sample name`)])]
all_umap[, Group:=ifelse(Group %in% top8_protocols, Group, "Other")]

top8_protocol_colors = protocol_colors
top8_protocol_colors[!names(top8_protocol_colors) %in% top8_protocols] = "lightgrey"
top8_protocol_colors = c(top8_protocol_colors, c(Other="lightgrey"))
top8_protocol_colors_dark = darken(top8_protocol_colors)

all_umap_plot = ggplot(all_umap[order(Group)], aes(x=V1, y=V2, color=Group, fill=Group, label=Sample)) + geom_point(size=0.5) +
  scale_fill_manual(values=top8_protocol_colors, name="Protocol") + scale_color_manual(values=top8_protocol_colors, name="Protocol") + guides(fill=guide_legend(ncol = 2, nrow=5)) +
  theme_cowplot(11) + theme(legend.position = c(0.025, 0.25), axis.text = element_blank(), axis.ticks = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")


############# Expression variance ################

spike_expr_df = melt(spike_expr_map_norm, id.vars = "Gene", variable.name = "Sample", value.name = "Count")
spike_expr_df[, Group:=metadata$Group[match(Sample, metadata$`Sample name`)]]

spike_expr_top8_df = spike_expr_df[Group %in% top8_protocols]
spike_expr_top8_df = spike_expr_top8_df[order(Count, decreasing = TRUE)]

spike_expr_top8_df_top100 = spike_expr_top8_df[, head(.SD, 100), by="Sample"]

spike_expr_stats = spike_expr_top8_df[, list(sd=sd(Count), mean=mean(Count), cv=sd(Count)/mean(Count)), by=c("Sample", "Group")]

expr_top100 = ggplot(spike_expr_top8_df_top100, aes(x=Count, y=Sample, fill=Group, color=Group)) +
  facet_grid(rows=vars(factor(Group, levels=spike_expr_stats[, mean(cv), by="Group"][order(V1)]$Group)), scales = "free_y", space = "free", switch = "y") +
  stat_density_ridges(jittered_points = TRUE, quantile_lines = TRUE, position="points_sina", point_size=0.1, size=0.1, alpha=0.5, point_alpha=0.15) +
  scale_fill_manual(values=protocol_colors) +
  scale_color_manual(values=protocol_dark_colors) +
  xlab("Expression [RPMM]") + ylab("") +
  theme_ridges(12) + theme(legend.position = "none", axis.text.y = element_blank(),
                           strip.text.y.left = element_text(angle = 0),
                           axis.title.x.bottom = element_text(hjust=0.5),
                           panel.spacing = unit(0, "pt"),
                           axis.line.x = element_line(),
                           axis.ticks.x = element_line(color="black"))

spike_expr_stats_boxplot = ggplot(spike_expr_stats, aes(x=factor(Group, levels=spike_expr_stats[, mean(cv), by="Group"][order(V1, decreasing = T)]$Group), y=cv, fill=Group, color=Group)) +
  geom_half_boxplot(outlier.colour = NA) +
  geom_half_point(size=0.5) +
  scale_y_continuous(name="Coefficient of variation", breaks=scales::pretty_breaks(6)) +
  coord_flip() +
  xlab("") + 
  scale_fill_manual(values=protocol_colors) +
  scale_color_manual(values=protocol_dark_colors) +
  theme_cowplot(11) + theme(legend.position = "none")


############# Expression variance on UMI ################

spike_expr_dedup = fread("data/supp_fig7_input/counts.dedup.csv")
spike_expr_dedup[, 2:ncol(spike_expr_dedup)] = as.data.table(as.matrix(spike_expr_dedup[,2:ncol(spike_expr_dedup)]) %*% diag(1/colSums(spike_expr_dedup[,2:ncol(spike_expr_dedup)])) * 1e6)

spike_expr_dedup_df = melt(spike_expr_dedup, id.vars = "Gene", variable.name = "Sample", value.name = "Count")
spike_expr_dedup_df[, Group:=metadata$Group[match(Sample, metadata$`Sample name`)]]

spike_expr_dedup_top8_df = spike_expr_dedup_df[Group %in% top8_protocols]
spike_expr_dedup_top8_df = spike_expr_dedup_top8_df[order(Count, decreasing = TRUE)]

spike_expr_dedup_top8_df_top100 = spike_expr_dedup_top8_df[, head(.SD, 100), by="Sample"]

spike_expr_dedup_stats = spike_expr_dedup_top8_df[, list(sd=sd(Count), mean=mean(Count), cv=sd(Count)/mean(Count)), by=c("Sample", "Group")]

expr_dedup_top100 = ggplot(spike_expr_dedup_top8_df_top100, aes(x=Count, y=Sample, fill=Group, color=Group)) +
  facet_grid(rows=vars(factor(Group, levels=spike_expr_dedup_stats[, mean(cv), by="Group"][order(V1)]$Group)), scales = "free_y", space = "free", switch = "y") +
  stat_density_ridges(jittered_points = TRUE, quantile_lines = TRUE, position="points_sina", point_size=0.1, size=0.1, alpha=0.5, point_alpha=0.15) +
  scale_fill_manual(values=protocol_colors) +
  scale_color_manual(values=protocol_dark_colors) +
  scale_x_continuous(limits=c(0, NA), expand = expansion(mult=0)) +
  xlab("Expression [CPM]") + ylab("") +
  theme_ridges(12) + theme(legend.position = "none", axis.text.y = element_blank(),
                           strip.text.y.left = element_text(angle = 0),
                           axis.title.x.bottom = element_text(hjust=0.5),
                           panel.spacing = unit(0, "pt"),
                           axis.line.x = element_line(),
                           axis.ticks.x = element_line(color="black"))


spike_expr_dedup_stats_boxplot = ggplot(spike_expr_dedup_stats, aes(x=factor(Group, levels=spike_expr_dedup_stats[, mean(cv), by="Group"][order(V1, decreasing = T)]$Group),
                                                                    y=cv, fill=Group, color=Group)) +
  geom_half_boxplot(outlier.colour = NA) +
  geom_half_point(size=0.5) +
  scale_y_continuous(name="Coefficient of variation", breaks=scales::pretty_breaks(6)) +
  xlab("") + 
  coord_flip() +
  scale_fill_manual(values=protocol_colors) +
  scale_color_manual(values=protocol_dark_colors) +
  theme_cowplot(11) + theme(legend.position = "none")



fig1defg = plot_grid(all_umap_plot,
                     distance_plot,
                     expr_top100,
                     spike_expr_stats_boxplot + theme(legend.position = "none"),
                     ncol=2, labels=c("d", "e", "f", "g"))

supp_fig7 = plot_grid(
  plot_grid(
    plot_grid(
      plot_grid(sample_stats_per_group_with_errorbars + theme(legend.position = "none"), sample_stats_spike_per_group_with_errorbars + theme(legend.position="none"), labels = c("a", "b")),
                                          sample_stats_per_group_legend, rel_heights = c(5,1), ncol=1),
                                plot_grid(detected_spikes_p, labels = "c"), nrow=1),
  fig1defg,
                      plot_grid(expr_dedup_top100, spike_expr_dedup_stats_boxplot, labels=c("h", "i")),
                      ncol=1,
                      rel_heights = c(1.4,2.5,1)
)

save_plot("results/figures/supp_fig7.pdf", supp_fig7, base_width = 210, base_height = 297, unit="mm")

setnames(sample_stats_mean_plot_df, c("Group", "value"), c("Protocol", "mean"))
fwrite(sample_stats_mean_plot_df[, c("variable", "Protocol", "mean", "sd"), with=F], "results/figures/supp_fig7a.txt", sep='\t')

setnames(sample_stats_spike_mean_plot_df, c("Group", "value"), c("Protocol", "mean"))
fwrite(sample_stats_spike_mean_plot_df[, c("variable", "Protocol", "mean", "sd"), with=F], "results/figures/supp_fig7b.txt", sep='\t')

setnames(detected_spikes, "Group", "Protocol")
fwrite(detected_spikes, "results/figures/supp_fig7c.txt", sep='\t')

setnames(all_umap, c("V1", "V2", "Group"), c("UMAP 1", "UMAP 2", "Protocol"))
fwrite(all_umap, "results/figures/supp_fig7d.txt", sep='\t')

setnames(dist_same_diff_protocols, c("value", "L1", "Protocol"), c("distance", "Protocol", "Comparison"))
fwrite(dist_same_diff_protocols[, c("distance", "Protocol", "Comparison"), with=F], "results/figures/supp_fig7e.txt", sep='\t')

setnames(spike_expr_top8_df_top100, "Group", "Protocol")
fwrite(spike_expr_top8_df_top100, "results/figures/supp_fig7f.txt", sep='\t')

setnames(spike_expr_stats, "Group", "Protocol")
fwrite(spike_expr_stats[, c("Sample", "Protocol", "cv"), with=F], "results/figures/supp_fig7g.txt", sep='\t')

setnames(spike_expr_dedup_top8_df_top100, "Group", "Protocol")
fwrite(spike_expr_dedup_top8_df_top100, "results/figures/supp_fig7h.txt", sep='\t')
setnames(spike_expr_dedup_stats, "Group", "Protocol")
fwrite(spike_expr_dedup_stats[, c("Sample", "Protocol", "cv"), with=F], "results/figures/supp_fig7i.txt", sep='\t')
