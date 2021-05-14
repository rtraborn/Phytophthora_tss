setwd("/home/rraborn/scratch/Phytop_tss/tsr")

library(TSRchitect)
library(TSRexploreR)
library(GenomicRanges)
library(DESeq2)
library(ChIPseeker)
library(ggseqlogo)
library(ggplot2)

load("/scratch/rraborn/Phytop_tss/tsr/PtSTRIPE_complete.RData")

#creating the annotation and assembly files
#update both paths below
Pt.annot <- "/scratch/rraborn/GoSTRIPES_Phytop/STRIPES/PtGENOME/GCF_000142945.1_ASM14294v1_genomic.gff"
Pt.assembly <- "/scratch/rraborn/GoSTRIPES_Phytop/STRIPES/PtGENOME/Phytop_genome.fa"

# Generate sample sheet
sample_sheet <- data.frame(
  sample_name=stringr::str_glue("PhytopSTRIPE_{seq_len(3)}"),
  file_1=NA, file_2=NA,
  condition=rep("Control", 3)
)


#writing the tss files to the workspace
tss.1 <- PtSTRIPE@tssCountData[[1]]
tss.2 <- PtSTRIPE@tssCountData[[2]]
tss.3 <- PtSTRIPE@tssCountData[[3]]

colnames(tss.1) <- c("seq","TSS", "strand", "score")
colnames(tss.2) <- c("seq","TSS", "strand", "score")
colnames(tss.3) <- c("seq","TSS", "strand", "score")

#making granges files from tss data frames
tss.1.gr <- makeGRangesFromDataFrame(tss.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.gr <- makeGRangesFromDataFrame(tss.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.3.gr <- makeGRangesFromDataFrame(tss.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

Pt.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr)
names(Pt.tss) <- c("control1", "control2", "control3")

#Creating the TSR explorer object
exp <- tsr_explorer(TSSs=Pt.tss, 
                    genome_annotation=Pt.annot, genome_assembly=Pt.assembly,
                    sample_sheet = sample_sheet
)

#Initial TSS processing

exp <- format_counts(exp, data_type="tss")

#Normalize TSSs
exp <- normalize_counts(exp, data_type = "tss", method = "DESeq2")

## TSS annotation
exp <- annotate_features(exp, data_type = "tss", feature_type="transcript")

##### Sequence analysis
## creating a truncated object for sequence analysis
## some intervals are too close to the edges of short scaffolds, so this was my workaround

plot_threshold_exploration(exp, samples="control1", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_Ptcontrol1.png")

plot_threshold_exploration(exp, samples="control2", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_Ptcontrol2.png")

plot_threshold_exploration(exp, samples="control3", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_Ptcontrol3.png")

## Correlation plot: replicates
plot_correlation(
  exp, data_type="tss",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="PhytopControl_correlation_matrix_tss.png")

### Genomic distribution analysis:
#### To update with new UTR-less annotation

plot_genomic_distribution(exp, data_type="tss", samples=c("control1", "control2","control3")) +
  scale_fill_viridis_d(direction=-1, name="Annotation")

ggsave(file="genomic_distribution_tss_phytopControl.png")

### promoter fraction plot
#### To update with new UTR-less annotation
plot_detected_features(exp, data_type="tss", samples=c("control1", "control2","control3")) +
  scale_fill_viridis_d(direction=-1)

ggsave(file="promoter_fraction_tss_Phytop_control.png")

### Density plot
plot_density(exp, data_type="tss", samples="control1")
ggsave(file="TSS_density_CDS_phytopControl1.png")

plot_density(exp, data_type="tss", samples="control2")
ggsave(file="TSS_density_CDS_phytopControl2.png")

plot_density(exp, data_type="tss", samples="control3")
ggsave(file="TSS_density_CDS_phytopControl3.png")

### TSS pileup heatmap
## need to update annotation/colour
plot_heatmap(
  exp, data_type="tss", samples="control1",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_control1.png")

plot_heatmap(
  exp, data_type="tss", samples="control2",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_control2.png")

plot_heatmap(
  exp, data_type="tss", samples="control3",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_control3.png")

### Sequence logo analysis- all three replicates
plot_sequence_logo(exp, samples="control1")
ggsave(file="control1_seq_logo.png")

plot_sequence_logo(exp, samples="control2")
ggsave(file="control2_seq_logo.png")

plot_sequence_logo(exp, samples="control3")
ggsave(file="control3_seq_logo.png")

### Dinucleotide frequency- all three replicates

plot_dinucleotide_frequencies(exp, samples="control1") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_control1.png")

plot_dinucleotide_frequencies(exp, samples="control2") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_control2.png")

plot_dinucleotide_frequencies(exp, samples="control3") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_control3.png")


### TSS Sequence colour map

plot_sequence_colormap(exp, samples="control1", rasterize=TRUE)
ggsave(file="sequence_colormap_control1.png")

plot_sequence_colormap(exp, samples="control2", rasterize=TRUE)
ggsave(file="sequence_colormap_control2.png")

plot_sequence_colormap(exp, samples="control3", rasterize=TRUE)
ggsave(file="sequence_colormap_control3.png")

### identify TSRs using clustering
exp <- tss_clustering(exp, threshold=3, n_samples=3, max_distance = 25)

# Associate TSSs with TSRs
exp <- associate_with_tsr(exp)

# Annotate TSRs
exp <- annotate_features(exp, data_type="tsr", upstream=250, downstream=100,
                         feature_type="transcript")

# Mark dominant TSS per TSR
exp <- mark_dominant(exp, data_type="tss")


# Calculate TSR metrics
exp <- tsr_metrics(exp)

###### Add tsr analysis from Standard Analysis documentation

## Correlation plot: replicates
plot_correlation(
  exp, data_type="tsr",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="control_correlation_matrix_tsr.png")

### Genomic distribution analysis:
#### To update with new UTR-less annotation

plot_genomic_distribution(exp, data_type="tsr", samples=c("control1", "control2","control4")) +
  scale_fill_viridis_d(direction=-1, name="Annotation")

ggsave(file="genomic_distribution_tsr_control.png")

### promoter fraction plot
#### To update with new UTR-less annotation
plot_detected_features(exp, data_type="tsr", samples=c("control1", "control2","control4")) +
  scale_fill_viridis_d(direction=-1)

ggsave(file="promoter_fraction_tsr_control.png")

### Density plot
plot_density(exp, data_type="tsr", samples="control1")
ggsave(file="tsr_density_CDS_control1.png")

plot_density(exp, data_type="tsr", samples="control2")
ggsave(file="tsr_density_CDS_control2.png")

plot_density(exp, data_type="tsr", samples="control4")
ggsave(file="tsr_density_CDS_control4.png")

### tsr pileup heatmap
## need to update annotation/colour
plot_heatmap(
  exp, data_type="tsr", samples="control1",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="tsr_pileup_control1.png")

plot_heatmap(
  exp, data_type="tsr", samples="control2",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="tsr_pileup_control2.png")

plot_heatmap(
  exp, data_type="tsr", samples="control4",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="tsr_pileup_control4.png")

### plot selected tsr metrics
#### TODO: to re-generate using custom script
plot_tsr_metric(exp, tsr_metrics=c("score", "width"), log2_transform=TRUE, samples="control1")
ggsave(file="plot_tsr_metrics_control1.png")

# Dinucleotide motifs by TSR shape

plot_sequence_logo(exp, dominant=TRUE, samples="control1",
                   data_conditions=conditionals(data_grouping=shape_class)
)
ggsave(file="dinucl_motif_plot_shape_control1.png")



plot_sequence_logo(exp, dominant=TRUE, samples="control2",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_control2.png")

plot_sequence_logo(exp, dominant=TRUE, samples="control4",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_control4.png")


### Plot selected gene track
#### TODO: Taylor check with Bob about rendering more clearly
gene_tracks(
  exp, feature_name="dp_gene18479",
  samples=c(TSS="antennae1", TSR="antennae1")
)

ggsave(file="gene_track_antennae_gene18479.png")