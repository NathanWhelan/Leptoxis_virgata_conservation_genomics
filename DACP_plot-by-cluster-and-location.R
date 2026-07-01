library(adegenet)
library(ape)
library(poppr)
library(na.tools)
library(devtools)
library(ggplot2)
library(cowplot)
library(ggstar)

rm(list=ls())

setwd("~/snails/L_virgata/population-genomics/virgata_m3M3n3_R80_maf025/")
data <- read.genepop(file="populations.haps.gen")

#Determine number of clusters
grp <- find.clusters(data, max.n.clust=25, n.pca = 300) ##values based on perlim run
#dapc
dapc <- dapc(data, grp$grp, n.pca = 10)

# Pull DAPC scores and metadata into a plotting data frame
# ind.coord columns 1 and 2 are the first two discriminant functions
# dapc$assign gives the cluster each individual was assigned to
plot_df <- data.frame(
  x        = dapc$ind.coord[, 1],
  y        = dapc$ind.coord[, 2],
  location = data$pop,
  cluster  = factor(dapc$assign)
)

# 11-color colorblind-friendly palette for sampling locations
# First 8 colors are Okabe-Ito, which is the gold standard for colorblind
# accessibility (safe for deuteranopia, protanopia, and tritanopia)
# Last 3 are from Paul Tol's muted palette, also colorblind safe
col_palette <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # pink
  "#000000",  # black
  "#999999",  # grey
  "#882255",  # wine
  "#44AA99",   # teal
  "white",
  "red"
)

# Integer codes for 10 distinct ggstar shapes, one per genetic cluster
# Run ggstar_shapes() to see what each number looks like
starshape_palette <- c(20, 15, 29, 25, 11,
                       14, 12, 13,  1,  3, 8, 6)

# Main DAPC scatter plot
# Color = sampling location, shape = assigned genetic cluster
# Using both aesthetics lets the reader cross-reference geography and genetics
main_plot <- ggplot(plot_df, aes(x = x, y = y, color = location,
                                 starshape = cluster, fill = location)) +
  
  # Draw individuals as stars colored by sampling location
  # Suppressing the color/fill legends here because the default ggstar legend
  # keys are stars - we want filled rectangles instead (see geom_point below)
  geom_star(size = 3, show.legend = c(fill = FALSE, color = FALSE)) +
  
  # This layer is invisible (alpha = 0) and only exists to create the
  # Sampling Location legend with filled square keys (shape 22 = filled square)
  # The override.aes in guides() below sets alpha back to 1 for the legend keys
  geom_point(aes(color = location), shape = 22, size = 4, alpha = 0,
             show.legend = TRUE) +
  
  # Color scale covers both the stars and the invisible points
  scale_color_manual(values = col_palette, name = "Sampling Location") +
  
  # Fill scale is hidden - fill colors are restored in the legend via
  # override.aes so we don't end up with a duplicate fill legend
  scale_fill_manual(values = col_palette, name = "Sampling Location",
                    guide = "none") +
  
  scale_starshape_manual(values = starshape_palette, name = "Cluster") +
  
  # Axis tick values are removed because discriminant function scores are
  # unitless - only the relative distances between points are meaningful
  labs(x = "Discriminant Function 1",
       y = "Discriminant Function 2") +
  
  # Force the Sampling Location legend to show filled colored squares
  # by overriding shape, fill, and alpha for the invisible geom_point layer
  guides(
    color = guide_legend(
      title = "Sampling Location",
      override.aes = list(shape = 22, size = 4, alpha = 1)
    ),
    starshape = guide_legend(title = "Cluster")
  ) +
  
  theme_classic() +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 12),
    axis.text         = element_blank(),        # remove unitless axis values
    axis.ticks        = element_line(linewidth = 0.75),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title        = element_text(size = 12)
  )

# DA eigenvalues screeplot
# Shows how much discrimination each DA axis captures
# Retained axes (used in the scatter plot) are black, discarded ones are grey
da_eig_df <- data.frame(
  axis       = 1:length(dapc$eig),
  eigenvalue = dapc$eig,
  retained   = c(rep("yes", dapc$n.da),
                 rep("no",  length(dapc$eig) - dapc$n.da))
)

da_scree <- ggplot(da_eig_df, aes(x = axis, y = eigenvalue, fill = retained)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("yes" = "black", "no" = "grey"), guide = "none") +
  labs(x = "", y = "", title = "DA eigenvalues") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 12),
        axis.text  = element_text(size = 12))

# PCA eigenvalues screeplot
# Shows cumulative variance explained as PCA axes are added
# Retained axes (fed into the DAPC step) are black, discarded ones are grey
# Plotted as a lollipop chart to match the default scatter.dapc() style
pca_eig_df <- data.frame(
  axis     = 1:length(dapc$pca.eig),
  cumvar   = 100 * cumsum(dapc$pca.eig) / sum(dapc$pca.eig),
  retained = c(rep("yes", dapc$n.pca),
               rep("no",  length(dapc$pca.eig) - dapc$n.pca))
)

pca_scree <- ggplot(pca_eig_df, aes(x = axis, y = cumvar, color = retained)) +
  geom_point(size = 1) +
  geom_segment(aes(xend = axis, yend = 0), linewidth = 0.5) +
  scale_color_manual(values = c("yes" = "black", "no" = "grey"), guide = "none") +
  labs(x = "PCA axis", y = "Cumulated variance (%)", title = "PCA eigenvalues") +
  ylim(0, 100) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 12),
        axis.text  = element_text(size = 12))

# Stack the two screeplots vertically with NULL spacers top and bottom
# so they appear vertically centred alongside the main plot
scree_col <- plot_grid(
  NULL,
  pca_scree,
  da_scree,
  NULL,
  nrow = 4,
  rel_heights = c(0.25, 0.25, 0.25, 0.25)
)

# Final figure - screeplots on the left, main scatter plot on the right

#pdf("DACP_virgata_m3M3n3_R80_maf025_6DA.pdf", height = 8, width = 11)
plot_grid(
  scree_col,
  main_plot,
  ncol = 2,
  rel_widths = c(0.25, 0.75)
)

#dev.off()
