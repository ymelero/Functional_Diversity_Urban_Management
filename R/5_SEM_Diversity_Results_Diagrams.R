library(DiagrammeR)

# RAO SEM DIAGRAM
rao_simplified <- grViz("
digraph SEM_RAO {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  Tipo [label = 'Area', shape = box, style = 'rounded,filled',
        colorscheme = gnbu3, fillcolor = 1]

  V1   [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  Resp [label = 'RaoQ', shape = box, style = rounded]

  Tipo -> V1  [color = peachpuff, penwidth = 1.2]
  V1   -> Resp [color = peachpuff, penwidth = 1.2]

  graph [nodesep = 0.2, ranksep = 0.3]
}
")

rao_simplified

library(DiagrammeRsvg)
library(rsvg)
svg_code <- export_svg(rao_simplified)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_RAO.pdf")

# Bind with boxplot
#install.packages("pdftools")
library(magick)
library(grid)
library(gridExtra)

# Read PDFs
img1 <- image_read_pdf("Results/Figures/SEM_results_RAO.pdf", density = 300)
img2 <- image_read_pdf("Results/Figures/Working/Fig_Predictions_RaoQ.pdf", density = 300)

# Convert to rasterGrob 
rao_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
rao_g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)

# Combine
grid.arrange(rao_g1, rao_g2, ncol = 2, widths = c(0.45, 3))
grid.text("(a)", x = 0.01, y = 0.90, gp = gpar(fontsize = 12))
#grid.text("(b)", x = 0.2, y = 0.90, gp = gpar(fontsize = 12))


# SHANNON SEM DIAGRAM
shannon_simplified <- grViz("
digraph SEM_SHANNON {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  Tipo    [label = 'Area',      shape = box,    style = 'rounded,filled',
           colorscheme = gnbu3, fillcolor = 1]

  C3node  [label = 'C3-grass',  shape = ellipse, style = filled, fillcolor = peachpuff]
  Vivnode [label = 'Viv-grass', shape = ellipse, style = filled, fillcolor = peachpuff]

  Resp    [label = 'Shannon',   shape = box, style = rounded]

  # Significant paths
  Tipo   -> Resp   [color = '#66c2a4', penwidth = 1.4]   # Positive
  Tipo   -> C3node [color = 'orange',  penwidth = 1.4]   # Negative
  C3node -> Resp   [color = 'orange',  penwidth = 1.4]   # Negative
  Vivnode-> Resp   [color = 'orange',  penwidth = 1.0]   # Negative, weaker effect

  graph [nodesep = 0.25, ranksep = 0.35]
}
")

shannon_simplified

# Read PDFs
img1 <- image_read_pdf("Results/Figures/Working/SEM_results_Shannon2.pdf", density = 900)
img2 <- image_read_pdf("Results/Figures/Working/Fig_Predictions_Shannon.pdf", density = 900)

# Convert to rasterGrob 
shannon_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
shannon_g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)

# Combine
grid.arrange(rao_g1, rao_g2, ncol = 2, widths = c(1.1, 3))
grid.text("(b)", x = 0.01, y = 0.90, gp = gpar(fontsize = 12))
#grid.text("(b)", x = 0.2, y = 0.90, g#grid.text("(b)", x = 0.2, y = 0.90, g#grid.text("(b)", x = 0.2, y = 0.90, gp = gpar(fontsize = 12))
