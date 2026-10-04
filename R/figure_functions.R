## ===========================================================================
## figure_functions.R
##
## Plotting helpers shared by the figure scripts (ggplot2).
##
##   venn3_counts()  region counts of three gene sets
##   plot_venn3()    three-set Venn diagram with equal circles (Fig. 5d,e)
## ===========================================================================

suppressPackageStartupMessages(library(ggplot2))

FIG_INK <- "#0b0b0b"; FIG_SEC <- "#52514e"; FIG_MUT <- "#8a8a86"; FIG_GRID <- "#e4e3df"

#' @param m data.frame or matrix with three 0/1 (or logical) membership columns
#' @return named vector: A, B, C (only), AB, AC, BC (exactly two), ABC
venn3_counts <- function(m) {
  a <- m[, 1] == 1; b <- m[, 2] == 1; c3 <- m[, 3] == 1
  c(A = sum(a & !b & !c3), B = sum(!a & b & !c3), C = sum(!a & !b & c3),
    AB = sum(a & b & !c3), AC = sum(a & !b & c3), BC = sum(!a & b & c3), ABC = sum(a & b & c3))
}

#' @param m     membership columns as for venn3_counts(); column names label the sets
#' @param title plot title (character or plotmath expression)
plot_venn3 <- function(m, title) {
  n <- venn3_counts(m); sets <- colnames(m)[1:3]
  reg <- data.frame(n = n, x = c(-1.05, 1.05, 0, 0, -0.62, 0.62, 0),
                    y = c(0.62, 0.62, -1.12, 0.78, -0.25, -0.25, 0.08))
  ctr <- data.frame(set = sets, x = c(-0.6, 0.6, 0), y = c(0.35, 0.35, -0.6))
  t <- seq(0, 2 * pi, length.out = 361)
  circ <- do.call(rbind, lapply(1:3, function(i)
    data.frame(set = ctr$set[i], x = ctr$x[i] + cos(t), y = ctr$y[i] + sin(t))))
  lab <- data.frame(set = sets, x = c(-1.45, 1.45, 0), y = c(1.45, 1.45, -1.78))
  ggplot() +
    geom_polygon(data = circ, aes(x, y, group = set), fill = "grey55", alpha = 0.25,
                 colour = "grey35", linewidth = 0.4) +
    geom_text(data = reg, aes(x, y, label = format(n, big.mark = ",")), size = 3.6, colour = FIG_INK) +
    geom_text(data = lab, aes(x, y, label = set), size = 4, colour = FIG_INK) +
    coord_equal(xlim = c(-1.9, 1.9), ylim = c(-1.95, 1.7), clip = "off") +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, colour = FIG_INK),
          plot.background = element_rect(fill = "white", colour = NA))
}

#' Two-set Venn diagram.
#' @param a,b      logical membership vectors over the same universe
#' @param sets     labels of the two sets
#' @param title    plot title
#' @param subtitle text below the diagram (e.g. the overlap test)
plot_venn2 <- function(a, b, sets, title = NULL, subtitle = NULL) {
  n <- c(sum(a & !b), sum(a & b), sum(!a & b))
  t <- seq(0, 2 * pi, length.out = 361)
  circ <- rbind(data.frame(set = sets[1], x = -0.55 + cos(t), y = sin(t)),
                data.frame(set = sets[2], x =  0.55 + cos(t), y = sin(t)))
  reg <- data.frame(n = n, x = c(-1, 0, 1), y = 0)
  lab <- data.frame(set = sets, x = c(-0.9, 0.9), y = 1.3)
  ggplot() +
    geom_polygon(data = circ, aes(x, y, group = set), fill = "grey55", alpha = 0.25,
                 colour = "grey35", linewidth = 0.4) +
    geom_text(data = reg, aes(x, y, label = format(n, big.mark = ",", trim = TRUE)), size = 3.6, colour = FIG_INK) +
    geom_text(data = lab, aes(x, y, label = set), size = 3.2, colour = FIG_INK, lineheight = 0.9) +
    coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.25, 1.6), clip = "off") +
    labs(title = title, caption = subtitle) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, colour = FIG_INK),
          plot.caption = element_text(hjust = 0.5, size = 8, colour = FIG_SEC, lineheight = 1.1),
          plot.background = element_rect(fill = "white", colour = NA))
}
