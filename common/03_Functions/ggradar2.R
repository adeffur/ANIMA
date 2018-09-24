ggradar2<-function (plot.data, font.radar = "Circular Air Light", values.radar = c("0%", 
                                                                         "50%", "100%"), axis.labels = colnames(plot.data)[-1], grid.min = 0, 
          grid.mid = 0.5, grid.max = 1, centre.y = grid.min - ((1/9) * 
                                                                 (grid.max - grid.min)), plot.extent.x.sf = 1, plot.extent.y.sf = 1.2, 
          x.centre.range = 0.02 * (grid.max - centre.y), label.centre.y = FALSE, 
          grid.line.width = 0.5, gridline.min.linetype = "longdash", 
          gridline.mid.linetype = "longdash", gridline.max.linetype = "longdash", 
          gridline.min.colour = "grey", gridline.mid.colour = "#007A87", 
          gridline.max.colour = "grey", grid.label.size = 7, gridline.label.offset = -0.1 * 
            (grid.max - centre.y), label.gridline.min = TRUE, label.gridline.mid = TRUE, 
          label.gridline.max = TRUE, axis.label.offset = 1.15, axis.label.size = 8, 
          axis.line.colour = "grey", group.line.width = 1.5, group.point.size = 6, 
          group.colours = NULL, background.circle.colour = "#D7D6D1", 
          background.circle.transparency = 0.2, plot.legend = FALSE, legend.title = "", plot.title = "", 
          legend.text.size = grid.label.size) 
{
  library(ggplot2)
  plot.data <- as.data.frame(plot.data)
  plot.data[, 1] <- as.factor(as.character(plot.data[, 1]))
  names(plot.data)[1] <- "group"
  var.names <- colnames(plot.data)[-1]
  plot.extent.x = (grid.max + abs(centre.y)) * plot.extent.x.sf
  plot.extent.y = (grid.max + abs(centre.y)) * plot.extent.y.sf
  if (length(axis.labels) != ncol(plot.data) - 1) 
    return("Error: 'axis.labels' contains the wrong number of axis labels")
  if (min(plot.data[, -1]) < centre.y) 
    return("Error: plot.data' contains value(s) < centre.y")
  if (max(plot.data[, -1]) > grid.max) 
    return("Error: 'plot.data' contains value(s) > grid.max")
  CalculateGroupPath <- function(df) {
    path <- df[, 1]
    angles = seq(from = 0, to = 2 * pi, by = (2 * pi)/(ncol(df) - 
                                                         1))
    graphData = data.frame(seg = "", x = 0, y = 0)
    graphData = graphData[-1, ]
    for (i in levels(path)) {
      pathData = subset(df, df[, 1] == i)
      for (j in c(2:ncol(df))) {
        graphData = rbind(graphData, data.frame(group = i, 
                                                x = pathData[, j] * sin(angles[j - 1]), y = pathData[, 
                                                                                                     j] * cos(angles[j - 1])))
      }
      graphData = rbind(graphData, data.frame(group = i, 
                                              x = pathData[, 2] * sin(angles[1]), y = pathData[, 
                                                                                               2] * cos(angles[1])))
    }
    colnames(graphData)[1] <- colnames(df)[1]
    graphData
  }
  CaclulateAxisPath = function(var.names, min, max) {
    n.vars <- length(var.names)
    angles <- seq(from = 0, to = 2 * pi, by = (2 * pi)/n.vars)
    min.x <- min * sin(angles)
    min.y <- min * cos(angles)
    max.x <- max * sin(angles)
    max.y <- max * cos(angles)
    axisData <- NULL
    for (i in 1:n.vars) {
      a <- c(i, min.x[i], min.y[i])
      b <- c(i, max.x[i], max.y[i])
      axisData <- rbind(axisData, a, b)
    }
    colnames(axisData) <- c("axis.no", "x", "y")
    rownames(axisData) <- seq(1:nrow(axisData))
    as.data.frame(axisData)
  }
  funcCircleCoords <- function(center = c(0, 0), r = 1, npoints = 100) {
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  plot.data.offset <- plot.data
  plot.data.offset[, 2:ncol(plot.data)] <- plot.data[, 2:ncol(plot.data)] + 
    abs(centre.y)
  group <- NULL
  group$path <- CalculateGroupPath(plot.data.offset)
  axis <- NULL
  axis$path <- CaclulateAxisPath(var.names, grid.min + abs(centre.y), 
                                 grid.max + abs(centre.y))
  axis$label <- data.frame(text = axis.labels, x = NA, y = NA)
  n.vars <- length(var.names)
  angles = seq(from = 0, to = 2 * pi, by = (2 * pi)/n.vars)
  axis$label$x <- sapply(1:n.vars, function(i, x) {
    ((grid.max + abs(centre.y)) * axis.label.offset) * sin(angles[i])
  })
  axis$label$y <- sapply(1:n.vars, function(i, x) {
    ((grid.max + abs(centre.y)) * axis.label.offset) * cos(angles[i])
  })
  gridline <- NULL
  gridline$min$path <- funcCircleCoords(c(0, 0), grid.min + 
                                          abs(centre.y), npoints = 360)
  gridline$mid$path <- funcCircleCoords(c(0, 0), grid.mid + 
                                          abs(centre.y), npoints = 360)
  gridline$max$path <- funcCircleCoords(c(0, 0), grid.max + 
                                          abs(centre.y), npoints = 360)
  gridline$min$label <- data.frame(x = gridline.label.offset, 
                                   y = grid.min + abs(centre.y), text = as.character(grid.min))
  gridline$max$label <- data.frame(x = gridline.label.offset, 
                                   y = grid.max + abs(centre.y), text = as.character(grid.max))
  gridline$mid$label <- data.frame(x = gridline.label.offset, 
                                   y = grid.mid + abs(centre.y), text = as.character(grid.mid))
  theme_clear <- theme_bw(base_size = 20) + theme(axis.text.y = element_blank(), 
                                                  axis.text.x = element_blank(), axis.ticks = element_blank(), 
                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(), legend.key = element_rect(linetype = "blank"))
  if (plot.legend == FALSE) 
    theme_clear <- theme_clear + theme(legend.position = "none")
  base <- ggplot(axis$label) + xlab(NULL) + ylab(NULL) + coord_equal() + 
    geom_text(data = subset(axis$label, axis$label$x < (-x.centre.range)), 
              aes(x = x, y = y, label = text), size = axis.label.size, 
              hjust = 1, family = font.radar) + scale_x_continuous(limits = c(-1.5 * 
                                                                                plot.extent.x, 1.5 * plot.extent.x)) + scale_y_continuous(limits = c(-plot.extent.y, 
                                                                                                                                                     plot.extent.y))
  base <- base + geom_text(data = subset(axis$label, abs(axis$label$x) <= 
                                           x.centre.range), aes(x = x, y = y, label = text), size = axis.label.size, 
                           hjust = 0.5, family = font.radar)
  base <- base + geom_text(data = subset(axis$label, axis$label$x > 
                                           x.centre.range), aes(x = x, y = y, label = text), size = axis.label.size, 
                           hjust = 0, family = font.radar)
  base <- base + theme_clear
  base <- base + geom_polygon(data = gridline$max$path, aes(x, 
                                                            y), fill = background.circle.colour, alpha = background.circle.transparency)
  base <- base + geom_path(data = axis$path, aes(x = x, y = y, 
                                                 group = axis.no), colour = axis.line.colour)
  base <- base + geom_path(data = group$path, aes(x = x, y = y, 
                                                  group = group, colour = group), size = group.line.width)
  base <- base + geom_point(data = group$path, aes(x = x, y = y, 
                                                   group = group, colour = group), size = group.point.size)
  if (plot.legend == TRUE) {
    base <- base + labs(colour = legend.title, size = legend.text.size)
  }
    
  base <- base + geom_path(data = gridline$min$path, aes(x = x, 
                                                         y = y), lty = gridline.min.linetype, colour = gridline.min.colour, 
                           size = grid.line.width)
  base <- base + geom_path(data = gridline$mid$path, aes(x = x, 
                                                         y = y), lty = gridline.mid.linetype, colour = gridline.mid.colour, 
                           size = grid.line.width)
  base <- base + geom_path(data = gridline$max$path, aes(x = x, 
                                                         y = y), lty = gridline.max.linetype, colour = gridline.max.colour, 
                           size = grid.line.width)
  if (label.gridline.min == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[1]), 
                             data = gridline$min$label, size = grid.label.size * 
                               0.8, hjust = 1, family = font.radar)
  }
  if (label.gridline.mid == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[2]), 
                             data = gridline$mid$label, size = grid.label.size * 
                               0.8, hjust = 1, family = font.radar)
  }
  if (label.gridline.max == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[3]), 
                             data = gridline$max$label, size = grid.label.size * 
                               0.8, hjust = 1, family = font.radar)
  }
  if (label.centre.y == TRUE) {
    centre.y.label <- data.frame(x = 0, y = 0, text = as.character(centre.y))
    base <- base + geom_text(aes(x = x, y = y, label = text), 
                             data = centre.y.label, size = grid.label.size, hjust = 0.5, 
                             family = font.radar)
  }
  if (!is.null(group.colours)) {
    colour_values = rep(group.colours, 100)
  }
  else {
    colour_values = rep(c("#FF5A5F", "#FFB400", "#007A87", 
                          "#8CE071", "#7B0051", "#00D1C1", "#FFAA91", "#B4A76C", 
                          "#9CA299", "#565A5C", "#00A04B", "#E54C20"), 100)
  }
  base <- base + theme(legend.key.width = unit(3, "line")) + 
    theme(text = element_text(size = 20, family = font.radar)) + 
    theme(legend.text = element_text(size = legend.text.size), 
          legend.position = "bottom") + theme(legend.key.height = unit(2, 
                                                                     "line")) + scale_colour_manual(values = colour_values) + 
    theme(text = element_text(family = font.radar)) + theme(legend.title = element_blank())
  if (plot.title != "") {
    base <- base + ggtitle(plot.title)
  }
  return(base)
}