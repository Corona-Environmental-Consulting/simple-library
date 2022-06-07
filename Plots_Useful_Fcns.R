#************************************************************************************************************#
#
# Corona library
# Contains custom: color palette, boxplot
# Author: Juliette Kaplan
# Reviewers:
#
#************************************************************************************************************

########################################################################################################################
########################################################################################################################
########################################################################################################################
#Custom Corona color palette, there are 10 colors here, and one palette called "sunset"

corona_hex <- c(`navy` = "#003A68",
                `pink` = "#CC3B7C",
                `green` = "#46963D",
                `orange` = "#ff9300",
                `periwinkle` = "#608FDF",
                `red` = "#870028",
                `purple` = "#641B9A",
                `tan` = "#D99777",
                `aqua` = "#48CBC2",
                `bubblegum` = "#EE86E6")

gilbert_hex <- c( `navy blue`  = "#00274C",
                  `tan`        = "#B88B73",
                  `yellow`     = "#FFD24F",
                  `orange`     = "#E58E1A",
                  `red`        = "#8B0E04",
                  `purple`     = "#5A2149",
                  `blue`       = "#3B6E8F",
                  `green`      = "#78A22F",
                  `dark green` = "#566C11",
                  `cool gray`  = "#A1A1A5")


#'  Corona Color Hex Codes
#'  If you need to manually set colors use this function
#'
#' @examples coronahexlist$navy
#'
#' @export
#'
corona_hex_list <- list(navy = "#003A68",
                pink = "#CC3B7C",
                green = "#46963D",
                orange = "#ff9300",
                periwinkle = "#608FDF",
                red = "#B91522",
                purple = "#641B9A",
                tan = "#D99777",
                aqua = "#48CBC2",
                bubblegum = "#EE86E6")

##############################################################################################################
## code stolen from: https://stackoverflow.com/questions/61674217/custom-discrete-color-scale-in-ggplot-does-not-respect-order
colorRamp_d <- function (colors, n,
                         bias = 1,
                         space = c("rgb", "Lab"),
                         interpolate = c("linear",
                                         "spline"),
                         alpha = FALSE){

  # PRELIMINARY STEPS ----------------
  if (bias <= 0)
    stop("'bias' must be positive")
  if (!missing(space) && alpha)
    stop("'alpha' must be false if 'space' is specified")
  colors <- t(grDevices::col2rgb(colors, alpha = alpha)/255)
  space <- match.arg(space)
  interpolate <- match.arg(interpolate)

  # CUT THE COLOR VECTOR ----------------------

  if (space == "Lab")
    colors <- grDevices::convertColor(colors, from = "sRGB", to = "Lab")
  interpolate <- switch(interpolate, linear = stats::approxfun,
                        spline = stats::splinefun)

  # RESPECT ORDER IF NCLASSES<NCOLORS
  if (n < nrow(colors)) colors <- colors[1:n,]

  if ((nc <- nrow(colors)) == 1L) {
    colors <- colors[c(1L, 1L), ]
    nc <- 2L
  }
  x <- seq.int(0, 1, length.out = nc)^bias
  palette <- c(interpolate(x, colors[, 1L]), interpolate(x,colors[, 2L]),
               interpolate(x, colors[, 3L]), if (alpha) interpolate(x,colors[, 4L]))
  roundcolor <- function(rgb) pmax(pmin(rgb, 1), 0)
  if (space == "Lab")
    function(x) roundcolor(convertColor(cbind(palette[[1L]](x),
                                              palette[[2L]](x), palette[[3L]](x), if (alpha)
                                                palette[[4L]](x)), from = "Lab", to = "sRGB")) *
    255
  else function(x) roundcolor(cbind(palette[[1L]](x), palette[[2L]](x),
                                    palette[[3L]](x), if (alpha)
                                      palette[[4L]](x))) * 255
}


colorRampPalette_d <- function (colors, ...){
  # n: number of classes
  function(n) {
    if (n>1){
      ramp <- colorRamp_d(colors, n, ...)
      x <- ramp(seq.int(0, 1, length.out = n))
      if (ncol(x) == 4L)
        rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
      else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
    } else {
      unname(colors[1])
    }
  }
}

##############################################################################################################
##############################################################################################################
## Function to extract corona colors as hex codes

corona_manual_colors <- function(...) {
  colors <- c(...)

  if (is.null(colors))
    return (corona_hex)

  corona_hex[colors]
}

gilbert_manual_colors <- function(...) {
  colors <- c(...)

  if (is.null(colors))
    return (gilbert_hex)

  gilbert_hex[colors]
}

## create palettes
corona_palettes <- list(
  `sunset` = corona_manual_colors("navy", "pink", "green", "orange", "periwinkle", "red", "purple", "tan", "aqua", "bubblegum"),
  `g_primary` = gilbert_manual_colors("navy blue", "tan"),
  `g_cool`  = gilbert_manual_colors("navy blue", "blue", "green", "dark green"),
  `g_warm`   = gilbert_manual_colors("yellow", "orange", "red", "tan"),
  `g_dark` = gilbert_manual_colors("navy blue", "purple", "blue", "dark green"),
  `g_all`  = gilbert_manual_colors("blue", "tan", "dark green", "orange", "purple", "yellow", "red", "green", "navy blue", "cool gray"),
  `g_babyshark`  = gilbert_manual_colors("red", "orange", "yellow", "green", "blue", "purple", "cool gray"),
  `g_mamashark`  = gilbert_manual_colors("blue", "dark green", "purple", "cool gray"),
  `g_narwhal`  = gilbert_manual_colors("red", "tan", "orange", "yellow", "green", "dark green", "blue", "navy blue", "purple", "cool gray")
)

## function to interpolate color palettes
cor_pal <- function(palette = "sunset", ...) {
  pal <- corona_palettes[[palette]]
  colorRampPalette_d(pal, ...)
}

##############################################################################################################
#' Corona color palette (color)
#'
#' @param palette "sunset" is currently the only option with Corona colors (and is the default).
#' Using Town of Gilbert colors, you can call:
#' "g_primary", "g_cool", "g_warm", "g_dark", "g_all", "g_babyshark", "g_mamashark", "g_narwhal".
#' To add more custom color palettes, the source code for this package must be modified.
#'
#' @examples ggplot(data) + geom_point(aes(x = Date, y = Result, color = Site)) + scale_color_corona()
#'
#' @export
#'
scale_color_corona <- function(palette = "sunset", ...) {
  pal <- cor_pal(palette = palette)
  ggplot2::discrete_scale("colour", paste0("corona_", palette), palette = pal, ...)
}

#' Corona color palette (fill)
#'
#' @param palette "sunset" is currently the only option with Corona colors (and is the default).
#' Using Town of Gilbert colors, you can call:
#' "g_primary", "g_cool", "g_warm", "g_dark", "g_all", "g_babyshark", "g_mamashark", "g_narwhal".
#' To add more custom color palettes, the source code for this package must be modified.
#'
#' @examples ggplot(data) + geom_bar(aes(x = Site, y = Result, color = Parameter)) + scale_fill_corona()

#' @export
#'
scale_fill_corona <- function(palette = "sunset", ...) {
  pal <- cor_pal(palette = palette)
  ggplot2::discrete_scale("fill", paste0("corona", palette), palette = pal, ..., )
}

#end color palette code
########################################################################################################################
########################################################################################################################
########################################################################################################################
#' Custom Corona plot formatting (includes the color palette)
#'
#' This formatting sets the text size, however feel free to override when needed using (here are options, pick one at a time to change the size of):
#' axis.text, axis.title, legend.text, legend.title = element_text(size = )
#'
#'@param color This determines whether the scale_color_corona and scale_fill_corona functions are called as 
#'part of the theme.  If set to FALSE, they will not be included.  Use if you need a continuous color or fill scale.
#'
#'@examples ggplot(data) + geom_point(aes(x = Date, y = Result, color = Site)) + corona_theme()
#'
#' @export
#'
corona_theme <- function(color = TRUE) {
  list(ggplot2::theme_bw(),
       ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                      panel.grid.minor.y = ggplot2::element_blank(),
                      axis.text = ggplot2::element_text(size = 9, colour = "black"),
                      axis.title = ggplot2::element_text(size = 10, colour = "black"),
                      legend.text = ggplot2::element_text(size = 9, colour = "black"),
                      legend.title = ggplot2::element_text(size = 10, colour = "black"),
                      strip.text = ggplot2::element_text(size = 9, colour = "black"),
                      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                      legend.background = ggplot2::element_rect(fill = "transparent"),
                      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
                      legend.key = ggplot2::element_rect(fill = "transparent"),
                      strip.background = element_rect(color= "black", fill = "#003A68"),
                      strip.text.x = element_text(color = "white", face = "bold"),
                      strip.text.y = element_text(color = "white", face = "bold")),
       if(color) 
         scale_color_corona(),
       if(color)
         scale_fill_corona())
}
#scale_color_corona

#end custom plot formatting
########################################################################################################################
########################################################################################################################
########################################################################################################################
#' Custom box plot function that uses the 5th, 25th, 50th, 75th, and 95th quantiles
#'
#' To use this function, you need to call a ggplot() first and give it the data and aesthetic arguments
#' i.e. ggplot(Data, aes(x = , y = , fill = ))
#'
#' @param outlier Default is that outlier points are included.
#' The n count will be placed under the minimum outlier point if outliers are included.
#' If outliers are not included, the n count will be placed under the 25th percentile whisker.
#' @param cdist The distance n count will be from either the minimum outlier point or the 25th percentile whisker.
#' @param hbar Horizontal bars on the whiskers.
#' @param count Display count below the box and whisker. Default displays count.
#' @param diamond Diamond can be equal to "mean" or "90th" to place the point at the mean or 90th percentile.
#' If set to FALSE, no diamond point will plot.
#' @param fill Can be called to manually set one fill color for the box and diamond. Use corona_hex_list to specify color.
#' @param alpha Can be used to manually set one alpha for the boxes (defaults to 0.8). If an alpha scale is required, set to NA
#' @param pointsize Scales the points (diamond and outlier)
#' @param fontsize Can be used to change the size of the count
#' @param width Can be used to change the width of the boxes
#' @param overlap When set to TRUE, the position_identity function is called on all elements to place all boxes on top of each other.
#' 
#'
#'@examples gacbox <- ggplot(WQ.TCP, aes(x=Site.Name, y= Calc.Result)) +
#'  corona_boxplot(cdist = 70, fill = corona_hex_list$navy, alpha = .8)
#'
#' @export
#'
corona_boxplot <- function (outlier = TRUE, cdist = 0, hbar = TRUE, count = TRUE, diamond = "mean", 
                            fill = NA, alpha = .8, pointsize = 1, fontsize = 9, width = .9, overlap = FALSE) {

  # makes the box and whiskers
  boxwhisker <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }

  # makes whiskers for errorbar call
  lowwhisker <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25))
    names(r) <- c("ymin", "ymax")
    r
  }
  hiwhisker <- function(x) {
    r <- quantile(x, probs = c(0.75, 0.95))
    names(r) <- c("ymin", "ymax")
    r
  }

  # makes outliers below 5th percentile and above 95th percentile
  outliers <- function(x) {
    subset(x, x < quantile(x,0.05) | quantile(x,0.95) < x)
  }

  #make the outlier points or don't and add the count to the corresponding place
  if (outlier == TRUE){
    ## Adding Count - below the minimum outlier point
    ncount <- function(x){
      return(c(y = min(x) - cdist, label = length(x)))
    }
  } else if (outlier == FALSE){
    ## Adding Count - below the 0.05 whisker
    ncount <- function(x){
      return(c(y = as.numeric(quantile(x, 0.05)) - cdist, label = length(x)))
    }
  }
  
  # make postion function to allow for dodge vs overlapping boxes
  if(overlap){
    position = ggplot2::position_identity()
  } else {
    position = ggplot2::position_dodge(0.9)
  }

  # Returns a list of stat_summary for the full plot
  list(
    if (hbar)
      ggplot2::stat_summary(fun.data = lowwhisker, geom = "errorbar", position = position, width = width),
    if (hbar)
      ggplot2::stat_summary(fun.data = hiwhisker, geom = "errorbar", position = position, width = width),
    if (is.na(fill) & is.na(alpha))
      ggplot2::stat_summary(fun.data = boxwhisker, geom = "boxplot", position = position, 
                            width = width),
    if (!is.na(fill) & !is.na(alpha))
      ggplot2::stat_summary(fun.data = boxwhisker, geom = "boxplot", position = position, 
                            alpha = alpha, width = width, fill = fill),
    if (is.na(fill) & !is.na(alpha))
      ggplot2::stat_summary(fun.data = boxwhisker, geom = "boxplot", position = position, 
                            alpha = alpha, width = width),
    if (!is.na(fill) & is.na(alpha))
      ggplot2::stat_summary(fun.data = boxwhisker, geom = "boxplot", position = position, 
                            width = width, fill = fill),
    if (diamond == "mean" & is.na(fill))
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, colour = "black",
                            size = 2 * pointsize, position = position),
    if (diamond == "90th" & is.na(fill))
      ggplot2::stat_summary(fun = quantile, fun.args = list(probs = 0.9),  geom = "point", shape = 23, colour = "black",
                            size = 2 * pointsize, position = position),
    if (diamond == "mean" & !is.na(fill))
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, colour = "black",
                            size = 2 * pointsize, position = position, fill = fill),
    if (diamond == "90th" & !is.na(fill))
      ggplot2::stat_summary(fun = quantile, fun.args = list(probs = 0.9),  geom = "point", shape = 23, colour = "black",
                            size = 2 * pointsize, position = position, fill = fill),
    if (outlier)
      ggplot2::stat_summary(fun = outliers, geom= "point", position = position, size = 1.5 * pointsize),
    if (count)
    ggplot2::stat_summary(fun.data = ncount, geom = "text", position = position,
                          size = fontsize / ggplot2::.pt),
    corona_theme(),
    ggplot2::xlab(""),
    scale_fill_corona())
}

#end boxplot function
########################################################################################################################
########################################################################################################################
########################################################################################################################
# Boxplot legend

#' Legend for Corona custom box plot
#'
#' This function alone will create a ggplot object that can be grouped with a box and whisker plot using draw_plot from cowplot.
#'
#' @param outlier Default is that outlier points are included.
#' The n count will be placed under the minimum outlier point if outliers are included.
#' If outliers are not included, the n count will be placed under the 25th percentile whisker.
#' @param hbar Horizontal bars on the whiskers.
#' @param diamond Diamond can be equal to "mean" or "90th" to place the point at the mean or 90th percentile.
#' If set to FALSE, no diamond point will plot.
#' @param alpha Sets the transparency of the box
#' @param fill Can be called to manually set one fill color for the box and diamond. Use corona_hex_list to specify color.
#' Default is Corona blue, fill = NA will make an unfilled box.
#' @param fontsize Can be adjusted to change font size - may have to iterate when you export the plot
#'
#'@param legendwidth Can be adjusted to change the width of the box.  Increase the number to increase the box width.
#'
#' @examples
#' boxlegend <- corona_boxplot_legend(fontsize = 9)
#' fullbox <- ggdraw(boxandwhisker) + draw_plot(boxlegend, x = .65, y = .6, width = .35, height = .4)
#'
#' @return
#' @export
#'
corona_boxplot_legend <- function(outlier = TRUE, hbar = TRUE, count = TRUE, diamond = "mean", alpha = .8,
                                  fill = corona_hex_list$navy, fontsize = 10, legendwidth = 1){

  midx = 1.2
  farx = 1.5
  
  if (diamond == "mean"){
    diamondy <- function(x) mean(x)
    diamondlab <- "Mean"
    diamondx <- farx
  }else if (diamond == "90th") {
    diamondy <- function(x) unname(quantile(x, .90))
    diamondlab <- "90th percentile"
    diamondx <- midx
  }

  if (hbar){
    whiskerx <- farx
  } else {
    whiskerx <- midx
  }

  if (outlier) {
    ny <- function(x) min(x) - 30
    lowlab <- "<5th percentile"
    hilab <- ">95th percentile"
  } else {
    ny <- function(x) unname(quantile(x, .05)) - 30
    lowlab <- ""
    hilab <- ""
  }

  if (count) {
    ctlab <- "Number of values"
  } else {
    ctlab <- ""
  }
  # Create data to use in the boxplot legend:
  sample_df <- data.frame(parameter = "test",
                          values = c(223,95,169,430,143,274,115,474,68,86,68,366,375,12,454,287,305,219,430,462,7,287,
                                     162,161,333,30,115,430,70,269,309,157,228,4,147,90,70,31,65,110,7,13,75,173,126,67,
                                     96,374,491,-18,-35,-22,525,540))

  label_df <- data.frame(y = c(ny(sample_df$values), -30, unname(quantile(sample_df$values, .05)),
                               unname(quantile(sample_df$values, .25)), unname(quantile(sample_df$values, .50)),
                               unname(quantile(sample_df$values, .75)), unname(quantile(sample_df$values, .95)),
                               525, diamondy(sample_df$values)),
                         x = c(midx, midx, whiskerx, farx, farx, farx, whiskerx, midx, diamondx),
                         label = c(ctlab, lowlab, "5th percentile",
                                   "25th percentile", "50th percentile", "75th percentile", "95th percentile",
                                   hilab, diamondlab))

  ggplot2::ggplot(sample_df, ggplot2::aes(x = parameter, y=values)) +
    corona_boxplot(cdist = 30, fill = fill, alpha = alpha, 
                   count = count, diamond = diamond, hbar = hbar, outlier = outlier, fontsize = fontsize) +
    ggplot2::geom_text(data = label_df, ggplot2::aes(x = x, y = y,label = label),
                       size = fontsize/ggplot2::.pt, vjust = 0.5, hjust = 0) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank()) +
    ggplot2::coord_cartesian(xlim = c(1, 2.6 / legendwidth))
}


#' @rdname corona_boxplot_legend
#' @examples
#' fullbox <- ggdraw(boxandwhisker) + draw_plot(chad(x = .65, y = .6, width = .35, height = .4))
#'
#' @export
#'
chad <- corona_boxplot_legend

########################################################################################################################
########################################################################################################################
########################################################################################################################
#' Creates GDrive filepath
#' The function generates the correct file path for mac or PC
#'
#' @param project_name Put the project name that matches the name in the G Drive
#'
#' @examples corona_gdrive("HBWS TCP Study")
#'
#' @export
#'
corona_gdrive <- function(project_name = "") {
  
  if(Sys.info()["sysname"] == "Windows" & Sys.info()["login"] == "Kara Mihalik") {
    googledrive <- paste("G:/.shortcut-targets-by-id/0B-vUH_NHow8fZk1wd05LemdMT0E/Corona/Projects/", project_name, sep = "")
  } else if(Sys.info()["sysname"] == "Darwin" & Sys.info()["user"] == "libbymckenna") {
    googledrive <- paste("/Volumes/GoogleDrive/.shortcut-targets-by-id/0B-vUH_NHow8fZk1wd05LemdMT0E/Corona/Projects/", project_name, sep = "")
  } else if(Sys.info()["sysname"] == "Darwin") {
    googledrive <- paste("/Volumes/GoogleDrive/My Drive/Corona/Projects/", project_name, sep = "")
  } else if(Sys.info()["sysname"] == "Windows") {
    googledrive <- paste("G:/My Drive/Corona/Projects/", project_name, sep = "")
  } else {
    googledrive <- rstudioapi::selectDirectory("Select Project Folder")
  }
  
  return(googledrive)
  }

#end set directories function
########################################################################################################################
########################################################################################################################
########################################################################################################################
#' Generate png
#'
#' @param plot The plot that you are exporting (What it is called in R)
#' @param filename The name you want the plot to export as
#' @param scale The multiplier you can change the defualt height to. The default height is 4 in,
#' add the multiplier to change it to what you want
#'
#'@examples corona_png("Historical Water Quality Data", wqdataplot, scale = 1.5)
#'
#' @export
#'
corona_png <- function(plot, filename, scale = 1){

  ggplot2::ggsave(plot = plot, filename = paste(filename, ".png", sep = ""),
                  width = 6.5, height = 4*scale, units = "in", dpi = 300, bg = "transparent")

}

#end png formatting
########################################################################################################################
########################################################################################################################
########################################################################################################################

