library(GGally)
library(MASS)
library(ggplot2)
library(viridis)

# https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
    
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
 
  return(dens$z[ii])
}



# http://genoweb.toulouse.inra.fr/~pmartin/pgpmartin/2018/11/14/nicer-scatterplot-in-gggally/

GGscatterPlot <- function(data, mapping, ..., 
                          method = "pearson") {
  
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method)
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  
  # # PCA
  # nonNull <- x!=0 & y!=0
  # dfpc <- prcomp(~x+y, df[nonNull,])
  # df$cols <- predict(dfpc, df)[,1]
  # # Define the direction of color range based on PC1 orientation:
  # dfsum <- x+y
  # colDirection <- ifelse(dfsum[which.max(df$cols)] < 
  #                          dfsum[which.min(df$cols)],
  #                        1,
  #                        -1)
  
  
  #Get 2D density for alpha
  # dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- tryCatch({get_density(df$x, df$y)},
                         error = function(e) { rep(1, nrow(df)) }
                         )
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, color = density)) +
    geom_point(shape=16, show.legend = FALSE) +
    scale_color_viridis(option="plasma") +
    # ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)),
      mapping = ggplot2::aes(x = xlabel, 
                             y = ylabel, 
                             label = lab),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    theme_minimal()
  
  return(pp)
}