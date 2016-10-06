.libPaths("c:/Rlibs")
Sys.setenv(R_LIBS="c:/Rlibs")
path <- "c:\\Users\\mtsai\\Documents\\R\\R-3.2.0\\bin\\x64;c:\\Rtools\\bin;c:\\Rtools\\gcc-4.6.3\\bin;"
Sys.setenv(PATH=path)
Sys.getenv("PATH")

require(datasets)
require(utils)
require(grDevices)
require(graphics)
require(stats)
require(grid)
require(gridExtra)

library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(GGally)
library(gtools)
library(readr)
library(readxl)
library(shiny)
library(shinyBS)
library(lazyeval)
##library(rdrop2)
library(Rcpp)
library(RcppArmadillo)
library(BH)
library(mrgsolve)


auc.calc <- function(TIME,CONC)
{
  auc.total <- 0
  n <- length(TIME)
  for(i in 1:(n - 1)) {
    auc.int <- 0.5 * (TIME[i + 1] - TIME[i]) * (CONC[i + 1] + CONC[i])
    auc.total <- auc.int + auc.total
  }
  return(auc.total)
}

auc.calc <- function(TIME,CONC)
{
  n <- length(unique(TIME))
  sum(sapply(1:(n-1), FUN=function(x) {0.5*(TIME[x+1]-TIME[x])*(CONC[x+1]+CONC[x])}))
}


tmax.calc <- function(TIME,CONC) TIME[which.max(CONC)]
                                      


slide_theme <- function (...) {
  theme(strip.text = element_text(size = rel(2)),
          axis.text = element_text(size = rel(2)),
          axis.title = element_text(size = rel(2.25)),
          legend.text = element_text(size = rel(2)),
          #panel.margin = unit(2, "lines"),
          legend.title = element_text(size = rel(2.25)),
          plot.title = element_text(size = rel(3)),...)
}

decode <- function(x, search, replace, default = NULL) {
 
    # build a nested ifelse function by recursion
    decode.fun <- function(search, replace, default = NULL)
        if (length(search) == 0L) {
            function(x) if (is.null(default)) x else rep(default, length(x))
        } else {
            function(x) ifelse(x == search[1L], replace[1L],
                                                decode.fun(tail(search,  -1L),
                                                           tail(replace, -1L),
                                                           default)(x))
        }
 
    return(decode.fun(search, replace, default)(x))
}


# define the summary function
bxplt.pct <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# sample data
#d <- data.frame(x=gl(2,50), y=rnorm(1000))

# do it
#ggplot(d, aes(x, y)) + stat_summary(fun.data = bxplt.pct, geom="boxplot")

# example with outliers
# define outlier as you want    
bxplt.out <- function(x) {
  subset(x, x < quantile(x,0.05) | quantile(x,0.95) < x)
}

# do it
#ggplot(d, aes(x, y)) + 
#  stat_summary(fun.data=bxplt.pct, geom="boxplot") + 
#  stat_summary(fun.y = bxplt.out, geom="point",alpha=0.05)

read.nonmem <- function(file, n=-1) {
    ## auxiliary function to split text lines at blanks
    my.split <- function(line, numeric=FALSE) {
        pieces <- unlist(strsplit(line, split=" +"))[-1]
        if( numeric )
            return(as.numeric(pieces))
        else
            return(pieces)
    }

    cat(sprintf("Reading NONMEM data from '%s'\n", file))
    lines <- readLines(file, n) # read file as text
    cat(sprintf("- %d lines\n", length(lines)))

    idx <- substring(lines,1,1)!="T"
    cat(sprintf("- %d tables\n", sum(!idx)))
    lines <- lines[idx] # strip lines starting with T (TABLE NO ...)

    ## do we have header lines??
    if( length(grep("^ +[A-Za-z]", lines[1]))>0 ) { # yes!
        data <- sapply(lines[lines!= lines[1]], my.split, numeric=TRUE)
        header <-  my.split(lines[1])
        cat(sprintf("- file has column names (%s)\n", paste(header,collapse=", ")))
    } else {                                        # no
        data <- sapply(lines, my.split, numeric=TRUE)
        header <- sprintf("Column%02d", 1:nrow(data)) # make fake column names
        cat("- file has NO header names - creating names Column01, Column02, ...\n")
    }
    cat(sprintf("- %d columns\n", nrow(data)))

    ## transpose data and make a data.frame
    df <- data.frame(data[1,])
    for( i in 2:nrow(data))
        df <- cbind(df, data[i,])

    ## set column and row names
    rownames(df) <- NULL
    colnames(df) <- header
    cat("ok.\n")
    return(df)
}


