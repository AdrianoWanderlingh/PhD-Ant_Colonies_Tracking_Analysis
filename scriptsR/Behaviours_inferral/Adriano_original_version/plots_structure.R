PLOT_OUTER_PDF <- TRUE
PLOT_INNER_PDF <- FALSE
############################
if (PLOT_OUTER_PDF==TRUE & PLOT_INNER_PDF==FALSE)
  {par(mfrow=c(3,5))
  pdf(...)
  plot(y ~ x, ....)
  }


for (i in 1:10)
  {############################
  if (PLOT_INNER_PDF==TRUE & PLOT_OUTER_PDF==FALSE){pdf(...)}
  for (j in 1:50)
    {
    if (PLOT_INNER_PDF==TRUE & PLOT_OUTER_PDF==FALSE)
      {
      plot(y ~ x, ....)
      }
    
    }
  if (PLOT_INNER_PDF==TRUE & PLOT_OUTER_PDF==FALSE){dev.off()}
  ############################
  }


if (PLOT_OUTER_PDF==TRUE & PLOT_INNER_PDF==FALSE==TRUE){dev.off()}
############################
