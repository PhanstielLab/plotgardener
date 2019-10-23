## Define function to plot Apa ####
plotApa <- function(apaPath = "data/output/upLoops_wt_apa/10000/gw/APA.txt", nl = 1, cols = c("white", "red"), zrange = NULL,
                    x = 0.5, y = 0.5, width=1, height=1, just=c("center", "center"), units = "npc", add=T){

  ## Error Checking ####
  if(width != height) warning("width must equal height for a square APA plot! Setting width and height to lowest value.")

  ## Load in required libraries
  library(grid)

  ## Read in dat from apaPath
  dat <- fread(apaPath)

  ## Remove brackets convert to numeric
  dat <- apply(dat, 2, function(x) gsub("\\[|\\]", "", x) %>% as.numeric())

  ## Set scale to 0; divide aggregate values by number of loops to get average loop strength
  dat <- dat - min(dat)
  dat <- dat/nl

  ## Define color function
  col_fun <- colorRampPalette(cols)

  ## Adjust scale to zrange ####
  if(is.null(zrange)){
    ## Convert to colorscale
    grps <- cut(dat, seq(min(dat), max(dat), length.out = 1000), include.lowest = T)
    cdat <- matrix(col_fun(1000)[grps], nrow = 21, ncol = 21)

    ## Create color values for legend scale
    scale_cols <- col_fun(1000)

  } else {
    ## Rescale data to zrange
    dat[dat >= zrange[2]] <- zrange[2]
    dat[dat <= zrange[1]] <- zrange[1]

    ## Convert to colorscale
    # cdat <- matrix(colourvalues::colour_values(dat, palette = colorRamp(c("white", "red"))((1:10)/10)), nrow=21, ncol=21)
    grps <- cut(dat, seq(zrange[1], zrange[2], length.out = 1000), include.lowest = T)
    cdat <- matrix(col_fun(1000)[grps], nrow = 21, ncol = 21)

    ## Create color values for legend scale
    scale_cols <- colourvalues::color_values(seq(zrange[1], zrange[2], zrange[2]*0.001), palette = colorRamp(cols)((1:10)/10))

    scale_cols <- col_fun(1000)

  }

  ## Plot apa as raster plot ####
  if(add == F) grid.newpage()

  ## Define viewports
  plotArea <- viewport(x = x, y = y, width = width, height = height, just = just, default.units = units, name = "plotArea")
  plotMatrix <- viewport(x = unit(x = 0.5, units = "snpc"), y = unit(x = 0.5, units = "snpc"),
                         width = unit(x = 1, units = "snpc"), height = unit(x = 1, units = "snpc"), name = "plotArea")

  ## Plot matrix ####
  pushViewport(plotArea)
  pushViewport(plotMatrix)
  grid.raster(cdat, interpolate = F)
  upViewport()
  upViewport()

  return(
    list(
      dat = dat,
      col_fun = col_fun,
      zrange = zrange
    )
  )
}
