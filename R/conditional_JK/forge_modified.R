### slightly modified forge function that can handle outputs from both FORDE and cFORDE.
# changes marked with comments '## change'

forge_modified <- function (params, n_synth, parallel = TRUE) 
{
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- variable <- j <- f_idx <- val <- . <- NULL
  
  ## change: cFORDE-Output-DTs have more columns than FORDE-Output-DTs. Select columns that are present in FORDE-Output
  omega <- params$forest[,c('f_idx','tree','leaf','cvg')]
  ### end change
  
  draws <- data.table(f_idx = omega[, sample(f_idx, size = n_synth, 
                                             replace = TRUE, prob = cvg)])
  omega <- merge(draws, omega, sort = FALSE)[, `:=`(idx, .I)]
  synth_cnt <- synth_cat <- NULL
  if (!is.null(params$cnt)) {
    fam <- params$meta[family != "multinom", unique(family)]
    
    ## change: Do not draw from cnt distr when condition for cnt variable is scalar (val is not NA)
    if(is.null(params$cnt$val)) {
      params$cnt[,val:=NA]
      params$cnt[,val := as.numeric(val)]
    }
    psi <- merge(omega, params$cnt[,c('f_idx','variable','min','max','mu','sigma','val')], by = "f_idx", sort = FALSE, 
                 allow.cartesian = TRUE)
    
    psi[!is.na(val), dat:= val]
    if (fam == "truncnorm") {
      psi[is.na(val), dat := truncnorm::rtruncnorm(.N, a = min,
                                            b = max, mean = mu, sd = sigma)]
    }
    else if (fam == "unif") {
      psi[is.na(val), dat := stats::runif(.N, min = min, max = max)]
    }
    ### end change
    
    synth_cnt <- dcast(psi, idx ~ variable, value.var = "dat")[, 
                                                               `:=`(idx, NULL)]
  }
  
  if (!is.null(params$cat)) {
    sim_cat <- function(j) {
      psi <- merge(omega, params$cat[variable == j], by = "f_idx", 
                   sort = FALSE, allow.cartesian = TRUE)
      psi[prob == 1, `:=`(dat, val)]
      psi[prob < 1, `:=`(dat, sample(val, 1, prob = prob)), 
          by = idx]
      setnames(data.table(unique(psi[, .(idx, dat)])[, 
                                                     dat]), j)
    }
    cat_vars <- params$meta[family == "multinom", variable]
    if (isTRUE(parallel)) {
      synth_cat <- foreach(j = cat_vars, .combine = cbind) %dopar% 
        sim_cat(j)
    }
    else {
      synth_cat <- foreach(j = cat_vars, .combine = cbind) %do% 
        sim_cat(j)
    }
  }
  x_synth <- cbind(synth_cnt, synth_cat)
  setcolorder(x_synth, params$meta$variable)
  setDF(x_synth)
  idx_factor <- params$meta[, which(class == "factor")]
  idx_logical <- params$meta[, which(class == "logical")]
  idx_integer <- params$meta[, which(class == "integer")]
  if (sum(idx_factor) > 0L) {
    x_synth[, idx_factor] <- as.data.frame(lapply(x_synth[, 
                                                          idx_factor, drop = FALSE], as.factor))
  }
  if (sum(idx_logical) > 0L) {
    x_synth[, idx_logical] <- as.data.frame(lapply(x_synth[, 
                                                           idx_logical, drop = FALSE], function(x) {
                                                             x == "TRUE"
                                                           }))
  }
  if (sum(idx_integer) > 0L) {
    x_synth[, idx_integer] <- as.data.frame(lapply(x_synth[, 
                                                           idx_integer, drop = FALSE], function(x) as.integer(x)))
  }
  if ("data.table" %in% params$input_class) {
    x_synth <- as.data.table(x_synth)
  }
  else if ("matrix" %in% params$input_class) {
    x_synth <- as.matrix(x_synth)
  }
  return(x_synth)
}
