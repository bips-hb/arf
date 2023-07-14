library("arf")
library("data.table")
library("doParallel")
library("foreach")
library("stringr")
library("truncnorm")

source("unoverlap_hyperrectangles.R") # routine for unoverlapping hyperrectangles
source("forge_modified.R") # slightly modified FORGE function that can handle outputs from both FORDE and cFORDE


### TODOs
# 1. Exception handling (check for input errors in cond vec)
# 2. Tidy up code
# 3. Find different name for attribute? Output params_cond of cforde has attribute 'prob_condition' which is not literally correct for lower-dim conditions (which should always be 0)
# 2. Accelerate (include row-wise calculations in one single call, change data.frames to data.tables, find faster way for unoverlapping?)
# 3. Allow more input options for cond vec? (e.g. logical NOT in input for cat and cnt range conditions, "< val","> val" instead of "(-Inf,val)", "(val,Inf)")
#    But already all possible logical statements can be formulated.
# 4. Adapt lik-function
# 5. Allow "mixed" columns for cnt conditions? (i.e. ranges and scalars)
#    Discuss! I don't think a valid cond density could always be calculated this way
#    Resulting hyperrectangles in lower-dimensional subspaces can then be completely differently oriented (-> no way to "weight" different zero-sets against each other)


### example
arf <- adversarial_rf(iris)
params_uncond <- forde(arf,iris)

# define conditional vector cond
# data.frame with colnames(cond) = colnames(data)
# resulting condition = row_1 OR ... OR row_nrow(cond) with row_i = col_i_1 AND ... AND col_i_ncol(cond)
# entry format:
# cnt: "(min,max)" for ranges
#      "value" (numeric or string) for scalar
#      expressions can be connected via logical OR "|", e.g. "5|6" or "(5,6)|(7,8)"
#      To calculate a valid cond density, ranges and scalars should not be mixed within one column (see remark in TODOs)
#      Single cond densities per row can still be calculated in "mixed" cases
#      NA equals "(-Inf,Inf)"
# cat: level as string
#      levels can be connected via logical OR "|"
#      NA equals all levels connected via logical OR, e.g. "Level1|Level2|Level3" for three levels with names "Level1","Level2","Level3"

cond <- data.frame(rbind(c("(5,7)",3,1.4,0.2,"virginica"),
                         c("(4,5)",4,2,1,"setosa")))
names(cond) <- names(iris)

# calculate cond density and sample 10 times (respecting the probabilities of entered "OR-ed" conditions)
params_cond <- cforde(params_uncond,cond)
forge_modified(params_cond,10)

# calculate cond density separately for each row in c and sample 10 times from each cond. density
# (this can be accelerated, see TODOs)
params_cond_array <- as.array(apply(cond, 1, function(x) {cforde(params_uncond,x)}))
rbindlist(apply(params_cond_array, 1, function(x) {forge_modified(x[[1]],10)}))

### conditional FORDE

cforde <- function(params_uncond,cond) {
  
  # store cond as data.frame if cond is atomic vector (will be the case if cforde is called row-wise via apply())
  if (is.atomic(cond)) {
    names_cond <- names(cond)
    cond <- transpose(data.frame(cond))
    names(cond) <- names_cond
  }
  
  meta <- params_uncond$meta
  forest <- params_uncond$forest
  cat <- params_uncond$cat
  cnt <- params_uncond$cnt
  cnt_cols <-meta[class != "factor", variable]
  cat_cols <-meta[class == "factor", variable]
  condition <- data.table(cond)
  
  # format c, calculate DNF and output disjoint hyperrectangles
  cond <- preprocess_cond(cond,params_uncond)
  
  # store resulting number of disjoint hyperrectangles
  nvols <- max(cond$volume_id)
  
  # store conditions for cat and cnt separately
  cat_conds <- cond[(cond$variable %in% cat_cols),c("volume_id","variable","val")]
  cnt_conds <- cond[(cond$variable %in% cnt_cols),c("volume_id","variable","min", "max")]
  cat_conds$variable <- factor(cat_conds$variable)
  cnt_conds$variable <- factor(cnt_conds$variable)
  rownames(cat_conds) <- NULL
  rownames(cnt_conds) <- NULL
  
  # identify for every volume which leaf matches the conditions
  relevant_leaves <- foreach(volume_id = 1:nvols, .combine = rbind) %dopar% {
    
    # store cat and cnt conditions for volume
    cat_conds_vol <- cat_conds[cat_conds$volume_id == volume_id,]
    cnt_conds_vol <- cnt_conds[cnt_conds$volume_id == volume_id,]
    
    # identify leaves that match cat conditions for volume
    if (length(cat_cols)>0) {
      cat_matches <- table(merge(cat,cat_conds_vol,by = c("variable","val"))$f_idx)
      relevant_leaves_volume_cat <- as.integer(names(cat_matches[cat_matches == nrow(cat_conds_vol)]))
      
      # only consider identified leaves for cnt
      cnt_subsetted <- cnt[cnt$f_idx %in% relevant_leaves_volume_cat,]
      cnt_subsetted$variable <- as.factor(cnt_subsetted$variable)
    } else {
      cnt_subsetted <- cnt
    }
    
    if(length(cnt_cols)>0) {
      # iterate through cnt conditions for volume
      for (cnt_cond in 1:nrow(cnt_conds_vol)) {
        
        # store variable name, min, max for condition
        cond_var <- cnt_conds_vol$variable[cnt_cond]
        cond_min <- cnt_conds_vol$min[cnt_cond]
        cond_max <- cnt_conds_vol$max[cnt_cond]
        
        # identify leaves that match condition for this variable
        cnt_matches <- cnt_subsetted[cnt_subsetted$variable == cond_var & !(cnt_subsetted$min >= cond_min & cnt_subsetted$min >= cond_max | cnt_subsetted$max <= cond_min & cnt_subsetted$max <= cond_max),]$f_idx
        
        # only consider identified leaves for next variable
        cnt_subsetted <- cnt_subsetted[cnt_subsetted$f_idx %in% cnt_matches]
      }
    } else {
      cnt_subsetted <- data.frame(f_idx = relevant_leaves_volume_cat)
    }
    
    # return identified leaves matching volume's conditions
    if (nrow(cnt_subsetted) > 0) {
      relevant_leaves_volume <- unique(cnt_subsetted$f_idx)
      data.frame(volume_id = volume_id, f_idx = relevant_leaves_volume)
    } else {
      NULL
    }
  }
  
  if(is.null(relevant_leaves)) {
    stop("Condition coverage equals 0 due to zero-probability values for categorical variables for the specified condition.")
  }
  
  # calculate new cat data.table
  if (length(cat_cols) > 0) {
    cat_new <- merge(cat_conds,merge(relevant_leaves,cat))
    cat_new <- cat_new[order(cat_new$volume_id,cat_new$f_idx,cat_new$variable),]
    rownames(cat_new)<-NULL
    names(cat_new) <- c("volume_id","variable","val","f_idx_uncond","cvg_factor")
    cat_new$f_idx <- rep(1:(nrow(cat_new)/length(cat_cols)),each = length(cat_cols))
    cat_new$prob <- 1
    cat_new <- cat_new[,c("f_idx","volume_id","f_idx_uncond","variable","val","prob","cvg_factor")]
    cat_new <- as.data.table(cat_new)
    cvg_factor_cat <- cat_new[,.(factor = prod(cvg_factor)),by = f_idx]
  } else {
    cat_new <- NULL
    cvg_factor_cat <- data.table()
    cvg_factor_cat$factor <- 1
  }
  
  # calculate new cnt data.table
  if (length(cnt_cols) > 0) {
    cnt_new <- data.table(merge(merge(relevant_leaves,cnt),cnt_conds,by=c("volume_id","variable")))
    cnt_new <- cnt_new[order(cnt_new$volume_id,cnt_new$f_idx),]
    rownames(cnt_new)<-NULL
    names(cnt_new) <- c("volume_id","variable","f_idx_uncond","min.x","max.x","mu","sigma","min.y","max.y")
    cnt_new$min <- pmax(cnt_new$min.x,cnt_new$min.y)
    cnt_new$max <- pmin(cnt_new$max.x,cnt_new$max.y)
    cnt_new[min == max, cvg_factor := 1]
    cnt_new[min != max ,cvg_factor := ptruncnorm(max, a=min.x, b=max.x, mean=mu,sd=sigma) - ptruncnorm(min, a=min.x, max.x, mean=mu,sd=sigma)]
    cnt_new[,c("min.x","max.x","min.y","max.y")] <- NULL
    cnt_new$f_idx <- rep(1:(nrow(cnt_new)/length(cnt_cols)),each = length(cnt_cols))
    cnt_new <- cnt_new[,c("f_idx","volume_id","f_idx_uncond","variable","min","max","mu","sigma","cvg_factor")]
    cvg_factor_cnt <- cnt_new[,.(factor = prod(cvg_factor)),by = f_idx]
  } else {
    cnt_new <- NULL
    cvg_factor_cnt <- data.table()
    cvg_factor_cnt$factor <- 1
  }
  
  # calculate new forest data.table
  forest_new <- merge(relevant_leaves,forest)
  forest_new <- forest_new[order(forest_new$volume_id,forest_new$f_idx),]
  rownames(forest_new)<-NULL
  names(forest_new) <- c("f_idx_uncond","volume_id","tree","leaf","cvg_arf")
  forest_new$f_idx <- as.integer(rownames(forest_new))
  cvg_new_unnormalized <- forest_new$cvg_arf * cvg_factor_cat$factor * cvg_factor_cnt$factor
  if (sum(cvg_new_unnormalized) == 0)  {
    warning("Condition coverage is numerically 0 due to highly unlikely ranges for continuous variables in the condition.")
    forest_new$cvg <- 0
  } else {
    forest_new$cvg <- cvg_new_unnormalized / sum(cvg_new_unnormalized)
  }
  forest_new <- forest_new[,c("f_idx","volume_id","f_idx_uncond","tree","leaf","cvg_arf","cvg")]
  forest_new <- as.data.table(forest_new)
  
  list(condition = condition, volumes = cond, prob_condition = sum(cvg_new_unnormalized)/max(params_uncond$forest$tree), cnt = cnt_new, cat = cat_new, forest = forest_new, meta = meta, input_class_x = params_uncond$input_class, input_class_c = class(c))
}

### preprocessing (formatting, DNF and unoverlapping hyperrectangles)

preprocess_cond <- function(cond, params_uncond) {
  meta <- params_uncond$meta
  cat <- params_uncond$cat
  cnt_cols <-meta[class != "factor", variable]
  cat_cols <-meta[class == "factor", variable]
  
  cond[cnt_cols][is.na(cond[cnt_cols])] <- "(-Inf,Inf)"
  for (cat_col in cat_cols) {
    lvls <-  levels(as.factor(cat[cat$variable == cat_col]$val))
    cond[cat_col][is.na(cond[cat_col])] <- paste(lvls,collapse="|")
  }
  
  cond <- as.data.frame(lapply(cond,str_replace_all," ",""))
  cond <- apply(cond,1,str_split,"\\|")
  cond <- lapply(cond,expand.grid)
  cond <- do.call(rbind,cond)
  cond <- cond[!duplicated(cond),]
  rownames(cond) <- NULL
  names(cond) <- meta[, variable]
  cols_check <- colSums(data.frame(lapply(cond[cnt_cols],str_detect,"\\(")))
  if (any(cols_check > 0 & cols_check < nrow(cond))){
    stop("Condition vector contains columns with both range and scalar entries. No valid conditional density can be calculated.")
  }
  scalar_cols <- names(which(cols_check == 0))
  cond[scalar_cols] <- lapply((cond[scalar_cols]),as.factor)
  factor_cols <- c(scalar_cols, cat_cols)
  lvls <- lapply(cond[factor_cols],levels)
  cond[factor_cols] <- lapply(cond[factor_cols],factorval2rng)
  cond <- format_rng(cond)
  if (!(all(cols_check == 0) | nrow(cond) == 1)) {
    cond <- unoverlap_hyperrectangles(cond)
  }
  cond$val <- NA
  if (length(cat_cols) > 0) {
    cond[cond$variable %in% cat_cols,] <- t(apply(cond[cond$variable %in% cat_cols,],1,catrng2val,lvls))
  }
  if (length(scalar_cols) > 0) {
    cond[cond$variable %in% scalar_cols,] <- t(apply(cond[cond$variable %in% scalar_cols,],1,scalarrng2val,lvls))
  }
  cond[,!names(cond) %in% c("variable","val")] <- lapply(cond[,!names(cond) %in% c("variable","val")], as.numeric)
  cond
}

catrng2val <- function(cat_row,lvls) {
  int <- mean(c(as.numeric(cat_row["min"]),as.numeric(cat_row["max"])))
  cat_row["val"] <- lvls[[cat_row["variable"]]][int]
  cat_row["min"] <- NA
  cat_row["max"] <- NA
  cat_row
}

scalarrng2val <- function(scalar_row,lvls) {
  int <- mean(c(as.numeric(scalar_row["min"]),as.numeric(scalar_row["max"])))
  scalar_row["val"] <- NA
  scalar_row["min"] <- lvls[[scalar_row["variable"]]][int]
  scalar_row["max"] <- lvls[[scalar_row["variable"]]][int]
  scalar_row
}

factorval2rng <- function(expr) {
  int <- as.numeric(expr)
  paste("(",int - 0.1,",",int + 0.1,")", sep="")
}

format_rng <- function(rng) {
  data.frame(
    volume_id = rep(1:nrow(rng), each=ncol(rng)),
    variable = factor(rep(names(rng), nrow(rng)), levels = names(rng)),
    min = as.numeric(as.vector(t(sapply(rng,str_extract,"(?<=\\().+(?=\\,)")))),
    max = as.numeric(as.vector(t(sapply(rng,str_extract,"(?<=\\,).+(?=\\))"))))
  )
}

