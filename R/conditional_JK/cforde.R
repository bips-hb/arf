library("ranger")
library("data.table")
library("doParallel")
library("foreach")
library("stringr")
library("truncnorm")

source("unoverlap_hyperrectangles.R") # routine for unoverlapping hyperrectangles
source("forge.R") # slightly modified FORGE function that can handle outputs from both FORDE and cFORDE
source("../forde.R")
source("../utils.R")
source("adversarial_rf.R")

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


# Calculate cond density and sample 10 times (respecting the probabilities of entered "OR-ed" conditions)
cond <- data.frame(rbind(c("(4,6)",NA,NA,NA,NA),
                         c("(5,6)",NA,NA,NA,"versicolor")))
names(cond) <- names(iris)
params_cond <- cforde(params_uncond,cond)
x_c_synth <- forge(params_cond,10)

# For imputation: Calculate cond density separately for each row in c and sample once from each cond. density
cond <- data.frame(iris)
cond$Species <- NA
params_cond_array <- foreach(cond_i = 1:nrow(cond), .combine = rbind) %dopar% {cforde(params_uncond,cond[cond_i,])}
x_c_synth <- foreach(params_cond_i = 1:nrow(params_cond_array), .combine = rbind) %dopar% {forge(params_cond_array[params_cond_i,],1)}

### conditional FORDE

cforde <- function(params_uncond,cond) {
  
  # store cond as data.frame if cond is atomic vector (will be the case if cforde is called row-wise via apply())
  input_class_cond <- class(cond)
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
  
  # format c, calculate DNF and output disjoint hyperrectangles
  volumes <- cond2volumes(cond,params_uncond)
  
  # store resulting number of disjoint hyperrectangles
  nvols <- max(volumes$volume_id)
  
  # store conditions for cat and cnt separately
  cat_conds <- volumes[variable %in% cat_cols,c("volume_id","variable","val")]
  cnt_conds <- volumes[variable %in% cnt_cols,c("volume_id","variable","min", "max","val")]
  cat_conds[,variable := factor(variable)]
  cnt_conds[,`:=` (variable = factor(variable),
                   val = as.numeric(val))]
  
  # subset cat according to cat_conds
  if (length(cat_cols)>0) {
    cat_subsetted <- merge(cat,cat_conds,by = c("variable", "val"),allow.cartesian=T)[,N:=.N,by=.(volume_id,f_idx)]
    cat_subsetted <- setorder(merge(cat_subsetted,cat_conds[,.N,by=volume_id],by=c("volume_id","N")),cols="volume_id","f_idx","variable","val")[,N:=NULL]
    relevant_leaves_cat <- unique(cat_subsetted[,.(volume_id,f_idx)])
  } else {
    relevant_leaves_cat <- data.table(volume_id = rep(1:4,each=max(forest$f_idx)), f_idx = forest$f_idx)
  }
  
  # calculate cnt_new according to cnt_conds
  if (length(cnt_cols)>0) {
      cnt_new <- merge(merge(relevant_leaves_cat,cnt,by="f_idx",allow.cartesian=T),cnt_conds,by=c("volume_id","variable"))
      cnt_new[!is.na(val),`:=` (min = min.x,
                                max = max.x)]
      cnt_new[is.na(val),`:=` (min = pmax(min.x, min.y),
                               max = pmin(max.x, max.y))]
      cnt_new[,cvg_factor := NA]
      cnt_new[,cvg_factor := as.numeric(cvg_factor)]
      cnt_new[!is.na(val), cvg_factor := dtruncnorm(val, a=min.x, b=max.x, mean=mu, sd=sigma)*(val != min.x)]
      cnt_new[is.na(val) ,cvg_factor := ptruncnorm(max, a=min.x, b=max.x, mean=mu, sd=sigma) - ptruncnorm(min, a=min.x, max.x, mean=mu,sd=sigma)]
      relevant_leaves <- merge(cnt_new[cvg_factor != 0,.N,by=.(volume_id, f_idx)],cnt_conds[,.N,by=volume_id])[,.(volume_id,f_idx)]
      if(is.null(relevant_leaves)) {
        stop("Condition coverage equals 0 due to zero-probability values for categorical variables for the specified condition.")
      }
      cnt_new <- merge(relevant_leaves,cnt_new,by=c("volume_id","f_idx"))
      cnt_new[,c("min.x","max.x","min.y","max.y") := NULL]
      names(cnt_new) <- c("volume_id","f_idx_uncond","variable","mu","sigma","val","min","max","cvg_factor")
      cnt_new[, f_idx := rep(1:(nrow(cnt_new)/length(cnt_cols)),each = length(cnt_cols))]
      setcolorder(cnt_new,c("f_idx","volume_id","f_idx_uncond","variable","min","max","mu","sigma","val","cvg_factor"))
      cvg_factor_cnt <- cnt_new[,.(factor = prod(cvg_factor)),by = f_idx]
  }  else {
    cnt_new <- NULL
    cvg_factor_cnt <- data.table()
    cvg_factor_cnt[, factor := 1]
  }
  
  # calculate cat_new
  if (length(cat_cols) > 0) {
    cat_new <- merge(relevant_leaves,cat_subsetted,by=c("volume_id","f_idx"))
    names(cat_new) <- c("volume_id","f_idx_uncond","variable","val","cvg_factor")
    cat_new[, `:=` (f_idx = rep(1:(nrow(cat_new)/length(cat_cols)),each = length(cat_cols)),
                          prob = 1)]
    setcolorder(cat_new,c("f_idx","volume_id","f_idx_uncond","variable","val","prob","cvg_factor"))
    cvg_factor_cat <- cat_new[,.(factor = prod(cvg_factor)),by = f_idx]
  } else {
    cat_new <- NULL
    cvg_factor_cat <- data.table()
    cvg_factor_cat[,factor := 1]
  }
    
  # calculate new forest
  forest_new <- merge(relevant_leaves,forest)
  setorder(forest_new, cols = "volume_id",  "f_idx")
  names(forest_new) <- c("f_idx_uncond","volume_id","tree","leaf","cvg_arf")
  forest_new[,f_idx := as.integer(rownames(forest_new))]
  cvg_new_unnormalized <- forest_new$cvg_arf * cvg_factor_cat$factor * cvg_factor_cnt$factor
  if (sum(cvg_new_unnormalized) == 0)  {
    warning("Condition coverage is numerically 0 due to highly unlikely ranges for continuous variables in the condition.")
    forest_new$cvg <- 0
  } else {
    forest_new[,cvg := cvg_new_unnormalized / sum(cvg_new_unnormalized)]
  }
  setcolorder(forest_new,c("f_idx","volume_id","f_idx_uncond","tree","leaf","cvg_arf","cvg"))
  
  list(condition = cond, volumes = volumes, prob_condition = sum(cvg_new_unnormalized)/max(params_uncond$forest$tree), cnt = cnt_new, cat = cat_new, forest = forest_new, meta = meta, input_class_x = params_uncond$input_class, input_class_cond = input_class_cond)
}

### preprocessing (formatting, DNF and unoverlapping hyperrectangles)

cond2volumes <- function(cond, params_uncond) {
  meta <- params_uncond$meta
  cat <- params_uncond$cat
  cnt_cols <-meta[class != "factor", variable]
  cat_cols <-meta[class == "factor", variable]
  
  cond <- data.table(cond)
  cond[,(cat_cols) := lapply(.SD,as.character),.SDcols = cat_cols]
  
  if (length(cnt_cols) > 0) {
    cond[,(cnt_cols) := lapply(.SD,function(col) replace(col, which(is.na(col)), "(-Inf,Inf)")),.SDcols = cnt_cols]
  }
  if (length(cat_cols) > 0) {
      cond[,(cat_cols) := mapply(function(colname,col) {
        lvls_str <- paste(levels(as.factor(cat[cat$variable == colname]$val)),collapse="|")
        list(replace(col,which(is.na(col)),lvls_str))
        },colname = colnames(.SD),
            col = .SD),
        .SDcols = cat_cols]
  }
  
  #cond <- lapply(cond,str_replace_all," ",""))
  cond <- apply(cond,1,str_split,"\\|")
  cond <- rbindlist(lapply(cond,expand.grid))
  cond <- cond[!duplicated(cond),]
  names(cond) <- meta[, variable]
  cols_check <- colSums(as.data.table(lapply(cond[,..cnt_cols],str_detect,"\\("))) #adapt
  if (any(cols_check > 0 & cols_check < nrow(cond))){
    stop("Condition vector contains columns with both range and scalar entries. No valid conditional density can be calculated.")
  }
  scalar_cols <- names(which(cols_check == 0))
  if (length(scalar_cols) > 0) {
    cond[,(scalar_cols) := lapply(.SD,as.factor),.SDcols = scalar_cols]
  }
  factor_cols <- c(scalar_cols, cat_cols)
  if (length(factor_cols) > 0) {
    lvls <- lapply(cond[,..factor_cols],levels)
    cond[,(factor_cols) := lapply(.SD,function(x){
      int <- as.numeric(x)
      paste("(",int - 0.1,",",int + 0.1,")", sep="")
      }),
      .SDcols = factor_cols
      ]
  }
  
  volumes <- data.table(
    volume_id = rep(1:nrow(cond), each=ncol(cond)),
    variable = factor(rep(names(cond), nrow(cond)), levels = names(cond)),
    min = as.numeric(unlist(transpose(lapply(cond,str_extract,"(?<=\\().+(?=\\,)")))),
    max = as.numeric(unlist(transpose(lapply(cond,str_extract,"(?<=\\,).+(?=\\))"))))
  )
  if (!(all(cols_check == 0) | nrow(cond) == 1)) {
    volumes <- unoverlap_hyperrectangles(volumes)
  }
  if (length(factor_cols) > 0) {
    volumes[variable %in% factor_cols, `:=` (val = apply(.SD,1,function(x){
      lvls[[x["variable"]]][round(as.numeric(x["max"]))]
    }), min = NA, max= NA)]
  }
  volumes
}

