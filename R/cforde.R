library("ranger")
library("data.table")
library("doParallel")
library("foreach")
library("stringr")
library("truncnorm")

source("unoverlap_hyperrectangles.R") # routine for unoverlapping hyperrectangles
source("forge.R") # slightly modified FORGE function that can handle outputs from both FORDE and cFORDE
source("forde.R") # handles missing values
source("utils.R")
source("adversarial_rf.R")


### conditional FORDE
cforde <- function(params_uncond, cond, row_mode = c("separate", "or"), stepsize = 200, parallel = T) {
  
  row_mode <- match.arg(row_mode)
  
  meta <- params_uncond$meta
  family <- meta[family != "multinom", unique(family)]
  forest <- params_uncond$forest
  cat <- params_uncond$cat
  cnt <- params_uncond$cnt
  cnt_cols <-meta[family != "multinom", variable]
  cat_cols <-meta[family == "multinom", variable]
  
  # format c, calculate DNF and output disjoint hyperrectangles
  volumes <- cond2volumes(cond,params_uncond, row_mode)
  
  if(nrow(volumes) == 0){
    return(NULL)
  }

  # store resulting number of disjoint hyperrectangles
  
  if(row_mode == "or") {
    nvols <- nvols_conditioned <- volumes[,max(c_idx)]
  } else {
    nvols <- nrow(cond)
    nvols_conditioned <- volumes[,uniqueN(c_idx)]
  }
  
  vols_conditioned <- volumes[,unique(c_idx)]
  
  # store conditions for cat and cnt separately
  cat_conds <- volumes[variable %in% cat_cols,c("c_idx","variable","val")][,variable := factor(variable)]
  cnt_conds <- volumes[variable %in% cnt_cols,c("c_idx","variable","min", "max","val")][,`:=` (variable = factor(variable),
                                                                                                   val = as.numeric(val))]

  # subset cat according to cat_conds ##############
  
  if (nrow(cat_conds) != 0) {
    cat_relevant <- cat[,.(.(f_idx)), by=.(variable,val)]
    setkey(cat_relevant, variable, val)
    setkey(cat_conds, variable, val)
    step_no <- ceiling(nvols_conditioned/stepsize)
    cat_relevant <- cat_conds[cat_relevant, on = .(variable, val),nomatch = NULL]
    setkey(cat_relevant,c_idx)
    
    if (parallel) {
      relevant_leaves_changed_cat <- foreach(step = 1:step_no, .combine = "rbind") %dopar% {
        index_start <- vols_conditioned[(step - 1)*stepsize + 1]
        index_end <- vols_conditioned[min(step * stepsize, nvols_conditioned)]
        cat_relevant[.(index_start:index_end), Reduce(intersect,V1),by = c_idx][,.(c_idx, f_idx = V1)]
      }
    } else {
      relevant_leaves_changed_cat <- foreach(step = 1:step_no, .combine = "rbind") %do% {
        index_start <- vols_conditioned[(step - 1)*stepsize + 1]
        index_end <- vols_conditioned[min(step * stepsize, nvols_conditioned)]
        cat_relevant[.(index_start:index_end), Reduce(intersect,V1),by = c_idx][,.(c_idx, f_idx = V1)]
      }
    }
    setorder(relevant_leaves_changed_cat)
    volumes_unchanged_cat <- (1:nvols)[!(1:nvols %in% cat_conds[,c_idx])]
    relevant_leaves_unchanged_cat <- data.table(c_idx = rep(volumes_unchanged_cat, each = nrow(forest) ), f_idx = rep(forest[,f_idx],length(volumes_unchanged_cat)))
    relevant_leaves_cat <- rbind(relevant_leaves_changed_cat, relevant_leaves_unchanged_cat)
    relevant_leaves_cat_list <- relevant_leaves_cat[,.(f_idx = .(f_idx)),by=c_idx]
  } else {
    relevant_leaves_cat <- data.table(c_idx = integer(), f_idx = integer())
  }
  
  
  if (nrow(cnt_conds) != 0) {
    cnt_conds_compact <- copy(cnt_conds)
    cnt_conds_compact[!is.na(val), `:=`(min = val, max = val)][,val := NULL]
    
    cnt_relevant <- cnt[,.(min = .(min), max = .(max)),by = .(variable)]
    cnt_relevant <- cnt_conds_compact[cnt_relevant,on = .(variable)]
    
    if (nrow(cat_conds) != 0) {
      cnt_relevant <- cnt_relevant[relevant_leaves_cat_list, on = .(c_idx)]
    } else {
      cnt_relevant[,f_idx := NA]
    }
    
    setkey(cnt_relevant,c_idx)
    step_no <- ceiling(nvols_conditioned/stepsize)
    
    if(parallel) {
      relevant_leaves_changed_cnt <- foreach(step = 1:step_no, .combine = "rbind") %dopar% {
        index_start <- vols_conditioned[(step - 1)*stepsize + 1]
        index_end <- vols_conditioned[min(step * stepsize, nvols_conditioned)]
        cnt_relevant[.(index_start:index_end), .(
          c_idx,
          variable,
          f_idx = Map(\(f_idx, min, max, i.min, i.max) {
            if (class(f_idx) != "logical") {
              rel_cat_min <- i.min[f_idx]
              rel_cat_max <- i.max[f_idx]
              rel_min <- f_idx[which(max > rel_cat_min)]
              rel_max <- f_idx[which(min <= rel_cat_max)]
            } else {
              rel_cat_min <- i.min
              rel_cat_max <- i.max
              rel_min <- which(max > rel_cat_min)
              rel_max <- which(min <= rel_cat_max)
            }
            intersect(rel_min,rel_max) 
            }, f_idx = f_idx, min = min, max = max, i.min = i.min, i.max = i.max))
        ][, Reduce(intersect,f_idx),by = c_idx][,.(c_idx, f_idx = V1)]
      }
    } else {
      relevant_leaves_changed_cnt <- foreach(step = 1:step_no, .combine = "rbind") %do% {
        index_start <- vols_conditioned[(step - 1)*stepsize + 1]
        index_end <- vols_conditioned[min(step * stepsize, nvols_conditioned)]
        cnt_relevant[.(index_start:index_end), .(
          c_idx,
          variable,
          f_idx = Map(\(f_idx, min, max, i.min, i.max) {
            if (class(f_idx) != "logical") {
              rel_cat_min <- i.min[f_idx]
              rel_cat_max <- i.max[f_idx]
              rel_min <- f_idx[which(max > rel_cat_min)]
              rel_max <- f_idx[which(min <= rel_cat_max)]
            } else {
              rel_cat_min <- i.min
              rel_cat_max <- i.max
              rel_min <- which(max > rel_cat_min)
              rel_max <- which(min <= rel_cat_max)
            }
            intersect(rel_min,rel_max) 
          }, f_idx = f_idx, min = min, max = max, i.min = i.min, i.max = i.max))
        ][, Reduce(intersect,f_idx),by = c_idx][,.(c_idx, f_idx = V1)]
      }
    }

    volumes_unchanged_cnt <- (1:nvols)[!(1:nvols %in% cnt_conds[,c_idx])]
    relevant_leaves_unchanged_cnt <- data.table(c_idx = rep(volumes_unchanged_cnt, each = nrow(forest) ), f_idx = rep(forest[,f_idx],length(volumes_unchanged_cnt)))
    relevant_leaves_cnt <- rbind(relevant_leaves_changed_cnt, relevant_leaves_unchanged_cnt)
  
    cnt_new <- merge(merge(relevant_leaves_cnt, cnt_conds, by = "c_idx", allow.cartesian = T, sort = F), cnt, by = c("f_idx", "variable"), all.x = T, allow.cartesian = T, sort = F)
    cnt_new[!is.na(val),`:=` (min = min.y,
                              max = max.y)]
    cnt_new[is.na(val),`:=` (min = pmax(min.x, min.y, na.rm = T),
                             max = pmin(max.x, max.y, na.rm = T))]
    cnt_new[,cvg_factor := NA_real_]
    if (family == "truncnorm") {
      cnt_new[!is.na(val), cvg_factor := dtruncnorm(val, a=min.y, b=max.y, mean=mu, sd=sigma)*(val != min.y)]
      cnt_new[is.na(val) & (min == min.y) & (max == max.y), cvg_factor := 1]
      cnt_new[is.na(val) & is.na(cvg_factor), cvg_factor := ptruncnorm(max, a=min.y, b=max.y, mean=mu, sd=sigma) - ptruncnorm(min, a=min.y, max.y, mean=mu,sd=sigma)] 
    } else if (family == "unif") {
      cnt_new[!is.na(val), cvg_factor := dunif(val, min=min.y, max=max.y)*(val != min.y)]
      cnt_new[is.na(val) & (min == min.y) & (max == max.y), cvg_factor := 1]
      cnt_new[is.na(val) & is.na(cvg_factor), cvg_factor := punif(max, min=min.y, max=max.y) - punif(min, min=min.y, max.y)]  
    }
    cnt_new[,c("min.x","max.x","min.y","max.y") := NULL]
    if (nrow(cat_conds) > 0) {
      relevant_leaves <- merge(relevant_leaves_cnt, relevant_leaves_cat, by = c("c_idx", "f_idx"))[,.(c_idx, f_idx = .I, f_idx_uncond =f_idx)]
    } else {
      relevant_leaves <- relevant_leaves_cnt[,.(c_idx, f_idx = .I, f_idx_uncond =f_idx)]
    }
    cnt_new <- setcolorder(merge(relevant_leaves, cnt_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = F),c("f_idx","c_idx","variable","min","max","val","cvg_factor"))[]
  } else {
    relevant_leaves <- relevant_leaves_cat[,.(c_idx, f_idx = .I, f_idx_uncond =f_idx)]
    cnt_new <- data.table(f_idx = integer(), cvg_factor = numeric(), c_idx = integer())
  }

  if(relevant_leaves[,uniqueN(c_idx)] < nvols_conditioned) {
    if(relevant_leaves[,uniqueN(c_idx)] == 0 & row_mode == "or") {
      stop("For all entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds=FALSE in forde().")
    } else {
      warning("For some entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds=FALSE in forde().")
      vols_impossible <- vols_conditioned[!(vols_conditioned %in% relevant_leaves[,unique(c_idx)])]
      relevant_leaves <- setorder(rbind(relevant_leaves, data.table(c_idx = vols_impossible, f_idx = NA_integer_, f_idx_uncond = NA_integer_)))
    }
    }
  
  cat_new <- merge(merge(relevant_leaves, cat_conds, by = "c_idx", allow.cartesian = T), cat, by.x = c("f_idx_uncond","variable", "val"), by.y = c("f_idx","variable", "val")) 
  cat_new[,`:=` (cvg_factor = prob, prob = 1)]
  setcolorder(cat_new,c("f_idx","c_idx","f_idx_uncond","variable","val","prob","cvg_factor"))
  
  
  forest_new <- merge(relevant_leaves,forest, by.x = "f_idx_uncond", by.y = "f_idx", all.x = T, sort = F)
  setnames(forest_new,"cvg","cvg_arf")
  
  cvg_new <- rbind(cat_new[,.(f_idx, c_idx, cvg_factor)], cnt_new[, .(f_idx, c_idx, cvg_factor)])
  
  if(nrow(cvg_new) > 0) {
    cvg_new <- merge(cvg_new, forest_new[,.(f_idx, cvg_arf)], sort = F)
    cvg_new[, `:=` (cvg = {
      if (any(cvg_factor == 0)) {
        -Inf
      } else {
        sum(log(c(cvg_arf[1], cvg_factor)))
      }
    })
    , by = f_idx]
    
    cvg_new <- unique(cvg_new[, .(f_idx, c_idx, cvg)])
    
    if(row_mode == "or") {
      if(cvg_new[,all(cvg == - Inf)]) {
        warning("All leaves have zero likelihood. This is probably because evidence contains an (almost) impossible combination. For categorical data, consider setting alpha>0 in forde().")
        cvg_new[, cvg := 1/.N]
      } else {
        cvg_new[, cvg := exp(cvg - max(cvg))]
        cvg_new <- cvg_new[cvg > 0,][,cvg := cvg / sum(cvg)]
      }
    } else {
      if(any(cvg_new[,all(cvg == -Inf), by = c_idx][,V1])) {
        warning("All leaves have zero likelihood for some entered evidence rows. This is probably because evidence contains an (almost) impossible combination. For categorical data, consider setting alpha>0 in forde().")
        cvg_new[, cvg := 1/.N, by = c_idx]
      } else {
        cvg_new[, cvg := exp(cvg - max(cvg)), by = c_idx]
        cvg_new <- cvg_new[cvg > 0,][,cvg := cvg / sum(cvg), by = c_idx]
      }
    }
  }
  
  forest_new <- merge(forest_new, cvg_new[,.(f_idx, cvg)], all.x = T, by = "f_idx")
  
  if (row_mode == "or") {
    if(forest_new[,all(is.na(f_idx))]) {
      forest_new[is.na(f_idx), cvg := 1/.N]
    } else {
      forest_new[is.na(f_idx), cvg := 0]
    }
  } else {
    forest_new[is.na(cvg), cvg := 1]
  }
  if(row_mode == "separate" & (nvols != nvols_conditioned)) {
    vols_unconditioned <- (1:nvols)[!(1:nvols) %in% vols_conditioned]
    forest_new_unconditioned <- copy(forest)
    forest_new_unconditioned <- rbindlist(replicate(length(vols_unconditioned), forest, simplify = F))
    forest_new_unconditioned[, `:=` (c_idx = rep(vols_unconditioned,each = nrow(forest)), f_idx_uncond = f_idx, cvg_arf = cvg)]
    forest_new <- rbind(forest_new, forest_new_unconditioned)
  }
  
  setorder(setcolorder(forest_new,c("f_idx","c_idx","f_idx_uncond","tree","leaf","cvg_arf","cvg")), c_idx, f_idx, f_idx_uncond, tree, leaf)
  list(condition = cond, volumes = volumes, cnt = cnt_new, cat = cat_new, forest = forest_new)
}


### preprocessing (formatting, DNF and unoverlapping hyperrectangles)

cond2volumes <- function(condition, params_uncond, row_mode) {
  n_row_cond <- nrow(condition)
  meta <- params_uncond$meta
  cat <- params_uncond$cat
  cnt_cols <-meta[family != "multinom", variable]
  cat_cols <-meta[family == "multinom", variable]
  
  cond <- copy(condition)
  cond <- setDT(cond)
  if (length(cat_cols > 0)) {
    cond[,(cat_cols) := lapply(.SD,as.character),.SDcols = cat_cols]  
  }
  
  cols_check_range <- cond[,sapply(.SD, function(x) sum((str_sub(x,,1) == "(") | is.na(x))), .SDcols = cnt_cols]
  cols_check_or <- cond[,sapply(.SD, function(x) sum(str_detect(x, "\\|")))]
  if (row_mode == "or") {
    if (any(cols_check_range > 0 & cols_check_range < nrow(cond))){
      stop("Condition vector contains columns with both range and scalar entries. No valid conditional density can be calculated.")
    }
    if (length(cnt_cols) > 0 & n_row_cond > 1) {
      cond[,(cnt_cols) := lapply(.SD,function(col) replace(col, which(is.na(col)), "(-Inf,Inf)")),.SDcols = cnt_cols]
    }
    if (length(cat_cols) > 0 & n_row_cond > 1) {
      cond[,(cat_cols) := mapply(function(colname,col) {
        lvls_str <- paste(levels(as.factor(cat[cat$variable == colname]$val)),collapse="|")
        list(replace(col,which(is.na(col)),lvls_str))
      },colname = colnames(.SD),
      col = .SD),
      .SDcols = cat_cols]
    }
    cond <- apply(cond,1,str_split,"\\|")
    cond <- rbindlist(lapply(cond,expand.grid))
    cond <- unique(cond)
    names(cond) <- meta[, variable]
  } else if (row_mode == "separate") {
    if (any(cols_check_or > 0,na.rm = T)) {
      stop("Please use the option row_mode = 'or' when including logical disjunctions.")
    }
  }
  
  cond[,c_idx := .I]
  
  volumes <- rbind(melt(cond[,.SD,.SDcols = c("c_idx",cat_cols)], id.vars = "c_idx", value.name = "val",),
                   melt(cond[,.SD,.SDcols = c("c_idx",cnt_cols)], id.vars = "c_idx", value.name = "val",), fill = T
                   )[!is.na(val),]
  if(row_mode == "or") {
    volumes[,c_idx:= .GRP, by = c_idx]
  }
  volumes[(variable %in% cnt_cols) & str_detect(val,"\\("),c("val", "min", "max") := cbind(c(NA_real_,transpose(strsplit(substr(val, 2, nchar(val) - 1), split = ","))))]
  volumes[, c("min", "max") := lapply(.SD,as.numeric), .SDcols = c("min","max")]
  setcolorder(volumes, c("c_idx", "variable", "min", "max"))
  
  if (row_mode == 'or' & nrow(cond) > 1 & !all(cols_check_range == 0)) {
    volumes <- unoverlap_hyperrectangles(volumes, cat_cols)
    volumes <- volumes[!(min == -Inf & max == Inf & is.na(val))]
  }
  setorder(volumes, c_idx)
  volumes
}

