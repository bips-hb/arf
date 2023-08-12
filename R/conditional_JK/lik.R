# (works with cforde, for tests)

lik <- function(x, params, arf, log = T) {
  
  cnt <- copy(params$cnt)

  if(is.null(cnt$val)) {
    cnt[,val:=NA]
    cnt[,val := as.numeric(val)]
  }
  cat <- copy(params$cat)
  forest <- copy(params$forest)[,cvg := cvg/sum(cvg)]
  
  cat_cols <-params$meta[class == "factor", variable]
  cnt_cols <-params$meta[class != "factor", variable]
  
  leaves <- as.data.table(predict(arf, x, type = 'terminalNodes')$predictions + 1)[,obs:=.I]
  n <- nrow(leaves)
  n_trees <- ncol(leaves) -1
  leaves <- melt(leaves, id="obs", variable.name = "tree", value.name= "leaf")[, tree := rep(seq_len(n_trees),each = n)]
  leaves <- merge(leaves,forest[,.(tree, leaf, f_idx)],sort=F, by=c("tree","leaf"), allow.cartesian=TRUE)[,.(obs,f_idx)]
  
  x <- setDT(x)[,obs:= .I]
  x_long <- suppressWarnings(melt(x,id.vars="obs" ,variable.name = "variable", value.name = "val"))
  
  if (length(cat_cols)>0) {
    cat_lik <- merge(merge(leaves,cat,sort=F,allow.cartesian = T, by = "f_idx"),x_long, sort=F, by=c("obs","variable","val"))[,.(obs,f_idx,prob)]
    setnames(cat_lik, "prob", "lik")
    leaves <- cat_lik[,.N,by=.(obs,f_idx)][N == length(cat_cols),.(obs,f_idx)]
  } else {
    cat_lik <- data.table()
  }
  if (length(cnt_cols)>0){
    cnt_lik <- merge(merge(leaves,cnt,sort=F,allow.cartesian = T, by = "f_idx"),x_long, sort=F, by=c("obs","variable"))
    setnames(cnt_lik,c("val.x","val.y"),c("val","obs_val"))
    cnt_lik <- cnt_lik[,lik:= is.na(val)*(obs_val!=min)*dtruncnorm(obs_val,a=min,b=max,mean=mu,sd=sigma) + (!is.na(val)&(obs_val==val))]
    cnt_lik <- cnt_lik[,.(obs,f_idx,lik)]
    # cnt_lik <- cnt_lik[,check:=any(lik==0),by=.(obs,f_idx)][check==F,-"check"]
    # 
    # if (nrow(leaves)*length(cnt_cols) > nrow(cnt_lik)) {
    #   leaves <- unique(cnt_lik[,.(obs,f_idx)])
    #   cat_lik <- merge(leaves,cat_lik,sort=F,allow.cartesian = T, by=c("obs","f_idx"))
    # }
    
  } else {
    cnt_lik <- data.table()
  }
  forest_cvg <- merge(leaves,forest, sort=F, by = "f_idx")$cvg
  
  lik <- rbind(cat_lik,cnt_lik)
  lik <- lik[,prod(lik),by = .(obs,f_idx)][,cvg := forest_cvg]
  lik <- lik[,crossprod(V1,cvg), by = obs][,V1]
  
  if (log) {
    log(lik)
  } else {
      lik
  }
}

