# data.table version (Faster)
# algorithm is not optimal (number of resulting hyperrectangles is not minimal)

unoverlap_hyperrectangles <- function(volumes, cat_cols) {
  
  scalar_cols <- volumes[!is.na(val) & !(variable %in% cat_cols),unique(as.character(variable))]
  volumes[variable %in% c(scalar_cols, cat_cols), val_fac := as.numeric(as.factor(val))]
  volumes[variable %in% c(scalar_cols, cat_cols), c("min", "max") := .(val_fac - 0.1, val_fac + 0.1)]
  
  d <- uniqueN(volumes[,variable]) 

  partitioned_space <- setorder(unique(data.table(variable = volumes[,variable], bound = volumes[,c(min,max)])))
  
  volumes_split <- volumes[,{
    var <- variable
    indices <- partitioned_space[(variable == var) & (bound %in% c(min,max)),,which=T]
    indices_min <- indices[1]
    indices_max <- indices[2]
    list(min = partitioned_space[indices_min:(indices_max-1),bound],max = partitioned_space[(indices_min+1):indices_max,bound],n_partials_var = indices_max - indices_min)
    },by = .(c_idx,variable)]
  
  volumes_split_cartesian <- volumes_split[,{
    lengths <- unique(.SD[,.(variable,n_partials_var)])[,n_partials_var]
    cartesian <- expand.grid(lapply(lengths,seq_len))
    subsetted <- .SD[as.vector(unlist(transpose(cartesian))) + rep(c(0,cumsum(lengths)[-d]),nrow(cartesian))]
    subsetted[, `:=` (n_partials_var = NULL, subvolume_id = rep(seq_len(.N/d),each=d), n_subvolumes = .N/d)]
    }, by = c_idx]
  
  
  volumes_split_cartesian_wide <- dcast(volumes_split_cartesian, c_idx + subvolume_id + n_subvolumes ~ variable, value.var =  c("min","max"))
  volumes_split_cartesian_wide[,N:=.N,by = eval(names(volumes_split_cartesian_wide)[4:(3+2*d)])][,N:=sum(N > 1), by = c_idx][,N:=N/n_subvolumes]
  setorder(volumes_split_cartesian_wide, N, n_subvolumes)

  volumes_unoverlapped_fragmented_wide <- unique(volumes_split_cartesian_wide,by = 4:(3+2*d))[,n_subvolumes_unoverlapped := .N, by = c_idx]
  
  volumes_unchanged <- volumes[c_idx %in% volumes_unoverlapped_fragmented_wide[n_subvolumes == n_subvolumes_unoverlapped,c_idx]]
  volumes_changed_wide <- volumes_unoverlapped_fragmented_wide[n_subvolumes != n_subvolumes_unoverlapped]
  volumes_changed <-  volumes_split_cartesian[c_idx %in% volumes_changed_wide[,c_idx] & subvolume_id %in% volumes_changed_wide[,subvolume_id],]
  volumes_changed[,c_idx := volumes_unchanged[,max(c_idx)] + .GRP, by = .(c_idx, subvolume_id)]
  volumes_unoverlapped <- rbind(volumes_unchanged[,1:4],volumes_changed[,1:4])
  volumes_unoverlapped[variable %in% c(cat_cols, scalar_cols), c("min", "max", "val_fac") := .(NA_real_, NA_real_, round(min))]
  volumes_unoverlapped[variable %in% c(scalar_cols,cat_cols), val := merge(volumes_unoverlapped, unique(volumes[variable %in% c(scalar_cols, cat_cols),.(variable, val_fac, val)]), by = c("variable", "val_fac"))[,val]]
  setcolorder(volumes_unoverlapped, names(volumes))  
  volumes_unoverlapped[,-c("val_fac")]
}

