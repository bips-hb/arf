# data.table version (Faster)
# algorithm is not optimal (number of resulting hyperrectangles is not minimal)

unoverlap_hyperrectangles <- function(volumes) {
  
  d <- length(unique(volumes[,variable]))

  partitioned_space <- setorder(unique(data.table(variable = volumes[,variable], bound = volumes[,c(min,max)])))
  
  volumes_split <- volumes[,{
    var <- variable
    indices <- partitioned_space[(variable == var) & (bound %in% c(min,max)),,which=T]
    indices_min <- indices[1]
    indices_max <- indices[2]
    list(min = partitioned_space[indices_min:(indices_max-1),bound],max = partitioned_space[(indices_min+1):indices_max,bound],n_partials_var = indices_max - indices_min)
    },by = .(volume_id,variable)]
  
  volumes_split_cartesian <- volumes_split[,{
    lengths <- unique(.SD[,.(variable,n_partials_var)])[,n_partials_var]
    cartesian <- expand.grid(lapply(lengths,seq_len))
    subsetted <- .SD[as.vector(unlist(transpose(cartesian))) + rep(c(0,cumsum(lengths)[-d]),nrow(cartesian))]
    subsetted[, `:=` (n_partials_var = NULL, subvolume_id = rep(seq_len(.N/d),each=d), n_subvolumes = .N/d)]
    }, by = volume_id]
  
  
  volumes_split_cartesian_wide <- dcast(volumes_split_cartesian, volume_id + subvolume_id + n_subvolumes ~ variable, value.var =  c("min","max"))
  volumes_split_cartesian_wide[,N:=.N,by = eval(names(volumes_split_cartesian_wide)[4:(3+2*d)])][,N:=sum(N > 1), by = volume_id][,N:=N/n_subvolumes]
  setorder(volumes_split_cartesian_wide, N, n_subvolumes)

  volumes_unoverlapped_fragmented_wide <- unique(volumes_split_cartesian_wide,by = 4:(3+2*d))[,n_subvolumes_unoverlapped := .N, by = volume_id]
  
  volumes_unchanged <- volumes[volume_id %in% volumes_unoverlapped_fragmented_wide[n_subvolumes == n_subvolumes_unoverlapped,volume_id]]
  volumes_changed_wide <- volumes_unoverlapped_fragmented_wide[n_subvolumes != n_subvolumes_unoverlapped]
  volumes_changed <-  volumes_split_cartesian[volume_id %in% volumes_changed_wide[,volume_id] & subvolume_id %in% volumes_changed_wide[,subvolume_id],1:4]
  volumes_unoverlapped <- rbind(volumes_unchanged,volumes_changed)[,volume_id := rep(seq_len(.N/d), each = d)]
  volumes_unoverlapped[]
  
}
