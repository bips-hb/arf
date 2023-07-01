### mostly plain copy from 'https://github.com/holub008/snippets/blob/master/overlapped_hyperrectangles/overlapped_hyperrectangles.Rmd' (author: Karl Holub)

## ----include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(dplyr)
library(ggplot2)
library(fuzzyjoin)

plot_rectangles <- function(volumes) {
# only allow 2 dimensional spaces for plotting
stopifnot(setequal(volumes$variable, c('x','y')))

# TODO can't remember the right tidyr encantation
volumes %>%
  group_by(volume_id) %>%
  summarize(
    ymin = max(ifelse(variable== 'y', min, -Inf)),
    ymax = max(ifelse(variable== 'y', max, -Inf)),
    xmin = max(ifelse(variable== 'x', min, -Inf)),
    xmax = max(ifelse(variable== 'x', max, -Inf)) 
  ) %>%
  mutate(volume_id = as.factor(volume_id)) %>%
  ggplot() +
  geom_rect(aes(xmin = xmin, xmax = xmax, 
                ymin = ymin, ymax = ymax,
                fill = volume_id), alpha = .5)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
build_fully_partitioned_space <- function(volumes) {
  volumes %>%
    mutate(bound = min) %>%
    select(variable, bound) %>%
    rbind(
      volumes %>%
        mutate(bound = max) %>%
        select(variable, bound),
      stringsAsFactors = FALSE
    )
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
generate_volumes_from_partitioned_space <- function(partitioned_space, id_starter = 1) {
  if (nrow(partitioned_space) == 0) {
    return(data.frame())
  }
  
  # pick an arbtirary first dimension
  dimension_of_interest <- partitioned_space$variable[1]
  dimension_bounds <- partitioned_space %>% 
    filter(variable== dimension_of_interest) %>%
    # this is a small optimization - equal bounds are redundant
    distinct() %>%
    arrange(bound)
  
  # there should always be 2 or more, since each bound corresponds to hyperrectangle edge
#  #stopifnot(nrow(dimension_bounds) > 1)
  
  # subspace meaning everything outside the dimension of interest
  partitioned_subspace <- partitioned_space %>% filter(variable!= dimension_of_interest)
  # recursively build ranges from the subspace before tacking on ranges for the dimension of interest in this stack frame
  subspace_volumes <- generate_volumes_from_partitioned_space(partitioned_subspace, id_starter = id_starter)
  
  # "expanded" by the dimension of interest, that is 
  expanded_volumes <- data.frame()
  for (bound_ix in 1:(nrow(dimension_bounds) - 1)) {
    # note that we are iterating on the sorted bounds
    lower_bound <- dimension_bounds$bound[bound_ix]
    upper_bound <- dimension_bounds$bound[bound_ix + 1]
    
    if (nrow(subspace_volumes) == 0) {
      # case this is the first dimension - there's nothing to add onto
      volume_id <- paste0(id_starter, '_', dimension_of_interest, '_', bound_ix)
      new_dimension_bounds <- list(volume_id = volume_id, 
                                   min = lower_bound, 
                                   max = upper_bound,
                                   variable= dimension_of_interest)
    }
    else {
      # case this is after the first dimension - create a new volume for each subspace volume with the new bounds added (cartesian product)
      new_dimension_bounds <- lapply(unique(subspace_volumes$volume_id), function(volume_id) {
        list(volume_id = paste0(volume_id, '_', dimension_of_interest, '_', bound_ix), # TODO this form of creating an ID could get costly in higher dimensions
             min = lower_bound, 
             max = upper_bound,
             variable= dimension_of_interest)
      }) %>% bind_rows() %>%
        rbind(subspace_volumes %>%
                mutate(volume_id = paste0(volume_id, '_', dimension_of_interest, '_', bound_ix)),
              stringsAsFactors= FALSE)
    }
    
    expanded_volumes <- rbind(expanded_volumes, new_dimension_bounds,
                              stringsAsFactors = FALSE)
  }
  
  return(expanded_volumes)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prune_uncovering_volumes <- function(new_volumes, original_volumes) {
    # we left join because not all new volumes belong to all old volumes
    # the range join prescribes that the original volumes contains the new volume
    original_to_new_volumes <- fuzzy_left_join(original_volumes, new_volumes,
                                         by = c('min' = 'min',
                                                'max' = 'max',
                                                'variable' = 'variable'),
                                         match_fun = c(`<=`, `>=`, `==`)) %>%
      # renaming some things in a reasonable way
      mutate(variable= variable.x) %>%
      select(-variable.x, -variable.y)
  
  covering_volumes <- data.frame()
  for (new_volume_id_to_check in unique(new_volumes$volume_id)) {
    volume <- new_volumes %>%
      filter(volume_id == new_volume_id_to_check)
    
    in_covering_space <- FALSE
    for (original_volume_id_to_check in unique(original_volumes$volume_id)) {
      original_volume_to_check <- original_to_new_volumes %>%
        filter(volume_id.x == original_volume_id_to_check)
      # here we make sure all dimensions are contained
      volume_dimensions_contained <- original_to_new_volumes %>%
        filter(volume_id.x == original_volume_id_to_check & 
               volume_id.y == new_volume_id_to_check) %>%
        pull(variable) %>% 
        setequal(original_volume_to_check$variable)
      
      if (volume_dimensions_contained) {
        in_covering_space <- TRUE
        break
      }
    }
    
    if (in_covering_space) {
      covering_volumes <- rbind(covering_volumes,
                             volume,
                             stringsAsFactors = FALSE)
    }
  }
  
  covering_volumes
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unoverlap_hyperrectangles <- function(volumes) {
partitioned_space <- build_fully_partitioned_space(volumes)
new_volumes <- generate_volumes_from_partitioned_space(partitioned_space)
prune_uncovering_volumes(new_volumes, volumes)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
compute_hypervolume <- function(volumes) {
  unoverlap_hyperrectangles(volumes) %>%
    group_by(volume_id) %>%
    summarize(
      hypervolume = Reduce(`*`, max - min)
    ) %>%
    summarize(
      total_hypervolume = sum(hypervolume)
    ) %>% 
    pull(total_hypervolume)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fuse_abutted_hyperrectangles <- function(volumes, original_volumes) {
  dimensionality <- n_distinct(volumes$variable)
  
  fused_volumes <- volumes
  fuses_possible <- TRUE
  fused_volume_unique_id <- 1
  colnames <- names(original_volumes)
  
  while (fuses_possible) {
    fuses_possible <- FALSE
    
    candidate_fuses <- fused_volumes %>%
      inner_join(fused_volumes, by = c('variable' = 'variable',
                                       'max' = 'min'),
                 suffix = c('.left', '.right')) %>%
      filter(volume_id.left != volume_id.right) %>%  # this should only happen if a range is of size 0
      mutate(
        max = max.right # since the left max (where the abuttment happens on the right min) must be less than the right max
      ) %>%
      distinct(variable, volume_id.left, volume_id.right, min, max) 
    
    # note this is a one to many maping, since the originals are overlapped
    current_volumes_to_original <- fused_volumes %>%
      fuzzy_inner_join(original_volumes, by = c('min' = 'min',
                                                'max' = 'max',
                                                'variable' = 'variable'),
                       match_fun = c(`>=`, `<=`, `==`)) %>%
      group_by(volume_id.x, volume_id.y) %>%
      filter(n_distinct(variable.x) == dimensionality) %>%
      summarize(
        volume_id = volume_id.x[1],
        original_volume_id = volume_id.y[1]
      ) %>%
      ungroup() %>%
      select(volume_id, original_volume_id)
    
    for (candidate_fuse_ix in seq_len(nrow(candidate_fuses))) {
      candidate_fuse <- candidate_fuses[candidate_fuse_ix, ]
      # subvolume because we ignore the dimension of the fuse
      subvolume_left <- fused_volumes %>%
        filter(volume_id == candidate_fuse$volume_id.left & variable!= candidate_fuse$variable)
      subvolume_right <- fused_volumes %>%
        filter(volume_id == candidate_fuse$volume_id.right & variable!= candidate_fuse$variable)
      
      # this case implies the volume has already been joined
      # meaning the candidate fuse may not be valid any longer - catch it next iteration
      if (nrow(subvolume_left) == 0 || nrow(subvolume_right) == 0) {
        next()
      }
      
      stopifnot(nrow(subvolume_left) == nrow(subvolume_right) && 
                  n_distinct(subvolume_left$variable) == n_distinct(subvolume_right$variable))
      
      dimension_matches <- subvolume_left %>%
        inner_join(subvolume_right, by = c('variable' = 'variable',
                                           'min' = 'min',
                                           'max' = 'max'))
      
      original_volume_counts <- current_volumes_to_original %>%
        filter(volume_id %in% c(candidate_fuse$volume_id.left, candidate_fuse$volume_id.right)) %>%
        group_by(original_volume_id) %>%
        count() %>%
        pull(n)
      
      if (nrow(dimension_matches) == dimensionality - 1 &&
          all(original_volume_counts == 2)) { #
        fuses_possible <- TRUE
        
        # add in the new volume
        fused_volume <- rbind(
          dimension_matches %>% select(min, max, variable),
          candidate_fuse %>% select(min, max, variable),
          stringsAsFactors = FALSE)
        fused_volume$volume_id <- paste0(candidate_fuse$volume_id.left, '_', 
                                         candidate_fuse$volume_id.right, '_', 
                                         as.character(fused_volume_unique_id))
        fused_volume_unique_id <- fused_volume_unique_id + 1
        fused_volumes <- rbind(fused_volumes,
                               fused_volume,
                               stringsAsFactors = FALSE)
        
        # clean up the old volumes
        fused_volumes <- fused_volumes %>%
          filter(volume_id != candidate_fuse$volume_id.left & volume_id != candidate_fuse$volume_id.right)
      }
    }
  }
  fused_volumes <- as.data.frame(fused_volumes)
  fused_volumes$volume_id <- rep(1:length(unique(fused_volumes$volume_id)),each = length(unique(fused_volumes$variable)))
  fused_volumes <- fused_volumes[colnames]
  fused_volumes <- fused_volumes[order(fused_volumes$volume_id,fused_volumes$variable),]
  rownames(fused_volumes) <- NULL
  return(fused_volumes)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unoverlap_hyperrectangles <- function(volumes) {
  partitioned_space <- build_fully_partitioned_space(volumes)
  new_volumes <- generate_volumes_from_partitioned_space(partitioned_space)
  solution <- prune_uncovering_volumes(new_volumes, volumes)
  fuse_abutted_hyperrectangles(solution, volumes)
}
