#' @keywords internal
#' @noRd

# Per-tree worker for adversarial_rf()'s post-hoc leaf pruning (enforces
# min_node_size w.r.t. real data). Base-R only (no arf/data.table internals), so
# mirai daemons need no package loaded. Returns the pruned child.nodeIDs for one
# tree: a list of two integer vectors (left, right children).
arf_prune_tree <- function(tree, child_nodeIDs, pred, min_node_size) {
  # Nodes to prune are leaves which contain fewer than min_node_size real samples
  out <- child_nodeIDs[[tree]]
  leaves <- which(out[[1]] == 0L)
  to_prune <- leaves[!(leaves %in% which(tabulate(pred[, tree]) >= min_node_size))]
  while (length(to_prune) > 0) {
    if (1 %in% to_prune) {
      # Never prune the root
      break
    }
    for (tp in to_prune) {
      # Find parent
      parent <- which((out[[1]] + 1L) == tp)
      if (length(parent) > 0) {
        # If node to prune (tp) is the left child of parent, replace left child with right child
        out[[1]][parent] <- out[[2]][parent]
      } else {
        # If node to prune (tp) is the right child of parent, replace right child with left child
        parent <- which((out[[2]] + 1L) == tp)
        out[[2]][parent] <- out[[1]][parent]
      }
    }
    # If both children of a parent are to be pruned, prune the parent in the next round
    # This happens if both children have been pruned
    to_prune <- which((out[[1]] + 1L) %in% to_prune)
  }
  out
}
