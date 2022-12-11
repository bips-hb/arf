# What to do when y is a column name?
# Problem: Sometimes pruning doesn't solve the no-variance problem b/c
# all samples in a leaf share the same value for a variable.

# Hexsticker: https://github.com/GuangchuangYu/hexSticker 

# Allow arbitrary queries for total evidence, marginals, MAP, etc.
# e.g. p(x1), p(x1, x2), p(x1 | x2). Problem: need to accommodate both
# equalities and inequalities on either side of the conditioning bar.
# Also conjunctions and disjunctions.

# Add conditional sampling option to forge

# Add methods for data imputation

# Other exponential family options, with normalization constant?

# Add KDE option?

# Idea: let each tree have its own synthetic dataset?

# Bayesian version with some prior? Could help with zero-density events when
# a level is unobserved, for instance, without being ruled out per se.

# Localized version, where we only regrow nodes that split too well

# Fully adversarial with alternating rounds at each node?

