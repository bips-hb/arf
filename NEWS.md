# arf 0.2.2
* Faster and vectorized conditional sampling
* Use min.bucket argument from ranger to avoid pruning if possible
* Option to sample NAs in generated data if original data contains NAs
* Stepsize in forge() to reduce memory usage
* Option for local and global finite bounds

# arf 0.2.0
* Vectorized adversarial resampling
* Speed boost for compiling into a probabilistic circuit
* Conditional densities and sampling
* Bayesian solution for invariant continuous data within leaf nodes
* New function for computing (conditional) expectations
* Options for missing data

# arf 0.1.3
* Speed boost for the adversarial resampling step 
* Early stopping option for adversarial training
* alpha parameter for regularizing multinomial distributions in forde
* Unified treatment of colnames with internal semantics (y, obs, tree, leaf)