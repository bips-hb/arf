dsofttruncnorm <- function(x,a,b,mean,sd, offset = 10^(-7), scale = 10^(-10)) {

  plogis(x,a + offset, scale) * dnorm(x, mean, sd)/sum(pnorm(c(b,a), mean, sd)*c(1,-1)) * plogis(x,b + offset,scale,lower.tail = F)
  
}

softquadratic_discrete <- function(x, val, prob) {
  -sum(plogis(x,val - 0.1, 10^(-10))* (-(x-val)^2 + prob) * plogis(x,val + 0.1, 10^(-10), lower.tail = F))
}

###

optim(c(3), softquadratic_discrete, val=c(3,4,5), prob = c(1,2,1.9),control = list(fnscale = -1), method = "SANN")
optim(c(-100), dsofttruncnorm, a = -1, b = 5, mean = 6, sigma = 1, control = list(fnscale = -1, maxit = 100000), method = "SANN")
softquadratic_discrete(c(3,0), 3, 1)


optimg(c(3.95), softquadratic_discrete, val=c(3,4,5), prob = c(1,2,1.9), method = "ADAM", verbose = T, full = T)

###

softlik_minimize <- function(x, params, arf, log = T) {
  -softlik(x, params, arf, log)
}

softlik(test, params, arf,log = F)
test <- data.frame(3.589176, 2.892007, 5.339714, 1.501191, 3)
names(test) <- names(iris)
x_ <- iris[1,]
lik(x_,params, arf, log= F)

start <- data.frame(0,0,0,0,0)
names(start) <- names(iris)

optimg(start, softlik_minimize, params= params, arf =arf, log = F, method = "ADAM")


softlik <- function(x, params, arf, log = T) {
  cnt <- copy(params$cnt)
  
  if(is.null(cnt$val)) {
    cnt[,val:=NA]
    cnt[,val := as.numeric(val)]
  }
  cat <- copy(params$cat)
  forest <- copy(params$forest)[,cvg := cvg/sum(cvg)]
  
  cat_cols <-params$meta[class == "factor", variable]
  cnt_cols <-params$meta[class != "factor", variable]
  
  # leaves <- as.data.table(predict(arf, x, type = 'terminalNodes')$predictions + 1)[,obs:=.I]
  # n <- nrow(leaves)
  # n_trees <- ncol(leaves) -1
  # leaves <- melt(leaves, id="obs", variable.name = "tree", value.name= "leaf")[, tree := rep(seq_len(n_trees),each = n)]
  # leaves <- merge(leaves,forest[,.(tree, leaf, f_idx)],sort=F, by=c("tree","leaf"), allow.cartesian=TRUE)[,.(obs,f_idx)]
  
  x_long <- data.table(variable = names(x), val = as.vector(unlist(x)))
  print(x_long)
  
  if (length(cat_cols)>0) {
    cat[,val:= as.numeric(as.factor(val))]
    cat_lik <- merge(cat, x_long, by = "variable", sort = F)
    cat_lik <- cat_lik[,softquadratic_discrete(val.y,val.x,prob), by = .(variable,f_idx)]
    setnames(cat_lik, "V1", "lik")

    #leaves <- cat_lik[,.N,by=.(obs,f_idx)][N == length(cat_cols),.(obs,f_idx)]
  } else {
    cat_lik <- data.table()
  }
  if (length(cnt_cols)>0){
    cnt_lik <- merge(cnt,x_long, sort=F, by= "variable")
    cnt_lik <- cnt_lik[,lik:= is.na(val.x)*dsofttruncnorm(val.y,a=min,b=max,mean=mu,sd=sigma) + (!is.na(val.x)&(val.y==val.x)), by = .(variable,f_idx)] ###########!!!!
    cnt_lik <- cnt_lik[,.(variable, f_idx,lik)]

    # cnt_lik <- cnt_lik[,check:=any(lik==0),by=.(obs,f_idx)][check==F,-"check"]
    # 
    # if (nrow(leaves)*length(cnt_cols) > nrow(cnt_lik)) {
    #   leaves <- unique(cnt_lik[,.(obs,f_idx)])
    #   cat_lik <- merge(leaves,cat_lik,sort=F,allow.cartesian = T, by=c("obs","f_idx"))
    # }
    
  } else {
    cnt_lik <- data.table()
  }
  
  lik <- rbind(cat_lik,cnt_lik)
  lik <- lik[,prod(lik),by = f_idx][,cvg := forest$cvg]
  lik <- as.vector(lik[,crossprod(V1,cvg)])
  
  if (log) {
    log(lik)
  } else {
    lik
  }
}

