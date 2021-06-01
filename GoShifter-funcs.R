#' shifter
#'
#' @param x vector of length n
#' @param n integer specifying how many places to shift elements in x
#'
#' @return vector x with each element shifted n places to the left
#' @export
#'
#' @examples
#' a = c(1,2,3,4,5)
#' shifter(a, 1)
#' shifter(a, 3)
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

#' my_GoShifter
#'
#' @param p vector of GWAS p-values for n SNPs
#' @param q vector of binary values corresponding to each of the n SNPs
#'
#' @return data.frame of p1, n1, p0, n0 values for the original SNP-annotation configuration (row 1) and for 1000 random permutations (shift) of q (n0 (n1) is the number of  SNPs with q=0 (q=1) and p1 (p0) is the mean -log10 p-val for SNPs (not) overlapping the annotation)
#' @export
#'
#' @examples
#' my_p <- runif(50, 0, 1)
#' my_q <- sample(c(0,1), 50, replace = T)
#' out <- my_GoShifter(p = my_p, q = my_q)
my_GoShifter <- function(p, q){
  
  n1 = length(which(q != 0))
  n0 = length(which(q == 0))
  
  if(n1 == 0){
    p1 = 0
  } else p1 = mean(-log10(p[which(q != 0)]))
  
  if(n0 == 0){
    p0 = 0
  } else p0 = mean(-log10(p[which(q == 0)]))
  
  initial_stat <- data.frame(shift = 0,  p1, n1, p0, n0)
  
  datalist = list()
  
  for (i in 1:1000) {
    
    shift_i = sample.int(length(q), 1)
    q_shift = shifter(q, shift_i)
    
    n1 <- length(which(q_shift != 0))
    n0 <- length(which(q_shift == 0))
    
    if(n1 == 0){
      p1 = 0
    } else p1 = mean(-log10(p[which(q_shift != 0)]))
    
    if(n0 == 0){
      p0 = 0
    } else p0 = mean(-log10(p[which(q_shift == 0)]))
    
    df = data.frame(shift = shift_i, p1, n1, p0, n0)
    
    datalist[[i]] <- df
  }
  
  shift_df = do.call(rbind, datalist)
  
  rbind(initial_stat, shift_df)
  
}

