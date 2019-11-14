#' Title
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
sampleCap = function(x, n = 500){
    n = min(n, length(unique(x)))
    out = sample(unique(x), n)
    if(is.factor(out)) out = as.character(out)
    out
}