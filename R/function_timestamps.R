
#' timestamp.mdy
#'
#' @return mmddyyyy timestamp for filenames
#' @export
#'
#' @examples
timestamp.mdy = function(){
    format(Sys.time(), "%m%d%y")    
}
