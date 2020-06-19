
#' read.table.fix_xlsx_dates
#'
#' attempt to load file f and repair any gene names that excel has converted to dates
#'
#' @param f file path or data.frame
#' @param sep delimeter character for read.table, default is ,
#' @param header header logical for read.table, default is FALSE
#'
#' @return data.frame with dates converted to gene names
#' @export
#'
#' @examples
read.table.fix_xlsx_dates = function(f, sep = ",", header = FALSE){
    if(is.character(f)){
        df = read.table(f, sep = sep, header = header)    
    }else if(is.data.frame(f)){
        df = f
    }else{
        stop("f must be file or data.frame")
    }
    todo = which(sapply(df, class) %in% c("factor", "character"))
    for(i in todo){
        df[[i]] = fix_xlsx_dates(df[[i]])
    }
    if(!is.null(rownames(df))) rownames(df) = fix_xlsx_dates(rownames(df))
    df
}

#' fix_xlsx_dates
#' 
#' repair any gene names that excel has converted to dates
#'
#' @param input character or factor to repair 
#' @param targets date strings to search for
#' @param fixes replacement character when targets found.  must parallel targets.
#'
#' @return same class as input with dates converted back to gene names
#' @export
#'
#' @examples
fix_xlsx_dates = function(input, targets = c("Mar", "Sep"), fixes = c("March", "Sept")){
    stopifnot(class(input) %in% c("factor", "character"))
    stopifnot(length(targets) == length(fixes))
    factor_mode = is.factor(input)
    out = input
    if(factor_mode) input = levels(input)
    for(i in seq_along(targets)){
        t.string = targets[i]
        f.string = fixes[i]
        hit_idx = which(grepl(paste0("^[0-9]+-", t.string, "$"), input))
        hit_str = input[hit_idx]
        fix_str = paste0(f.string, sub(paste0("-", t.string), "", hit_str))
        input[hit_idx] = fix_str
    }
    if(factor_mode) levels(out) = input
    out
}

# f = "/slipstream/home/joeboyd/R/SF_AutoImmune_ssv/data/RNAseq/DESeq2_outputs/DF4_NvsM40_up+down.csv"
# sep = ","
# header = TRUE
