#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @return The altered object.
'%na%' <- function(a, b) {
	if (is.null(a) || length(a) == 0) return(b)
	if (is.null(b) || length(b) == 0) return(a)
	return(ifelse(is.na(a), b, a))
}

#' Uses b if a is NULL
#' @param a,b objects to be tested or coerced.
#' @return An un-null object.
'%null%' <- function(a, b) {
	if (is.null(a)) return(b)
	return (a)
}





