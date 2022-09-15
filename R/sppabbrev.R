# x = scientific names
# x=gsub(" ", "_", meta$Accepted_sciname)
# output genome species abbreviates (e.g., "homSap" for Homo sapiens)
sppAbbrev <- function(x) {
	tmp <- str_split(x, " |_")
	paste0(tolower(substring(sapply(tmp, "[[", 1), 1, 3)),
		   toupper(substring(sapply(tmp, "[[", 2), 1, 1)),
		   substring(sapply(tmp, "[[", 2), 2, 3))
}
