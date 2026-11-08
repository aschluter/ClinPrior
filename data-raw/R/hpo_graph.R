# Build the HPO graph (g1) from hp.obo.
#
# Same result as the original createHPOigraph(): one directed edge
# parent -> child for every "is_a" line, for the term and each of its
# alt_ids, plus, for obsolete terms, an edge from the parent of the first
# term that names the obsolete term in "consider" (as a child) to the
# obsolete term. Vertex order is kept (it defines the column order of HPOadj).
buildHPOgraph <- function(obo_file) {
  lines <- readLines(obo_file, warn = FALSE)
  n     <- length(lines)
  # The original loop stops before the last line
  lines <- lines[-n]
  idx   <- seq_along(lines)
  
  # Lines inside a [Term] block: the state switches on at "[Term]" and off
  # at each blank line (as the "bol" flag of the original loop)
  event   <- ifelse(lines == "[Term]", 1L, ifelse(lines == "", 0L, NA_integer_))
  seen    <- cumsum(!is.na(event))
  state   <- c(0L, event[!is.na(event)])[seen + 1L]
  in_term <- state == 1L
  block   <- cumsum(lines == "[Term]")
  
  is_id    <- in_term & grepl("^id: HP:[0-9]{7}|^alt_id: HP:[0-9]{7}", lines, perl = TRUE)
  is_isa   <- in_term & grepl("^is_a: HP:[0-9]{7}", lines, perl = TRUE)
  is_cons  <- in_term & grepl("^consider: HP:[0-9]{7}", lines, perl = TRUE)
  
  ids <- data.frame(block = block[is_id], line = idx[is_id],
                    term = gsub("alt_id: |id: ", "", lines[is_id]),
                    stringsAsFactors = FALSE)
  isa <- data.frame(block = block[is_isa], line = idx[is_isa],
                    parent = regmatches(lines[is_isa],
                                        regexpr("HP:[0-9]{7}", lines[is_isa])),
                    stringsAsFactors = FALSE)
  cons <- data.frame(block = block[is_cons], line = idx[is_cons],
                     consider = regmatches(lines[is_cons],
                                           regexpr("HP:[0-9]{7}", lines[is_cons])),
                     stringsAsFactors = FALSE)
  
  # is_a edges: each is_a line with every id/alt_id seen before it in the block
  m <- merge(isa, ids, by = "block", suffixes = c(".isa", ".id"))
  m <- m[m$line.id < m$line.isa, ]
  m <- m[order(m$line.isa, m$line.id), ]
  descendant <- cbind(m$parent, m$term)            # parent, child
  
  # Obsolete terms: first id of the block, paired with each "consider" term
  first_id <- ids[!duplicated(ids$block), c("block", "term")]
  ob <- merge(cons, first_id, by = "block")
  ob <- ob[order(ob$line), ]
  pos <- match(ob$consider, descendant[, 2])
  ensemble <- rbind(descendant, cbind(descendant[pos, 1], ob$term))
  
  igraph::graph_from_edgelist(ensemble, directed = TRUE)
}