# #converting ucsc id to gene symbol
# gene_symbol <- read.delim(system.file("extdata", "gene_symbol.txt", package = "RTDetect"),
#                           header=TRUE, comment.char="#")
# gene_symbol <- bind_rows(gene_symbol, data.frame(kgID=NA, geneSymbol=NA))
# 
# #insSite
# txs <- RT$txs
# gene_syms <- lapply(txs, function(x) gene_symbol$geneSymbol[gene_symbol$kgID %in% x])
# gene_syms <- lapply(gene_syms, unique)
# unlist(gene_syms) %>% unique()
# 
# 
# 
# RT$rt$gene_syms <- .txs2genesym(RT$rt$txs)
# unlist(.txs2genesym(RT$insSite$txs)) %>% length()
# 
# .txs2genesym(RT$insSite$txs) %>% as_tibble()
# 
# 
# lapply(RT$rt$gene_syms, function(x) y %in% x)
# lapply(X = RT$rt$gene_syms, FUN = function(x) RT$rt[x %in% RT$rt$gene_syms])
# 
# 
# i <- lapply(unique(unlist(RT$rt$gene_syms)), function(gs)
#     sapply(RT$rt$gene_syms, function(x) gs %in% x))
# l <- lapply(i, function(i) RT$rt[i])
# names(l) <- c('a','b','3','4','5')
# ll <- lapply(l, function(i) list(rt=l[[1]]))
# 
# RT$rt$gene_symbol <- .txs2genesym(RT$rt$txs)
# RT$insSite$gene_symbol <- .txs2genesym(RT$insSite$txs)
# l_gene_symbol <- unique(c(unlist(.txs2genesym(RT$rt$txs)), unlist(.txs2genesym(RT$insSite$txs))))
# 
# #RT GRangesList by gene
# rt.gr.idx <- lapply(l_gene_symbol, function(gs) 
#     sapply(RT$rt$gene_syms, function(x) gs %in% x))
# rt.grlist <- setNames(lapply(rt.gr.idx, function(i) list(rt=RT$rt[i])), l_gene_symbol)
# #names(rt.grlist) <- l_gene_symbol
# 
# #InsSite GRangesList by gene
# insSite.gr.idx <- lapply(l_gene_symbol, function(gs) 
#     sapply(RT$insSite$gene_syms, function(x) gs %in% x))
# #insSite.grlist <- lapply(insSite.gr.idx, function(i) RT$insSite[i])
# insSite.grlist <- setNames(lapply(insSite.gr.idx, function(i) list(insSite=RT$insSite[i])), l_gene_symbol)
# #names(insSite.grlist) <- l_gene_symbol
# 
# ll <- pc(rt.grlist, insSite.grlist)
# ll <- lapply(ll, function(x) setNames(x, c('RT', 'insSite')))
# #ll <- regroup(c(rt.grlist, insSite.grlist), l_gene_symbol)
# 
# rt <- lapply(rt.grlist, function(l) list(rt=l[[1]]))
# ins <- lapply(insSite.grlist, function(l) list(insSite=l[[1]]))
# 
# 
# 
# class(GRangesList(insSite.grlist))
# 
