run_enrichGO <- function(res,
                       lfc_threshold = 1,
                       padj_threshold = 0.05,
                       keyType = 'ENTREZID',
                       OrgDb = org.Hs.eg.db,
                       ont = 'all') {
  #Genes Up
  genes_up <- res %>% 
    filter(logFC > lfc_threshold, adj.P.Val < padj_threshold) %>% 
    pull(!!sym(keyType)) %>% 
    unique()
  
  #Genes down
  genes_down <- res %>% 
    filter(logFC < -lfc_threshold, adj.P.Val < padj_threshold) %>% 
    pull(!!sym(keyType)) %>% 
    unique()
  
  #Enrichment
  go_up <- NULL
  go_down <- NULL
  
  if(length(genes_up)>0) {
    go_up <- enrichGO(gene = genes_up,
                      OrgDb = OrgDb,
                      keyType = keyType,
                      ont = ont)
  }
  
  if(length(genes_down)>0) {
    go_down <- enrichGO(gene = genes_down,
                        OrgDb = OrgDb,
                        keyType = keyType,
                        ont = ont)
  }
  return(list(up = go_up, down = go_down))
}
