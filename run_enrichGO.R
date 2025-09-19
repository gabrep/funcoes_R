###################################################
##GO de DEGs a partir do pacote limma, função topTable##
###################################################

##Pacotes utilizados##
  #"clusterProfiler"
  #"org.Hs.eg.db" 
  #"dplyr"


# Com esta função, o enriquecimento por Gene Ontology é realizado para uma lista
#de resultados (DEGs a partir do pacote limma).
# A função separa entre resultados Up e Down, de acordo com os cut-offs de lfc_threshold e
#padj_threshold
# Se houver mais de uma lista de DEGs de diferentes datasets, pode-se criar uma lista
#contendo cada resultado, por exemplo:
#   list(res.1 = DEGs_dataset_1, res.2 = DEGs_dataset_2)
#keyType pode ser especificado ao chamar a função com lapply


#Parametros da função enrichGO:
run_enrichGO <- function(res,
                       lfc_threshold = 1,
                       padj_threshold = 0.05,
                       keyType = 'ENTREZID',
                       OrgDb = org.Hs.eg.db,
                       ont = 'all') {
  
  #Separar os resultados em Up e Down
  #Genes Up
  genes_up <- res %>% 
    filter(logFC > lfc_threshold, adj.P.Val < padj_threshold) %>% 
    pull(!!sym(keyType)) %>% #!!sym indica que a string de entrada deve ser entendida como um nome de coluna
    unique()
  
  #Genes down
  genes_down <- res %>% 
    filter(logFC < -lfc_threshold, adj.P.Val < padj_threshold) %>% 
    pull(!!sym(keyType)) %>% 
    unique()
  
  #Enrichment
  #criar objeto para armazenar resultados
  go_up <- NULL
  go_down <- NULL
  
  #GO de genes up regulados
  if(length(genes_up)>0) {
    go_up <- enrichGO(gene = genes_up,
                      OrgDb = OrgDb,
                      keyType = keyType,
                      ont = ont)
  }
  
  #GO de genes down regulados
  if(length(genes_down)>0) {
    go_down <- enrichGO(gene = genes_down,
                        OrgDb = OrgDb,
                        keyType = keyType,
                        ont = ont)
  }
  
  #armazenar resultados
  return(list(up = go_up, down = go_down))
}
