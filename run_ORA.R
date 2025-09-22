###################################################
##ORA de DEGs a partir do pacote limma, função topTable##
###################################################

##Pacotes utilizados##
  #"clusterProfiler"
  #"org.Hs.eg.db" 
  #"dplyr"


# Com esta função, o enriquecimento de Over Representation Analysis é realizado para uma lista
#de resultados (DEGs a partir do pacote limma). Pode-se realizar ORA via Gene Ontology, KEGG e WikiPathways
# A função separa entre resultados Up e Down, de acordo com os cut-offs de lfc_threshold e
#padj_threshold
# Se houver mais de uma lista de DEGs de diferentes datasets, pode-se criar uma lista
#contendo cada resultado, por exemplo:
#   list(res.1 = DEGs_dataset_1, res.2 = DEGs_dataset_2)



#Parametros da função enrichGO:
run_ORA <- function(res,
                       lfc_threshold = 1,
                       padj_threshold = 0.05,
                       GO = T,
                       KEGG = T,
                       WIKIPATHWAYS = T,
                       keyType = 'ENTREZID',
                       OrgDb = org.Hs.eg.db,
                       organismKEGG = 'hsa', organismWP = 'Homo sapiens',
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
  
  #############Enrichment############
  #############GO############
  if(GO == T) {
    
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
  }
  
  
  #############KEGG################
  if(KEGG == T){
    kegg_up <- NULL
    kegg_down <- NULL
    
    if (length(genes_up) > 0) {
      kegg_up <- enrichKEGG(genes_up, 
                            organism = organismKEGG, 
                            keyType = 'kegg')
    }
    
    if (length(genes_down) > 0) {
      kegg_down <- enrichKEGG(genes_down, 
                            organism = organismKEGG, 
                            keyType = 'kegg')
    }
  }

  
  ############REACTOME############
  if(WIKIPATHWAYS == T){
    wikipathways_up <- NULL
    wikipathways_down <- NULL
    
    if (length(genes_up) > 0) {
      wikipathways_up <- enrichWP(genes_up, 
                            organism = organismWP)
    }
    
    if (length(genes_down) > 0) {
      wikipathways_down <- enrichWP(genes_down, 
                              organism = organismWP)
    }
  }
  
  
  #guardar resultados
  return(list(go_up = go_up, go_down = go_down,
              kegg_up = kegg_up, kegg_down = kegg_down,
              wp_up = wikipathways_up, wp_down = wikipathways_down))
}
