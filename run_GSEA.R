###################################################
##Gene set enrichment analysis
#de genes de acordo com rankeamento por seu valor de expressão
#extraido a partir do pacote limma, função topTable##
###################################################

##Pacotes utilizados##
#"clusterProfiler"
#'msigdbr'
#'dplyr'



run_GSEA <- function(res, #dataframe de resultado de DEGs extraído com limma::topTable
                     genes = 'SYMBOL', #coluna dos resultados que traz o identificador dos genes
                     logfc = 'logFC', #coluna dos resultados que traz o valor de logFC
                     bases = c('HALLMARK', 'KEGG', 'REACTOME', 'WIKIPATHWAYS') #selecionar as bases a se utilizar para rodar o GSEA
){
  ##Download das informações para enriquecimento, de acordo com as bases selecionadas
  
  CPs <- NULL
  for (i in seq_along(bases)) {
    #HALLMARKS sao uma subcategoria separada, nao pertencem a categoria C2 do msigdb
    if (bases[i] == 'HALLMARK') { 
      H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      CPs[[bases[i]]] <- H
      } else {
        ## kegg, reactome e wikipathways sao subcategorias dentro da categoria C2.
          # para acessar as subcategorias é preciso o prefixo 'CP:' antes do nome da subcategoria
        subcat <- paste0('CP:', bases[i])
        CPs[[bases[i]]] <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = subcat)
      }
  }
    CPs_db <- dplyr::bind_rows(CPs)
    CPs_db <- dplyr::select(CPs_db, gs_name, gene_symbol)
  
    #remover os genes duplicados, mantendo apenas os de valores mais extremos (maior valor para logFC positivo, menor para logFC negativo)
    
    res <- res[order(abs(res[[logfc]]), decreasing=T),]
    res <- res[!duplicated(res[[genes]]),] %>% na.omit()
    
    
  ### Criar rank de valores de expressão ##
    ord <- order(res[[logfc]], decreasing = T) #ordena os genes do maior logFC para o menor
    res_ord <- res[ord,] #dataframe de resultados com a ordem estabelecida
    
    # Remover NAs
    res_ord <- res_ord[!is.na(res_ord[[genes]]),] #remove os genes sem simbolo mantendo a ordem
    
    #criar o vetor de ranks
    ranks <- res_ord[[logfc]]
    names(ranks) <- res_ord[[genes]]
    

  ### GSEA para os termos selecionados
  res.gsea <- GSEA(geneList = ranks, TERM2GENE = CPs_db)
  res.gsea <- res.gsea@result
  res.gsea$database <- sub('_.*', '', res.gsea$ID)
  return(res.gsea)
}  
