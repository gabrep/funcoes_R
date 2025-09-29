###################################################
##Gene set enrichment analysis
#de genes de acordo com rankeamento por seu valor de expressão
#extraido a partir do pacote limma, função topTable##
###################################################

##Pacotes utilizados##
#"clusterProfiler"
#'msigdbr'
#'dplyr'



run_GSEA <- function(res, #dataframe de resultado de DEGs
                     org = 'Homo sapiens',
                     genes = 'SYMBOL', #coluna dos resultados que traz o identificador dos genes
                     bases = c('HALLMARK', 'KEGG', 'REACTOME', 'WIKIPATHWAYS') #selecionar as bases a se utilizar para rodar o GSEA
){
  #Converter nomes de colunas de resultados vindos de DESEq para manter padrao igual ao limma
  if (class(res) == 'DESeqResults') {
    res <- as.data.frame(res)
    res$SYMBOL <- rownames(res) 
    
    if ('log2FoldChange' %in% colnames(res)) {
      res <- dplyr::rename(res, logFC = log2FoldChange)
    }
    
    if ('padj' %in% colnames(res)) {
      res <- dplyr::rename(res, adj.P.Val = padj)
    }
   
  }
  
  if (any(colnames(res) %in% c('padj', 'log2FoldChange', 'ENTREZ'))) {
    res <- as.data.frame(res)
    
    #Renomear
    if ('log2FoldChange' %in% colnames(res)) {
      res <- dplyr::rename(res, logFC = log2FoldChange)
    }
    
    if ('padj' %in% colnames(res)) {
      res <- dplyr::rename(res, adj.P.Val = padj)
    }
    
    if ('ENTREZ' %in% colnames(res)) {
      res <- dplyr::rename(res, ENTREZ = ENTREZID)
    }
    
  } 
  
  ##Download das informações para enriquecimento, de acordo com as bases selecionadas
  
  CPs <- NULL
  
  if (org == 'Homo sapiens') {
    
    for (i in seq_along(bases)) {
      #HALLMARKS sao uma subcategoria separada, nao pertencem a categoria C2 do msigdb
      if (bases[i] == 'HALLMARK') { 
        H <- msigdbr::msigdbr(species = org, category = "H")
        CPs[[bases[i]]] <- H
        } else {
          ## kegg, reactome e wikipathways sao subcategorias dentro da categoria C2.
            # para acessar as subcategorias é preciso o prefixo 'CP:' antes do nome da subcategoria
          subcat <- paste0('CP:', bases[i])
          CPs[[bases[i]]] <- msigdbr::msigdbr(species = org, category = "C2", subcategory = subcat)
        }
    }
  }
  
  if (org == 'Mus musculus') {
    CPs <- msigdbr::msigdbr(species = 'Mus musculus', category = "C2", subcategory = 'CP')
  }
    CPs_db <- dplyr::bind_rows(CPs)
    CPs_db <- dplyr::select(CPs_db, gs_name, gene_symbol)
  
    #remover os genes duplicados, mantendo apenas os de valores mais extremos (maior valor para logFC positivo, menor para logFC negativo)
    
    res <- res[order(abs(res$logFC), decreasing=T),]
    res <- res[!duplicated(res[[genes]]),] %>% na.omit()
    
    
  ### Criar rank de valores de expressão ##
    ord <- order(res$logFC, decreasing = T) #ordena os genes do maior logFC para o menor
    res_ord <- res[ord,] #dataframe de resultados com a ordem estabelecida
    
    # Remover NAs
    res_ord <- res_ord[!is.na(res_ord[[genes]]),] #remove os genes sem simbolo mantendo a ordem
    
    #criar o vetor de ranks
    ranks <- res_ord$logFC
    names(ranks) <- res_ord[[genes]]
    

  ### GSEA para os termos selecionados
  res.gsea <- GSEA(geneList = ranks, TERM2GENE = CPs_db)
  res.gsea <- res.gsea@result
  res.gsea$database <- sub('_.*', '', res.gsea$ID)
  return(res.gsea)
}  
