
# funcoes_R

Repositório com funções em R para uso em análises de bioinformática, desde manipulação de tabelas a criação de gráficos

<!-- badges: start -->
## Descrição

Este repositório contém scripts de funções criadas para uso pessoal ou em múltiplos projetos. Cada função está em seu próprio arquivo `.R` dentro do repositório e pode ser carregada via `source()`.
<!-- badges: end -->


## Funções incluídas

- **run_ORA**: run_ORA(res,  
                       lfc_threshold = 1,  
                       padj_threshold = 0.05,  
                       enrich=c('GO', 'KEGG', 'WIKIPATHWAYS'),  
                       keyType = 'ENTREZID',  
                       OrgDb = org.Hs.eg.db,  
                       organismKEGG = 'hsa', organismWP = 'Homo sapiens',  
                       ont = 'all')    

  
  Realiza o enriquecimento de Over Representation Analysis de DEGs. São aceitos resutados obtidos por DESeq2 (results(DESeq(dds))) ou limma (topTable(fit2)).
  A função filtra e traz resultados separados entre genes UP e DOWN  
  Pode ser inserido DEGs de diferentes estudos/datasets em uma lista. 
  Exemplo de aplicação para mais de um dataset ao mesmo tempo:  
  lapply(lista_degs_datasets, function(i) {run_enrichGO(res = i})
  
  
- **run_GSEA**: function(res, 
                     genes = 'SYMBOL',  
                     logfc = 'logFC',   
                     bases = c('HALLMARK', 'KEGG', 'REACTOME', 'WIKIPATHWAYS'))  
    
  Realiza enriquecimento de Gene Set Analysis de dados de expressao. São aceitos resutados obtidos por DESeq2 (results(DESeq(dds))) ou limma (topTable(fit2)).
  A função recebe a lista de genes com valor de expressão e gera o rank de genes para utilização no GSEA.  
  Pode ser escolhido realizar GSEA tradicional com uma das bases disponíveis (Hallmark, KEGG, Reactome ou Wiki Pathways), ou manter mais de uma para utilização de todos os termos em uma única análise
    
     
- **color_dend**: color_dend(exp, grupos, cores_dend, col_branches)    
  Gera dendograma de distancias euclidianas com personalização de cores para cada grupo.  
  col_branches deve ser utilizado como T/F para colorir os ramos do dendrograma (mesmas cores de cores_dend)

