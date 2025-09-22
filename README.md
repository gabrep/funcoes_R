
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
                       GO = T,
                       KEGG = T,
                       WIKIPATHWAYS = T,
                       keyType = 'ENTREZID',
                       OrgDb = org.Hs.eg.db,
                       organismKEGG = 'hsa', organismWP = 'Homo sapiens',
                       ont = 'all')  

  Realiza o enriquecimento de Over Representation Analysis de DEGs geradas a partir do pacote limma e selecionadas pela função _topTable_.
  Pode ser inserido DEGs de diferentes estudos/datasets em uma lista. 
  A função filtra e traz resultados separados entre genes UP e DOWN  
  Exemplo de aplicação para mais de um dataset ao mesmo tempo:  
  lapply(lista_degs_datasets, function(i) {run_enrichGO(res = i})
    
- **color_dend**: color_dend(exp, grupos, cores_dend)  
  Gera dendograma de distancias euclidianas com personalização de cores para cada grupo.  

