
# funcoes_R

Repositório com funções em R para uso em análises de bioinformática, desde manipulação de tabelas a criação de gráficos

<!-- badges: start -->
## Descrição

Este repositório contém scripts de funções criadas para uso pessoal ou em múltiplos projetos. Cada função está em seu próprio arquivo `.R` dentro do repositório e pode ser carregada via `source()`.
<!-- badges: end -->


## Funções incluídas
- **run_enrichGO**
  Realiza o enriquecimento por Gene Ontology de DEGs geradas a partir do pacote limma e selecionadas pela função _topTable_.
  Pode ser inserido DEGs de diferentes estudos/datasets em uma lista. 
  A função filtra e tráz resultados separados entre genes UP e DOWN

