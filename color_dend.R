###################################################
##Gera dendograma personalizado para os resultados##
###################################################

# Esta função depende do pacote dendextend

# Variaveis:
  ##exp = matrix com valores de expressão, o nome das colunas será utilizada como identificador no dendrograma
  ##grupos = dataframe e coluna de pdata que contem o agrupamento (ex. pdata$treatment)
  ##cores_dend = personalizar cores


color_dend <- function(exp, grupos, cores_dend, col_branches) {
  
  
  require(dendextend)
  #cores genericas para o dendograma
  cores_dend <- c('firebrick1', 'deepskyblue2', 'gold', 'violetred', 'purple', 'chocolate', 'green', 'black')
  
  #conferir grupos presentes nos metadados
  grupos_unicos <- unique(grupos)
  
  
  #checar se o numero de cores é suficiente
  if (length(cores_dend) < length(grupos_unicos)) {
    stop('Número de cores insuficientes')
  }
  
  # calculo de distancia e clustering
  dist <- dist(t(exp))
  hc <- hclust(dist)
  dend <- as.dendrogram(hc)
  
  #dendograma
  for (i in seq_along(grupos_unicos)) {
    dend <- dend %>% color_labels(labels = colnames(exp[, which(grupos == grupos_unicos[i])]),
    col = cores_dend[i])
      
  }
  if (col_branches == TRUE) {
  dend <- color_branches(dend, k=length(grupos_unicos), col = cores_dend)
  }
  plot(dend)
  
  #legenda
  legend('topright', legend = grupos_unicos, col = cores_dend, pch=20)
}