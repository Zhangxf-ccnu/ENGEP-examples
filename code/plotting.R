plot_exp_pattern <- function(n,locations,gene_expressions,pt_size=1,n_col=3) {
  library(ggplot2)
  library(gridExtra)
  df <- mtcars[1:n, ]
  plots <- lapply(1:n, function(i) {
    df <- locations
    df$value <- gene_expressions[,i]
    ggplot(df, aes(X0, X1, color = value)) +
      geom_point(size = pt_size) +
      scale_color_viridis_c() +
      coord_equal()+
      theme(legend.position="none")+
      ggtitle(colnames(gene_expressions)[i])+
      theme(panel.grid =element_blank()) +    
      theme(axis.text = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())
  })
  grid.arrange(grobs = plots, nrow = ceiling(n / n_col))
}

plot_score<-function(maxcor,threshold){
  df <- data.frame(correlation = maxcor)
  ggplot(data=df,aes(x=correlation)) + 
    geom_histogram(binwidth=0.05,fill="#69b3a2", 
                   color="#e9ecef")+
    geom_vline(xintercept = threshold)+ 
    theme_bw()+
    labs(x="correlation between genes and archetypes",y="")
}