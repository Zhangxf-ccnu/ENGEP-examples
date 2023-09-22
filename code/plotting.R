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


#####

plot_score<-function(maxcor,density_estimation,x_value_at_min_density){
  hist_data <- data.frame(
    x = maxcor  
  )
  line_data <- data.frame(
    x = density_estimation$x,  
    y = density_estimation$y  
  )
  mmax = x_value_at_min_density+0.05
  mmin = x_value_at_min_density-0.05
  ggplot() +
    geom_histogram(data = hist_data, aes(x = x,y=..density..),
                   bins = 20, fill = "#20B2AA", color = "#e9ecef") +
    geom_line(data = line_data, aes(x = x, y = y),linetype = 1,size=1.3,color = '#FF7F50') +
    geom_vline(xintercept = mmin, linetype = "dashed",color = "#FF7F50",size=1.3) +
    geom_vline(xintercept = mmax, linetype = "dashed",color = "#FF7F50",size=1.3) +
    geom_point(data = data.frame(x = x_value_at_min_density, y = c(0)), aes(x, y), color = "blue", size = 3, shape = 1) +
    labs(
      x = "Likelihood scores",
      y = "Density"
    ) +
    theme_minimal() +#theme_bw
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    )
  
}


