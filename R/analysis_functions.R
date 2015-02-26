#' Bi-partite network analysis tools
#'
#' This function analyzes a bi-partite network, such as a Transcription factor to gene network derived from the PANDA algorithm.
#'
#' @param net1 starting network, a genes by transcription factors data.frame with scores for confidence in the existence of edges between
#' @param net2 final network, a genes by transcription factors data.frame with scores for confidence in the existence of edges between
#' @keywords keywords
#' @export
#' @examples
#' data(yeast.panda)
#' t.matrix <- transformation.matrix(yeast.panda$cell.cycle, yeast.panda$stress.response)
#' hcl.heatmap.plot(t.matrix, method="pearson")
transformation.matrix <- function(network.1, network.2, by.genes=F,standardize=T, remove.diagonal=T){
  if(is.list(network.1)&&is.list(network.2)){
    if(by.genes){
      net1 <- t(network.1$reg.net)
      net2 <- t(network.2$reg.net)
    } else {
      net1 <- network.1$reg.net
      net2 <- network.2$reg.net
    }
  } else if(is.matrix(network.1)&&is.matrix(network.2)){
    if(by.genes){
      net1 <- t(network.1)
      net2 <- t(network.2)
    } else {
      net1 <- network.1
      net2 <- network.2
    }
  } else {
    stop("Networks must be lists or matrices")
  }
  #gene.trans.matrix <- svd(net2)$v %*% diag(1/svd(net2)$d) %*% t(svd(net2)$u) %*% net1
  tf.trans.matrix   <- svd(t(net2))$v %*% diag(1/svd(t(net2))$d) %*% t(svd(t(net2))$u) %*% t(net1)
  if (standardize){
    tf.trans.matrix <- apply(tf.trans.matrix, 1, function(x){
      #   x.zero <- (x-mean(x))
      x/sum(abs(x))
    })
  }
  if (remove.diagonal){
    diag(tf.trans.matrix) <- 0
  }
  # Add column labels
  colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
  tf.trans.matrix
}

#' Sum of squared off-diagonal mass
#'
#' This function calculates the off-diagonal sum of squared mass for a transition matrix
#'
#' @param tm a transition matrix for two bipartite networks
#' @keywords keywords
#' @export
#' @examples
#' data(yeast.panda)
#' t.matrix <- transformation.matrix(yeast.panda$cell.cycle, yeast.panda$stress.response)
#' ssodm(t.matrix)
ssodm <-  function(tm){
  diag(tm)<-0
  sort(apply(tm,1,function(x){t(x)%*%x}))
}

#' Transformation matrix plot
#'
#' This function plots a hierachically clustered heatmap and corresponding dendrogram of a transaction matrix
#'
#' @param net1 starting network, a genes by transcription factors data.frame with scores for confidence in the existence of edges between
#' @param method distance metric for hierarchical clustering.  Default is "Pearson correlation"
#' @keywords keywords
#' @export
#' @examples
#' data(yeast.panda)
#' t.matrix <- transformation.matrix(yeast.panda$cell.cycle, yeast.panda$stress.response)
#' hcl.heatmap.plot(t.matrix, method="pearson")
hcl.heatmap.plot <- function(x, method="pearson"){
  if(method=="pearson"){
    dist.func <- pearson.dist
  } else {
    dist.func <- dist
  }
  # x <- as.matrix(scale(mtcars))
  x <- scale(x)
  dd.col <- as.dendrogram(hclust(dist.func(x)))
  col.ord <- order.dendrogram(dd.col)

  dd.row <- as.dendrogram(hclust(dist.func(t(x))))
  row.ord <- order.dendrogram(dd.row)

  xx <- x[col.ord, row.ord]
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  df$Var1 <- xx_names[[1]]
  df$Var1 <- with(df, factor(Var1, levels=Var1, ordered=TRUE))
  mdf <- melt(df)


  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)

  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  ### Set up a blank theme
  theme_heatmap <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  ### Create plot components ###
  # Heatmap
  p1 <- ggplot(mdf, aes(x=variable, y=Var1)) +
    geom_tile(aes(fill=value)) + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

  # Dendrogram 1
  p2 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_none + theme(axis.title.x=element_blank())

  # Dendrogram 2
  p3 <- ggplot(segment(ddata_y)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    coord_flip() + theme_none

  ### Draw graphic ###

  grid.newpage()
  print(p1, vp=viewport(0.80, 0.8, x=0.400, y=0.40))
  print(p2, vp=viewport(0.73, 0.2, x=0.395, y=0.90))
  print(p3, vp=viewport(0.20, 0.8, x=0.910, y=0.43))
}

#' Principal Components plot of transformation matrix
#'
#' This function plots the first two principal components for a transaction matrix
#'
#' @param tm a transition matrix for a bipartite network
#' @param title The title of the plot
#' @param clusters A vector indicating the colors to be plotted for each node
#' @param alpha A vector indicating the level of transparency to be plotted for each node
#' @keywords keywords
#' @export
#' @examples
#' data(yeast.panda)
#' t.matrix <- transformation.matrix(yeast.panda$cell.cycle, yeast.panda$stress.response)
#' p.values <- runif(ncol(t.matrix))  # Generate a uniform random to simulate p.values
#' clusters <- kmeans(t.matrix,3)$cluster # Color the nodes according to cluster membership
#' pca.plot(t.matrix,title="PCA Plot of Transition - Cell Cycle vs Stress Response", clusters=clusters,alpha=p.values)
pca.plot <-  function(tm, title="PCA Plot of Transition", clusters=1, alpha=1){
  tm.pca <- princomp(tm)
  odsm <- ssodm(tm)
  odsm.scaled <- 2*(odsm-mean(odsm))/sd(odsm)+4
  scores.pca <- as.data.frame(tm.pca$scores)
  scores.pca <- cbind(scores.pca,'node.names'=rownames(scores.pca))
  ggplot(data = scores.pca, aes(x = Comp.1, y = Comp.2, label = node.names)) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_text(size = odsm.scaled, alpha=alpha, color=clusters) +
    ggtitle(title)
}

#' This function plots the Off diagonal mass of an observed Transition Matrix compared to a set of null TMs
#'
#' @param tm.obs The observed transition matrix
#' @param tm.null A list of null transition matrices
#' @keywords keywords
#' @export
#' @examples
#' example1
ssodm.plot <- function(tm.obs, tm.null, sort.by.sig=F, rescale=F, plot.title=NA, highlight.tfs=NA){
  if(is.na(plot.title)){
    plot.title <- "SSODM observed and null"
  }
  num.iterations <- length(tm.null)
  # Calculate the off-diagonal squared mass for each transition matrix
  null.SSODM <- lapply(tm.null,function(x){
    apply(x,1,function(y){t(y)%*%y})
  })
  null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
  null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))

  ssodm <- apply(tm.obs,1,function(x){t(x)%*%x})

  # Get p-value (rank of observed within null ssodm)
#   p.values <- sapply(1:length(ssodm),function(i){
#     1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
#   })
  p.values <- 1-pnorm(sapply(1:length(ssodm),function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }))

# Process the data for ggplot2
  combined.mat <- cbind(null.ssodm.matrix, ssodm)
  colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")


  if (rescale){
    combined.mat <- t(apply(combined.mat,1,function(x){
      (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
    }))
    x.axis.order <- rownames(tm.null[[1]])[order(p.values)]
    x.axis.size  <- 10 # pmin(15,7-log(p.values[order(p.values)]))
  } else {
    x.axis.order <- rownames(tm.null[[1]])
    x.axis.size  <- pmin(15,7-log(p.values))
  }

  null.SSODM.melt <- melt(combined.mat)[,-1][,c(2,1)]
  null.SSODM.melt$TF<-rep(rownames(tm.null[[1]]),num.iterations+1)

  ## Plot the data
  ggplot(null.SSODM.melt, aes(x=TF, y=value))+
    geom_point(aes(size=1,color=factor(Var2),alpha = .5*as.numeric(factor(Var2)))) +
    scale_color_manual(values = c("blue", "red")) +
    scale_x_discrete(limits = x.axis.order ) +
    theme(legend.title=element_blank(),axis.text.x = element_text(colour = 1+x.axis.order%in%highlight.tfs, angle = 90, hjust = 1, size=x.axis.size,face="bold")) +
    ylab("Sum of Squared Off-Diagonal Mass") +
    ggtitle(plot.title)
}

#' Calculate p-values for a tranformation matrix
#'
#' This function calculates the significance of an observed transition matrix given a set of null transition matrices
#'
#' @param tm.obs The observed transition matrix
#' @param tm.null A list of null transition matrices
#' @param method one of 'z-score' or 'non-parametric'
#' @keywords keywords
#' @export
#' @examples
#' example1
calculate.tm.p.values <- function(tm.obs, tm.null, method="z-score"){
  num.iterations <- length(tm.null)
  # Calculate the off-diagonal squared mass for each transition matrix
  null.SSODM <- lapply(tm.null,function(x){
    apply(x,1,function(y){t(y)%*%y})
  })
  null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
  null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))

  ssodm <- apply(tm.obs,1,function(x){t(x)%*%x})

  # Get p-value (rank of observed within null ssodm)
  if(method=="non-parametric"){
    p.values <- sapply(1:length(ssodm),function(i){
      1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
    })
  } else if (method=="z-score"){
    p.values <- pnorm(sapply(1:length(ssodm),function(i){
      (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
    }))
  } else {
    print('Undefined method')
  }
  p.values
}

#' This function generates n null transformation matrices from a folder of networks
#'
#' @param path a path to the folder containing the network files (and only the network files)
#' @param keyword a string to be used as a key for distinguishing group 1 from group 2
#' @keywords keywords
#' @export
#' @examples
#' example1
load.null.tms <- function(dir.path, keyword="sub1"){
  perm.filenames <- list.files(dir.path)
  null.networks.list <- split(lapply(file.path(dir.path,perm.filenames), function(x){
    file.to.regnet(x)
  }),grepl(keyword, perm.filenames))

  # Calculate the transition matrices between those "null" networks
  tm.null <- lapply(1:length(null.networks.list[[1]]), function(x){
    transformation.matrix(null.networks.list[[1]][[x]], null.networks.list[[2]][[x]])
  })
  tm.null
}
#' This function reads in panda C output and converts to a reg.net mxn matrix
#'
#' @param path a path to the folder containing the network files (and only the network files)
#' @param keyword a string to be used as a key for distinguishing group 1 from group 2
#' @keywords keywords
#' @export
#' @examples
#' example1
file.to.regnet <- function(file.name){
  library(reshape2)
  library(ggplot2)
  reg.net.melt <- read.table(file.name, header=F)
  reg.net <- dcast(reg.net.melt, V1 ~ V2, value.var='V4')
  rownames(reg.net) <- reg.net[,1]
  reg.net <- as.matrix(reg.net[,-1])
}
