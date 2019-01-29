#####################################################################
#############       get command line options        #################
#####################################################################


library(optparse)

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default="results/counts.txt", help="normalized and filtered counts tabulated file from filter_for_plots.py", metavar="character"),
  make_option(c("-m", "--method_of_clustering"), type="character", default="hierarchical", help="method of clustering chose between: hierarchical, 
                kmean or fuzzy_kmean
                default is hierarchical", metavar="character"),
  make_option(c("-n", "--opt_clust_number"), type="character", default="dynamic", help="determination of optimal number of clusters, 
                If --method_of_clustering=hierarchical the options are: fixed, height or dynamic.
                If--method_of_clustering=kmean or --method_of_clustering=fuzzy_kmean the options
                are: fixed, average_silhouette_width, calinsky_criterion, 
                gap_statistic or affinity_propogation.
                default is dynamic.", metavar="character"),
  make_option(c("-k", "--number_of_clusters"), type="integer", default=5, help="number of clusters to be used. choose an integer of two or more. 
                To be used if --det_hclust_number=fixed or --det_kclust_number=fixed.
                default is 5", metavar="integer"),
  make_option(c("-H", "--height_in_dendrogram"), type="double", default=1.5, help="height in dendrogram to use for determination of cluster number.
               default is 1.5", metavar="double"),
  make_option(c("-q", "--membership_min"), type="double", default=0.2, help="minimal membership value of fuzzy kmean cluster to be included in 
                the cluster. should be a value between 0.0 and 1.0. Only to 
                be used in combination with -m = fuzzy_kmean. Best to first run with a low value,
	        and adjust based on results if needed.
                default is o.2", metavar="double"),
  make_option(c("-r", "--correlation_min"), type="double", default=0.5, help="minimal correlation value to a kmean cluster to be included in 
                the cluster, if too low it will be excluded from all clusters!. should be a value
                between 0.0 and 1.0. Only to be used in combination with -m = kmean or -m = pam. 
                Best to first run with a low value, and adjust based on results if needed.
                default is 0.5", metavar="double"),
  make_option(c("-c", "--colour_of_heatmaps"), type="character", default=c("white","green","green4","violet","purple"), help="colors of the heatmaps. needs to be a list of two or more colors. 
                To see a list of all 657 colour names in R: colors()
                default is the colourblind friendly: 
                c(\"white\",\"green\",\"green4\",\"violet\",\"purple\")", metavar="character"),
  make_option(c("-p", "--plots_output_file"), type="character", default="results/plots.pdf", help="name of multipage pdf file where to output all the plots.", metavar="character"),
  make_option(c("-o", "--clusters_output_file"), type="character", default="results/plots.pdf", help="name of multipage pdf file where to output all the plots.", metavar="character")
  )
# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list , add_help_option = TRUE, description = "\nThis script produces multiple plots and clusters from a RNAseq
output from the filter_for_plots.py. 
plots include dendrograms of the samples and the genes, a heatmap, an elbow plot,
a heatmap with clusters specicied by colorbar, and finally plots and heatmaps of the
different clusters.
Because the 'correct' method for clustering RNAseq data is a matter of perspective; 
it is the one that allows the researcher to make the most out of her data. Gene expression
data is also full of noise which can make clustering tricky when using algorithms optimised
for chunkier data. With that in mind it's good to try several methods and compare them.
So for that multiple ways of clusting can be choosen also the number of clusters can be manually
choosen or calculated by a number of given methods. When making use of the snakemake pipeline,
it is essential to remove or rename the output file (plots.pdf), change the arguments you want to 
change in the config.yaml, and rerun the pipeline with 'snakemake --use-conda' (to be sure only
the last script reruns, start of with 'snakemake -np').", epilogue = "Pease feel free to mail me, bliek@uva.nl, if you would like more info or have suggestions for optimization.\n\n");
opt = parse_args(opt_parser);

i <- opt$input_file
m <- opt$method_of_clustering
n <- opt$opt_clust_number
k <- opt$number_of_clusters
h <- opt$height_in_dendrogram
c <- opt$colour_of_heatmaps
q <- opt$membership_min
r <- opt$correlation_min
p <- opt$plots_output_file
o <- opt$clusters_output_file

######################################################################
##############          import needed packages         ###############
######################################################################


library(gplots)
library(dendextend)
library(dynamicTreeCut)
#library("pheatmap")
library(cluster)
library(vegan)
library(apcluster)
#library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(colorRamps)
library(e1071)
library(tidyr)
library(ggplot2)

######################################################################
###############              read the data           #################
######################################################################

data <- read.delim(i, header=T, row.names="genes")
z <- as.matrix((data))

# scaling the data to avoid clustering based on expression level
scaledata <- t(scale(t(z))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]

# open pdf for the plots
pdf(file=p)

#mean/variance calculations and plot
z_var <- apply(z, 1, var)
z_mean <- apply(z, 1, mean)
plot(log2(z_mean), log2(z_var), pch='.')


##############################################################
############    hierarchical dendrogram      #################
##############################################################

# dendrogram of the samples
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
TreeC = as.dendrogram(hc, method="average")
plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")

# dendrogram of the genes
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
TreeR = as.dendrogram(hr, method="average")
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

#############################################################
##############        create Heatmap       ##################
#############################################################

# colors of the heatmap
my_palette <- colorRampPalette(c)(100)
#(c("magenta", "black", "green"))(n = 299)

# plot the heatmap
heatmap.2(z,
          Rowv=as.dendrogram(hr), 
          Colv=NA,
          col=my_palette,
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          dendrogram = "row",
          main = "Heatmap.2",
          trace = "none")


#################################################################
#######   Elbow plot ( the sum of squared error (SSE))    #######
#################################################################

wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


##############################################################
##############        Do the clustering        ###############
##############################################################



if(m == "kmean" || m == "fuzzy_kmean" || m == "pam"){
  if(n == "average_silhouette_width"){
    # method 1: average silhouette width
    sil <- rep(0, 20)
    #repeat k-means for 1:20 and extract silhouette:
    for(i in 2:20){
      k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
      ss <- silhouette(k1to20$cluster, dist(scaledata))
      sil[i] <- mean(ss[, 3])
    }
    # Plot the  average silhouette width
    plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
    abline(v = which.max(sil), lty = 2)
    
    k <- which.max(sil)
  }
  else if(n == "calinsky_criterion"){
    # method 2: Calinsky criterion
    fit <- cascadeKM(scaledata, 1, 20, iter = 100)
    plot(fit, sortg = TRUE, grpmts.plot = TRUE)
    calinski.best <- as.numeric(which.max(fit$results[2,]))
    cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
    
    k <- calinski.best
  }
  else if(n == "gap_statistic"){
    # method 3: Gap statistic
    set.seed(13)
    gap <- clusGap(scaledata, kmeans, 20, B = 100, verbose = interactive())
    plot(gap, main = "Gap statistic")
    abline(v=which.max(gap$Tab[,3]), lty = 2)
    
    k <- which.max(gap$Tab[,3])
  }
  else if(n == "affinity_propogation"){
    # method 4: Affinity propogation
    d.apclus <- apcluster(negDistMat(r=2), scaledata)
    cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
    #uncomment the next line for the heatmap, it takes a long time to run
    #heatmap(d.apclus,cexRow=0, cexCol=0)
    
    k <- length(d.apclus@clusters)
  }
  else if(n == "fixed"){
    k <- k
  }else{
    print("your choosen method of determining the optimal number of clusters is unclear.")
    print("When using --method_of_clustering=kmean or --method_of_clustering=fuzzy_kmean ")
    print("the options are:fixed, average_silhouette_width, calinsky_criterion,")
    print("gap_statistic or affinity_propogation.")
  }
  print(paste0("Dataset is clustered in ",k, " clusters "))
  if(m == "kmean"){
    set.seed(20)
    clusts <- kmeans(scaledata, centers=k, nstart = 1000, iter.max = 20)
    clusts <- clusts$cluster
  }
  else if(m == "fuzzy_kmean"){
    mestimate<- function(df){
      N <-  dim(df)[[1]]
      D <- dim(df)[[2]]
      m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
      return(m.sj)
    }
    mes <- mestimate(scaledata)
    fcm_results <- cmeans(scaledata,centers=k,m=mes)
    fcm_plotting_df <- data.frame(scaledata)
    fcm_plotting_df$gene <- row.names(fcm_plotting_df)
    fcm_plotting_df$cluster <- fcm_results$cluster
    fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster),function(row){
      clust <- fcm_plotting_df$cluster[row]
      fcm_results$membership[row,clust]
    })
    clusts <- fcm_plotting_df$cluster
    names(clusts) <- fcm_plotting_df$gene
    # get the core data
    fcm_centroids <- fcm_results$centers
    fcm_centroids_df <- data.frame(fcm_centroids)
    fcm_centroids_df$cluster <- row.names(fcm_centroids_df)
    centroids_long <- tidyr::gather(fcm_centroids_df,"sample",'value', 1:ncol(scaledata))
    #start with the input data
    fcm_plotting_df <- data.frame(scaledata)
    #add genes
    fcm_plotting_df$gene <- row.names(fcm_plotting_df)
    #bind cluster assinment
    fcm_plotting_df$cluster <- fcm_results$cluster
    #fetch the membership for each gene/top scoring cluster
    fcm_plotting_df$membership <- sapply(1:length(fcm_plotting_df$cluster),function(row){
      clust <- fcm_plotting_df$cluster[row]
      fcm_results$membership[row,clust]
    })
    # filter out genes don't really belong to any of the clusters
    selection <- fcm_plotting_df %>% filter(membership > q)
  }
  else if(m == "pam"){
    pam.Pclusts <- pam(scaledata, k=k)
    clusts <- pam.Pclusts$'clustering'
  }
}else if(opt$method_of_clustering == "hierarchical"){
  if(n == "dynamic"){
    clusts <- cutreeDynamic(hr, distM = as.matrix(as.dist(1-cor(t(scaledata)))), method = "hybrid")
    names(clusts) <- rownames(scaledata)
    clust <- as.data.frame(clusts)
    colnames(clust) <- "cluster"
    k = length(unique(clusts))
    print(paste0("Dataset is clustered in ",k, " clusters "))
  }
  else if(n == "fixed"){
    clusts = cutree(hr, k=k)
    k = length(unique(clusts))
    print(paste0("Dataset is clustered in ",k, " clusters "))
  }
  else if(n == "height"){
    clusts = cutree(hr, h=h)
    k = length(unique(clusts))
    print(paste0("Dataset is clustered in ",k, " clusters "))
  }
  else{
    print("your choosen method of determining the optimal number of clusters is unclear.\
          When using --method_of_clustering=hierarchical the options are: fixed, height or dynamic.")
  } 
  }else{
    print("your choosen method of clustering is unclear. 
          It should be one of the following: hierarchical, kmean or fuzzy_kmean")
}


###################################################################
############        Heatmap with clusterbar          ##############
###################################################################

mycolhc <- rainbow(length(unique(clusts)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(clusts)]

hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete")
	  
heatmap.2(z,
          Rowv=as.dendrogram(hr), 
          Colv=NA,
          col=my_palette,
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          dendrogram = "row",
          main = paste0("Heatmap.2 with colour bar indicating the clusters (",k,")"),
          trace = "none",
          RowSideColors=mycolhc,
          key = FALSE)


###############################################################################
########   calculate the cluster 'cores' aka centroids and plot them   ########
###############################################################################

# get the cluster centroids (the average profile of expression for each of the clusters)
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
clusUniq <- unique(clusts)
kClustcentroids <- sapply(levels(factor(clusts)), clust.centroid, scaledata, clusts)
head(kClustcentroids)

Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

# In case of use of pam, the cluster centroids wil be the medoids (the expression profiles of the genes that fit the clusters best)
if (m == "pam"){
  pClustcentroids <- t(pam.Pclusts$"medoids")
  Kmolten <- melt(pClustcentroids)
  colnames(Kmolten) <- c('sample','cluster','value')
} 
	  
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression of the samples",color = "Cluster")
p1

# in case of pam, change names of the clusters to numbers in stead of the medoids gene names
if (m == "pam"){
  colnames(pClustcentroids) <- clusUniq
  Kmolten <- melt(pClustcentroids)
  colnames(Kmolten) <- c('sample','cluster','value')
}
###########################################################################
#############      heatmaps and plots of the clusters      ################
###########################################################################

# open list to collect the correlation values in case of kmean, fuzzy_kmean and pam
cors <- c()

# loop though the clusters and create expression profile plots and heatmap for each of them
if(m == "fuzzy_kmean"){
  for(i in 1:k){
    print(i)
    cluster_plot_df <- dplyr::filter(selection, cluster == i) %>%
      dplyr::select(.,1:ncol(scaledata),membership,gene) %>%
      tidyr::gather(.,"sample",'value',1:ncol(scaledata)) 
    #order the dataframe by score
    cluster_plot_df <- cluster_plot_df[order(cluster_plot_df$membership),]
    #set the order by setting the factors using forcats
    cluster_plot_df$gene = forcats::fct_inorder(cluster_plot_df$gene)
    
    #subset the cores by cluster
    core <- dplyr::filter(centroids_long, cluster == i)
    
    g1 <-  ggplot(cluster_plot_df, aes(x=sample,y=value)) + 
      geom_line(aes(colour=membership, group=gene)) +
      scale_colour_gradientn(colours=c('blue1','red2')) +
      #this adds the core 
      geom_line(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
      geom_point(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
      xlab("Sample") +
      ylab("Expression") +
      labs(title= paste0("Cluster ",i," Expression of the Samples"),color = "Score")
    print(g1)
    
    group <- selection %>% filter(selection$cluster==i)
    rownames(group) <- group$gene
    group <- as.matrix(group[,1:ncol(scaledata)])
    if(nrow(group) > 1){
      hrc <- hclust(as.dist(1-cor(t(group), method="pearson")), method="complete")
      #pheatmap(group, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = F)  #cellheight = 8
      heatmap.2(group,
                Rowv=as.dendrogram(hrc), 
                Colv=NA,
                col= magenta2green(100),
                labRow=rownames(group),
                scale="row",
                margins = c(7, 7),
                cexCol = 0.7,
                dendrogram = "row",
                main = paste0("Cluster ",i, " consisting ", nrow(group), " genes"),
                trace = "none",key = F)
    }
  }
}else if(m == "kmean" || m == "hierarchical" || m == "pam"){
  for(i in clusUniq){
    K2 <- (scaledata[clusts==i,])
    if(nrow(K2) > 1){
      #calculate the correlation with the core
      core <- Kmolten[Kmolten$cluster==i,]
      corscore <- function(x){cor(x,core$value)}
      score <- apply(K2, 1, corscore)
      # add the correlations to the list
      cors <- c(cors, score)
      #get the data frame into long format for plotting
      K2molten <- melt(K2)
      colnames(K2molten) <- c('gene','sample','value')
      #add the score
      K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
      colnames(K2molten) <- c('gene','sample','value','score')
      #order the dataframe by score
      #to do this first create an ordering factor
      K2molten$order_factor <- 1:length(K2molten$gene)
      #order the dataframe by score
      K2molten <- K2molten[order(K2molten$score),]
      #set the order by setting the factors
      K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
      
      # Everything on the same plot
      if(m == "kmean" || m =="pam"){
        p2 <- ggplot(K2molten, aes(x=sample,y=value)) + 
          geom_line(aes(colour=score, group=gene)) +
          scale_colour_gradientn(colours=c('blue1','red2')) +
          #this adds the core 
          geom_line(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
          geom_point(data=core, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
          xlab("Samples") +
          ylab("Expression") +
          labs(title= paste0("Cluster ",i, " consisting ", nrow(K2), " genes"),color = "Score")
      }
      else if(m == "hierarchical"){
        p2 <- ggplot(K2molten, aes(x=sample,y=value)) + 
          geom_line(color="grey", aes(color="grey", group=gene)) +
          #this adds the core 
          geom_line(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
          geom_point(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
          xlab("Samples") +
          ylab("Expression") +
          labs(title= paste0("Cluster ",i, " consisting ", nrow(K2), " genes"),color = "Score")
      }
      print(p2)
      hrc <- hclust(as.dist(1-cor(t(K2), method="pearson")), method="complete")
      #pheatmap(group, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, cluster_cols = F)  #cellheight = 8
      heatmap.2(K2,
                Rowv=as.dendrogram(hrc), 
                Colv=NA,
                col= my_palette,
                labRow=rownames(K2),
                scale="row",
                margins = c(7, 7),
                cexCol = 0.7,
                dendrogram = "row",
                main = paste0("Cluster ",i, " consisting ", nrow(K2), " genes"),
                trace = "none",key = F)
    }
  }
}
dev.off() 

#########################################################################################
#############       create output file containing the clusters          #################
#########################################################################################

# create data frame containing cluster and if available correlation values or memberships
if(m == "fuzzy_kmean"){
  # fuzzy_kmeans
  lijst <- data.frame(selection$gene, selection$cluster, selection$membership )
  colnames(lijst) <- c("gene", "clusters", "membership")
}else if(m == "kmean" || m == "pam"){
  # kmean and pam
  lijst <- merge(as.data.frame(clusts), as.data.frame(cors), by="row.names", sort=FALSE)
  colnames(lijst) <- c("gene", "clusters", "correlation")
  lijst <- lijst %>% filter(correlation > r)
}else{
  # hierarchical
  lijst <- data.frame(rownames(as.data.frame(clusts)), clusts)
  colnames(lijst) <- c("gene", "clusters")
}

# Write the dataframe as a tab-delimited file
write.table(lijst, file=o, sep = "\t",quote=F,row.names=F)
