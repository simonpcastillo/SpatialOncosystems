#This function takes as input a three-clolumn dataframe (x, y, class) and compute a cluster detection based on
#user-defined subset of classes, neighbourhooud radio and size. In adition for each cluster.

#The function returns a dataframe of n x m ; with n being the number of clusters detected in the sample. and m equals to 3 + length(neighbourhood.classes).
# The columns display the sample id, the cluster id, the abundnace of each cell class in the cluster,and the area of the cluster.


scanClusters <- function(df=df,neighbourhood.radius = 20, neighbourhood.size = 3, neighbourhood.classes = c("p", "n", "l", "s")) {
  pacman::p_load(fpc, factoextra, dbscan, sp)

  df <-df[df$class %in% neighbourhood.classes,]

  # Compute DBSCAN using fpc package
  db <- fpc::dbscan(df[,c("x", "y")], eps = neighbourhood.radius, MinPts = neighbourhood.size)

  df$cluster <- db$cluster
  g <- unique(df$cluster)

  df.cluster <- data.frame()

  for (i in g) { #g is a vector storing each cluster id
    box.coords <- df[df$cluster==i,c("x", "y")]
    box.hpts <- chull(x=box.coords$x, y=box.coords$y)
    box.hpts <- c(box.hpts, box.hpts[1])
    box.chull.coords <- box.coords[box.hpts,]
    chull.poly <- Polygon(box.chull.coords, hole=F)
    chull.area <- chull.poly@area
    abundances<-c()
    for (class in neighbourhood.classes){
      abx<- nrow(df[df$cluster==i & df$class == class,])
      abundances<- append(abundances, abx)
    }
    names(abundances)<- paste0('n.',neighbourhood.classes)
    abundances <- t(data.frame(abundances))
    df0 <- data.frame(sample.id = sample.id, cluster = i, abundances, chull.area = chull.area)
    df.cluster <- rbind(df.cluster, df0)
  }

  return(df.cluster)
}
