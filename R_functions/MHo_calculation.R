#This small routine compute the Morisita-Horn overlap index between two cell types (in this case positive 'p' and lymphocyes 'l')
# This function runs for a set of files/saples stored in path.to.files
#This function returns a dataframe with two columns: sample= sample id as the complete path to the file and mh=Morisita=Horn index.



path.to.files = "/folder/parentalFolder"


.mhWSIdata <- function(path.to.files, coloc=c('p','l')){

  pacman::p_load(parallel, foreach,spaa, sp, IDmining, spatstat)

  if(length(coloc)!= 2)stop('colocalisation argument must include two cell classes')

# Morisita-Horn overlap function
  #x1 and x2 are the abunadnce per sampling unit. Note the length of x1 and x2 must be equal.
  .mh<-function(x1, x2){
    if(length(x1)!=length(x2)) stop('length of x1 and x2 differs')
    mhdf<-data.frame(x1=as.numeric(x1), x2=as.numeric(x2))
    mhdf$x1x2 = mhdf$x1*mhdf$x2
    mhdf$sqx1<-mhdf$x1*mhdf$x1
    mhdf$sqx2<-mhdf$x2*mhdf$x2
    mh<-2*sum(mhdf$x1x2)/((sum(mhdf$sqx1)/(sum(mhdf$x1))^2 + (sum(mhdf$sqx2)/(sum(mhdf$x2))^2))*sum(mhdf$x1)*sum(mhdf$x2))
    return(mh)
  }

  path <- paste0(path.to.files,"/*")
  filespos <- Sys.glob(path)
  mh_df<-data.frame()

  foreach(k = 1:length(filespos)) %do%{ #
    tryCatch({
      print(k)
      datacsv <- read.csv(filespos[k], sep= ",", header=T)
      df<-datacsv[datacsv$class %in% coloc,]
      m <- 2
      scaleQ <- 1:500

      #This computes Morisita dispersion index at different spatial scales
      mh_d <- MINDEX_SP(df[,c('x','y')], scaleQ, m, 2,
                          Wlim_x = c(min(datacsv$x),max(datacsv$x)),
                          Wlim_y = c(min(datacsv$y),max(datacsv$y)))

      #Get the value of distance at which Morisita dispersion is maximim
      deltamaxI <- mh_d[mh_d$m2 == max(mh_d$m2),]$Delta

      #Create point spatial objects
      pp_1<- as.ppp(df[df$class==coloc[1],c('x','y')],owin(xrange = range(df$x), yrange =range(df$y)))
      pp_2<- as.ppp(df[df$class==coloc[2],c('x','y')],owin(xrange = range(df$x), yrange =range(df$y)))

      #Count the number of cells of each type defined in the argument of the function for each quadrant
      count1<-quadratcount(pp_1, (max(df$x)-min(df$x))/sqrt(deltamaxI^2/2),(max(df$y)-min(df$y))/sqrt(deltamaxI^2/2))
      count2<-quadratcount(pp_2, (max(df$x)-min(df$x))/sqrt(deltamaxI^2/2),(max(df$y)-min(df$y))/sqrt(deltamaxI^2/2))

      mh_df<- rbind(mh_df, data.frame(sample=filespos[k],mh=.mh(count1, count2)))

    },
    error= function(e){})
  }
 return(mh_df)
}






