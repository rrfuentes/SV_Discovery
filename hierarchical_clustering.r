Data = read.table("E:/SV_Paper/SV_stat/Merging/NB_DUP_mergesam_intersect_uniq.txt")

names(Data) = c("CHR", "START", "END", "CLUST","LEN")
library(data.table)
library(ggplot2)
library('Rcpp')
library('inline')
library(fastcluster)
library(igraph)
library(RColorBrewer)

setDT(Data)

Data[, .N, by=.(CLUST)][order(N)]

outfile="E:/SV_Paper/SV_stat/Merging/NB_DUP_cluster.txt"

if(file.exists(outfile)) file.remove(outfile)

cat("CHR\tSTART\tEND\tLEN\tINI_CLUST\tNEW_CL\n", file=outfile)

cppFunction("
   NumericVector distint_cpp(NumericMatrix x) { 
            // START END LEN
            int n = x.nrow(); // number of rows
            int L1,L2,r,l;
            NumericVector out(n*(n-1)/2); // result numeric vector
            int k = 0;
            for (int i=0; i<n-1; ++i) {
              for (int j=i+1; j<n; ++j) {  //compute only half of the matrix
                L1 = x(i,2);
                L2 = x(j,2);
                r = x(i,0) - x(j,0); //non-intersecting region (start)
                l = x(i,1) - x(j,1); //non-intersectiong region (end)
                if(l<0) l = -l;
                if(r<0) r = -r;
                out[k++] = (double)(r+l)/(L1 + L2); //get distance between two intervals
              }
            }
            return out;
    }
")

cppFunction("
 SEXP buildg_cpp(NumericMatrix x) { 
            // START END LEN
            int n = x.nrow(); // number of rows
            int L1,L2,r,l;
            double dist;
            std::map< int,  std::vector<int> > vec;
            for (int i=0; i<n-1; ++i) {
              for (int j=i+1; j<n; ++j) {
                L1 = x(i,2);
                L2 = x(j,2);
                r = x(i,0) - x(j,0); //non-intersecting region (start)
                l = x(i,1) - x(j,1); //non-intersectiong region (end)
                if(l<0) l = -l;
                if(r<0) r = -r;
                dist = double(r+l)/(L1 + L2); 
                //DETERMINES interval(vertices) that are connected
                if(((r+l)<=10) && (dist < 0.3)){ //allows grouping of shorts events with 70% RO and <11bp distance
                  vec[i+1].push_back(j+1);
                }else if(dist < 0.1){
                  vec[i+1].push_back(j+1);
                }
              }
            }
            return wrap(vec);
      }
")

distint <- function(interval) {
  stopifnot(is.numeric(interval), is.matrix(interval), ncol(interval) == 3)
  
  res <- distint_cpp(interval) # use Rcpp to calculate the distances
  
  # return the result similar to the one of dist()
  structure(res, class='dist', Size=nrow(interval), Diag=FALSE, Upper=FALSE)
}

build_gr <- function(interval) {
  
  res <- buildg_cpp(interval) # use Rcpp to calculate the distances
  res
}

general_hist <- numeric()

for(clust_id in unique(Data$CLUST)){ 

  #clust_id = 10 #49677 2426 1794
  
  cat( file=stderr(),  "Processing vcluster ", clust_id, "\n")
  cl_data = Data[ Data$CLUST == clust_id ,] 

  N = nrow(cl_data)  # number of intervals in the cluster

  if(N>1){
  
    intervals =  as.matrix( cl_data[,c("START","END","LEN"), wi=F])
    
    ###CONNECTED COMPONENTS CLUSTERING###
  
    edges = build_gr(intervals) #get edges based on 0.90 min.RO (Reciprocal Overlap) or 0.1 max boundary error 
    L = sapply(edges,length)
    v1 = rep(as.integer(names(edges)), times=L)
    v2 = unlist(edges)
    edge_list = matrix(c(v1,v2),ncol=2)
    if(length(v2)==0){ #no connected vertices or multiple intervals without 90% RO
      final_cl = as.data.frame( cl_data[,c("START","END","LEN"), wi=F])
      final_cl[ , "NEW_CL"] = 1
      chrom = cl_data$CHR[1] #get chrom of the cluster
      final_cl[ , "CHR" ] = chrom
      final_cl[ , "COMP_ID" ] = seq(1,nrow(cl_data))
      write.table(final_cl[ , c("CHR","START","END","LEN","COMP_ID","NEW_CL")], outfile, append=T,col.names = FALSE, row.names = FALSE,quote = FALSE)
      next
    }
    
    G = graph.edgelist( edge_list, directed = F) #should be undirected since RO is used
    
    
    comp = components(G) #get connected components
    #plot(G, vertex.color=comp$membership, vertex.label=NA)
    #barplot(comp$csize)
    
    #PLOT WITHOUT selection
    #plot(intervals[,1],1:nrow(intervals), type="n")
    #plot( x = c( min(intervals[ ,1]), max(intervals[ ,2])), y = c(1, nrow(intervals)), type="n", xlab="Position", ylab="Dummy index")
    #segments(intervals[,1], 1:nrow(intervals), intervals[,2], 1:nrow(intervals), col = comp$membership)
    
    #PLOT WITH selection
    #selected_components = which( comp$csize > 2 )  #filter based on member count
    #to_plot = which( comp$membership %in% selected_components) #get index of selected component based on component list
    #to_plot = which( comp$membership >= 3179 & comp$membership <= 3179)
    #plot( x = c( min(intervals[ to_plot,1]), max(intervals[ to_plot,2])), y = c(1, length(to_plot)), type="n", xlab="Position", ylab="Dummy index")
    #segments(intervals[to_plot,1], 1:length(to_plot), intervals[to_plot,2], 1:length(to_plot), col = comp$membership[to_plot])
    
    tmp_intervals = cbind.data.frame(intervals, id = 1:nrow(intervals)) #add ID to intervals
    memb = data.frame( id = 1:length(comp$membership), cl_id = comp$membership) #make a data frame for with col:interval idx, cl_id
    #dim(tmp_intervals) and length(comp$membership) may not be equal if some node did not intersect to at least 1 interval
    
    #Merge tmp_intervals and cluster IDs
    intervals_cl = merge( tmp_intervals, memb, by="id", all.x=T)
  
    M = comp$no # number of components
    #intervals_cl[is.na(intervals_cl$cl_id),] #show rows with 'NA' or no cluster/component
    
    #merging intervals & comp(components) may miss some events (not 1event-to-1comp_member); some with cl_id='NA'
    if(nrow(intervals) != length(comp$membership)) intervals_cl[is.na(intervals_cl$cl_id), "cl_id"] = seq(M+1, M+sum(is.na(intervals_cl$cl_id))) #assign numbers to intervals w/o cl_id
    rm(tmp_intervals)
    
    #cmpnt = 236
    
    spreadhist = vector(mode='double', max(intervals_cl$cl_id)) #include the clusters with only 1 interval
    for(cmpnt in  1:comp$no){  
        #SELECTION of component that will require further clustering    
        selected_inter = which( intervals_cl$cl_id == cmpnt)
        cl_elem = as.matrix( intervals_cl[selected_inter,c("START","END","LEN"),]) #select rows with cl_id=cmpnt
        leftmost_bkpt=min(cl_elem[,1]) #get start bkpt of leftmost interval
        rightmost_bkpt=max(cl_elem[,2]) #get end bkpt of rightmost interval
        spread = (rightmost_bkpt - leftmost_bkpt)/ min(cl_elem[,3]) #length of region covered by the component
        spreadhist[cmpnt] = spread
        if(spread > 1.5 && comp$csize[cmpnt]>8){
             #length covered by the cluster over the length of shortest interval
             #0.9^8 = 0.43 or 43% min overlap between fastest intervals in the component
             
             ###HIERARCHICAL CLUSTERING###
            d2 = distint(cl_elem)
            H = hclust( as.dist(d2), method="complete")  
            a = cutree(H, h= 0.1) #set threshold for cluster
            #plot(H)
            #table(a)
            cl_info = data.frame( row_num =H$order, new_cl = a) #create info table
          
            final_cl = data.frame(cl_elem)
            final_cl[ cl_info$row_num, "TREE_ORDER"] = 1:comp$csize[cmpnt] 
            final_cl[ , "NEW_CL"] = as.integer( cl_info$new_cl ) #Add Column
          
        }else{ 
          final_cl = data.frame(cl_elem)
          final_cl[ ,"NEW_CL"] = 1
        }
        
        chrom = cl_data$CHR[1] #get chrom of the cluster
        final_cl[ , "CHR" ] = chrom
        final_cl[ , "COMP_ID" ] = paste(clust_id,cmpnt,sep="_")
        
        #cols = colorRampPalette(brewer.pal(11, "Set3"))
        #G = ggplot(final_cl, aes(x=START, xend=END, y=c(1:nrow(final_cl)), yend =c(1:nrow(final_cl)), col=as.character(NEW_CL) )) + geom_segment()
        #G = G + scale_color_manual(values=cols(26)) 
        #G = G + labs(title = paste0("Component ", cmpnt))
        #G
        
        #BREAKPOINTS
        #ggplot(final_cl, aes(x=START)) + geom_density() #plot based on START
        #ggplot(final_cl, aes(x=END)) + geom_density() #plot based on END
        #ggplot(final_cl, aes(x=LEN)) + geom_density() #plot based on LEN
        
        write.table(final_cl[ , c("CHR","START","END","LEN","COMP_ID","NEW_CL")], outfile, append=T,col.names = FALSE, row.names = FALSE,quote = FALSE)
    }
    general_hist = append(general_hist, spreadhist)
    
    if(nrow(intervals) != length(comp$membership)){
      offset = M+1
      for(cmpnt in  offset:max(intervals_cl$cl_id)){
        spreadhist[cmpnt] = 1
        final_cl = intervals_cl[intervals_cl$cl_id == cmpnt, c("START","END","LEN"),]
        final_cl[ , "NEW_CL"] = 1
        chrom = cl_data$CHR[1] #get chrom of the cluster
        final_cl[ , "CHR" ] = chrom
        final_cl[ , "COMP_ID" ] = paste(clust_id,cmpnt,sep="_")
        write.table(final_cl[ , c("CHR","START","END","LEN","COMP_ID","NEW_CL")], outfile, append=T,col.names = FALSE, row.names = FALSE,quote = FALSE)
      }
    }
  
    
  } else {
    final_cl = as.data.frame( cl_data[,c("START","END","LEN"), wi=F])
    final_cl[ , "NEW_CL"] = 1
    chrom = cl_data$CHR[1] #get chrom of the cluster
    final_cl[ , "CHR" ] = chrom
    final_cl[ , "COMP_ID" ] = clust_id #paste(clust_id,"1",sep="_")
    write.table(final_cl[ , c("CHR","START","END","LEN","COMP_ID","NEW_CL")], outfile, append=T,col.names = FALSE, row.names = FALSE,quote = FALSE)
  }

}

#barplot(general_hist,ylim=c(0,10),ylab="Spread")
#hist(  general_hist,50)
#hist( log(general_hist[general_hist > 1],2) , 100, main=NA,xlab=bquote(log[2]~spread),ylim=c(0,20000))
#sum( general_hist == 1)
#hist( general_hist[ general_hist > 1 ], 50)

#plot(x = comp$csize, y= general_hist[1:comp$no ], log="x",ylab="Spread",xlab="Component Size")
#abline( h = 1.5, lty=2)
#abline( v = 8, lty=2)


