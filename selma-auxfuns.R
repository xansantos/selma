

bootsamplerGN<-function(X){
  
  require( data.table )
  
  DT <- copy( X )
  
  setDT( DT )
  
  setnames( DT, c( "haul", "string", "l", "n_1", "n_2", "q_1", "q_2") )
  
  haul_ev <- DT[, unique( haul ) ] 
  
  haul_bv <- sample(  haul_ev, length(  haul_ev ),replace = T ) # vector with resampled hauls id
  
  DT_bbb <- lapply( seq_along( haul_bv ), 
                    
                    FUN = function( i ){
    
                    haul_b <- DT[ haul == haul_bv[ i ] ]
                                              
                    haul_string_b <- split( haul_b, haul_b$string )
    
                    string_ev <- names( haul_string_b )
                                              
                    string_bv<-if( length( string_ev )==1 ) {
                      
                                                            string_ev
                      
                                                            }else {
                      
                                                            sample(string_ev, length(string_ev), replace=TRUE)
                                                            
                                                                }
  
                    haul_string_bb<-haul_string_b[string_bv]
    
    
    lapply(haul_string_bb, function(x){
      
      x[,':=' (haul=i,
               
               n_1=if(sum(n_1)==0) 0 else {as.integer(as.vector(rmultinom(1, sum(n_1), prob = n_1)))}, 
               
               n_2=if(sum(n_2)==0) 0 else {as.integer(as.vector(rmultinom(1, sum(n_2), prob = n_2)))}
               
               )
        ]
      
                                      })
    
    x<-rbindlist(haul_string_bb)
    
    x[, string:=NULL]
 
    x
    
                                      }
              )
  
  out_dt <- rbindlist(DT_bbb)
  
  setDF(out_dt)
  
  return(out_dt)
  
}



f_bci<-function(x,q){
  
  xm<-x[1]
  
  xci<-quantile(x,probs=c((1-q/100)/2,1-((1-q/100)/2)),na.rm = TRUE)
  
  res<-c(xci[1],Mean.Val=xm,xci[2])
  
  return(res)
  
}
    

pooldata<-function(DT, NA_report=FALSE){
  
  X <- copy(DT)
  
  setDF(X)
  
  n.cols <- grep("^n_", names(X))
  
  q.cols <- grep("^q_", names(X))
  
  n.names <- names( X )[ n.cols ]
  
  N.names <- toupper( n.names )

  q.names <- names( X )[ q.cols ]

  
  ##handle nas
  na.report <- which(is.na(X), arr.ind = TRUE)
  if (length(na.report) > 0) {
    nacols <- unique(na.report[, 2])
    for (i in nacols) {
      if (i %in% q.cols)
        X[na.report[na.report[, 2] == i, 1], i] <- 1
    }
    
    if (NA_report) {
      if (all(nacols %in% q.cols)) {
        cat("Sampling ratio NAs found in rows:", na.report[, 1], "\n")
      } else {
        stop("NA found in count variables")
      }
    }
  } else if (NA_report) {
    cat("No NAs detected\n")
  }
  
  tmp_1<-X[,n.cols, drop = FALSE]
  
  tmp_1[, N.names] <- X[, n.cols, drop = FALSE] / X[, q.cols, drop = FALSE]
  
  
  tmp_2<-apply(tmp_1,2,sum)
  
  tmp_q<-tmp_2[grep("^n_",names(tmp_2))]/tmp_2[grep("^N_",names(tmp_2))]
  
  
  names(tmp_q)<-q.names
  
  pooled<-aggregate(X[, n.cols, drop = FALSE], by=list(X$l), FUN=sum)
  
  for (qn in q.names) pooled[[qn]] <- tmp_q[qn]
  
  names(pooled)<-c("l",n.names,q.names)
  
  pooled
  
}


