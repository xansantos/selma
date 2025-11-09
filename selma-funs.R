# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  SELMA: Functions for estimating smooth relative SELectivity curves
#          using weighted polynomial basis via Model Averaging
#
#  Author: Juan Santos
#  
#  Description:
#  
#  This script contains a collection of functions developed for the analyses
#  presented in Milanelli et al. (in review). These tools are provided solely
#  for transparency and reproducibility purposes and are not intended as a
#  general-purpose R package. 
#  
#  You are free to get the code and adapt it as you wish. 
#
#  For additional details, please refer to the associated README file, or
#  contact me at juan.santos@thuenen.de
#  
#  Last update: 25 October 2025
#  
#  #  Dependencies: 
#  pkgs = data.table, MuMIn, doMC, foreach
#  
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#X=X_ple
#ncores=4
#lengthRange=c(19,45)
#typePredict="response"
#Rank="AICc"
#sampler="bootsamplerGN"
#exportfuns=c("sampler","pooldata","glmAvg")
#type_boot="comp"
#B=100

# core selma function
selma<-function(dt_in,
                lengthRange = NULL,
                report = c( "simple", "full", "boot" ),
                ranking = "AIC",
                typePredict = "response",
                visualReport = T){
  
  ## libraries and aux function:
  require( MuMIn )
  
  require( data.table )  
  
  nseq <- function(X, len = length(X)) seq(min(X, na.rm = TRUE),max(X, na.rm=TRUE), length = len)
  
  ## copy data and transform it to a data.table object:
  
  X <- copy( dt_in ) 
  
  setDT( X )
  
  ## data transformations in place: raised data and  
  X[ ,':=' ( q = q_1 / q_2, N_1 = n_1 / q_1, N_2 = n_2 / q_2 ) ] # calculate raised data
  
  X[ ,':=' ( ptest = N_1 / ( N_1 + N_2 ) , catch = N_1 + N_2 ) ] # empirical catch proportion in test, total catch
  
  ## vectorized aggregated catch data
  catch <- X[[ "catch" ]]
  
  ptest <- X[[ "ptest " ]]
  
  l <- X[[ "l" ]]
  
  ## prediction grid of lengths
  
  if( is.null( lengthRange ) ) lengthRange <- range( l )
  
  l_grid <- seq( lengthRange[ 1 ], lengthRange[ 2 ], 0.5)
  
  ## full model fit
  full_model <- glm( cbind( n_1, n_2 ) ~ l + I( l^2 ) + I( l^3 ) + I( l^4 ) + offset( log( q ) ), data = X, family = binomial )
  
  ## automatic submodels fit based on the full model 
  options( na.action = "na.fail" )
  
  ranking_models <- dredge( full_model, fixed = "offset(log(q))", rank = ranking )
  
  ## list of models for model averaging
  ranking_models_ma <- get.models( ranking_models, subset = delta <= 10 ) # models with delta > 10 excluded
  
  
  ranking_models_n <- length(ranking_models_ma)
  
  ## new data table for average prediction
  model_matrix <- setDT( data.frame( sapply( 1:4, function( i ) I( l_grid^( i ) ) ) ) )
  
  setnames( model_matrix, c( "l", "l^2", "l^3", "l^4" ) )
  
  newdata <- as.data.table( lapply( lapply( model_matrix, mean ), rep, nrow( model_matrix ) ) )
  
  newdata[,':='( l = l_grid, q = 1) ] 
  
  
  if(length( ranking_models_ma ) ==1 ){ # if the list of models consist of only one
    
                                       models_container <- get.models( ranking_models_ma, subset = 1 )[[1]]
    
                                       ccl <-  predict( models_container, newdata, type = typePredict )
    
    
                                       } else{
    
                                              models_averaged <- model.avg(ranking_models_ma)
    
                                              ccl <-  predict( models_averaged, newdata,type = typePredict )
    
                                              }
  
  
  
  if ( report == "full" ) {
    out <- list(
      ranking_models = as.data.frame( ranking_models ),
      models_averaged = if ( length( ranking_models_ma ) > 1 ) {
        models_averaged
      } else {
        NA_character_
      },
      ccl = ccl
    )
    
  } else if  (report == "simple" ) {
    out <- list(
      ccl = ccl
    )
    
  } else if ( report == "boot" ) {
    out <- list(
      typePredict = typePredict,
      l = l,
      p_l = ptest,
      catch = catch,
      L = l_grid,
      ccl = ccl
    )
  }
  
  
  if(visualReport==T & typePredict=="response"){
    
    par( cex = 1.3, cex.axis = 1,cex.main = 1.3 )
    
    plot( ccl ~ l_grid,type = "n", 
          bty = "n", 
          col = 2,
          ylim = c( 0, 1 ),
          xlim = lengthRange,
          ylab = "Catch share in TEST gear",
         xlab = "fish length (cm)"
         )
    
    points( ptest ~ l,
            pch = 21,
            col = "darkblue",
            cex = 3*catch / max( catch ),
            bg="darkgrey"
           )
    
    lines( ccl ~ l_grid,
           type = "l",
           lwd = 3,
           col = 1
           )
    
    abline( h = .5,
           lwd = 2,
           lty = 2,
           col = 2
           )  
    
  }
  
  
  
  return(out)
  
  
  
}

selma_boot<-function(DT_in,
                     lengthRange,
                     typePredict="response",
                     ranking="AICc",
                     sampler,
                     exportfuns=c("sampler","pooldata","selma"),
                     type_boot="semiparametric",
                     ncores=4,
                     B=1000){
                     
  
  require(data.table)
  
  cat("This process is using", ncores, "cores from your computer \n", "be pacient and wait for the results" )
  
  require(doParallel)
  
  
  require(foreach)
  
  
  sampler<-get(sampler)
  
  set.seed(999)
  
  X<-copy(dt_in)
  
  # Create and register a cluster that works on any OS
  cl <- makeCluster(ncores)
  
  registerDoParallel(cl)
  
  
  ## generate boostrap distributions
  r <- foreach(b = 1:B, 
               .export=exportfuns,
               .packages = c("data.table", "MuMIn")) %dopar% {
    
    if( b > 1)  Xb <- sampler( X ) else Xb <- copy(X)
    
    Xo <- setDT( pooldata( Xb ) ) 
    
    tryCatch(suppressWarnings(selma(X = Xo,
                                    lengthRange = lengthRange,
                                     report = "boot",
                                     ranking = ranking,
                                     typePredict = typePredict,
                                     visualReport = F
                                    )
                              ),
             error = function(e) {cat("")})
    
                                                               }
  
  

  
  ## Extract bootstrap distributions 
  b.ccl <- do.call("cbind", lapply( r, function( b ) b$ccl ) )
  
  b.crl <- b.ccl / ( 1 - b.ccl )

  ## extract empirical information
  l <- r[[ 1 ]]$l
  
  L <- r[[ 1 ]]$L
  
  p_l <- r[[ 1 ]]$p_l
  
  n_l <- as.integer(r[[ 1 ]]$catch)
  
  dt_out <- data.table( l = l,
                        p_l = p_l,
                        n_l =  n_l 
                        )
  
  ## percentile confidence intervals
  ccl.ci<-t(apply(b.ccl,1,function(x) f_bci(x,q=95)))
  
  crl.ci<-t(apply(b.crl,1,function(x) f_bci(x,q=95)))
  
  ccl.ci<-if(r[[1]]$typePredict=="link") apply(ccl.ci,2,boot::inv.logit) else ccl.ci
  
 # crl.ci<-if(r[[1]]$typePredict=="link") apply(crl.ci,2,boot::inv.logit) else crl.ci
  
  ccl.ci <- data.frame( l = L, ccl.ci )
  
  crl.ci <- data.frame( l = L,crl.ci )
  
  names( ccl.ci ) <- c( "l", "lci", "mu", "uci")
  
  names( crl.ci ) <- c( "l", "lci", "mu", "uci" )
  
  return(list( B = B, dt_out = dt_out, ccl.ci = ccl.ci, crl.ci = crl.ci ) ) 
                                    }













