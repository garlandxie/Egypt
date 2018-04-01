############### Phylogenetic Eigenvector Regression - Egyptian Seed Traits ############################
# This code is used to calculate phylogenetic eigenvector regression for seed traits in Egyptian plants
# Code developed by Garland Xie
#######################################################################################################

# Load libraries ----------------------------------------------------------------------------------------------
library(ape)
library(geiger)
library(here)
library(picante)
library(splancs)

# Load datasets -----------------------------------------------------------------------------------------------


# Custom functions --------------------------------------------------------------------------------------------

# PVR: Phylogenetic Eigenvectors Regression

# arguments
# @x: An object of class PVR (created by the PVRdecomp function) or class PSR (requiered by the "PSR" method).
# @phy:An object of class phylo that contains an ultrametric phylogeny.
# @trait:A vector, data frame or matrix that contains traits sets (for data frames and matrices, each column must represent a trait set).
# @envVar:A vector, data frame or matrix that contains environmental variables. Used to estimates the variation of a trait set that is explained by phylogeny and by environment.
# @method: Character string. A name for the eigenvectors selection method. It can be "moran", "stepwise", "psr" or "sequential".
# @weights: Weighting matrix based on phylogentic distances used in the "moran" method. If no weights matrix is provided, weights will be set to max(D) - Dij, where D is the phylogenetic distance matrix. 
# @scaled: Logical. Should the phylogenetic distances be scaled into the range of 0 to 1. Default is FALSE.
# @significance: Logical. Should the eigenvectors selected by the "moran" method be selected by the significance of residuals autocorrelation. If FALSE the eigenvectors will be selected by Moran's I values.
# @alterrnative: The alternative hypothesis used to compute the significance level when selecting eigenvectors by the "moran" method.
# @MI.treshold: Minimum residuals Moran's I value used to select eigenvectors when significance is FALSE.
# @sig.treshold: The significance treshold used to select eigenvectors by the "moran" method.
# @psr.treshold: The minimum acumulate R2 gain treshold used to select eigenvectors by the "PSR" method.
# @accvalue.treshold: Relative accumulated eigenvalue treshold use to select the eigenvectors by the "sequential" method.
# @dots: Parameters passed to the stepwise regression used in the "AIC" method

# return
# @value: A PVR class object

PVR <- function(x, phy = NULL, trait = NULL, envVar = NULL, method = "moran", 
                weights = NULL, scaled = FALSE, significance = TRUE, alternative = "two.sided",
                sig.treshold = 0.05, MI.treshold = 0.05, psr.treshold = 0.01, accvalue.treshold = 0.9, ...){
  
  if(method == "moran" | method == "Moran" | method == "MoransI" | method == "Moransi"){
    
    pvr <- x@Eigen$vectors
    Naxis <- ncol(pvr)
    if(is.null(weights)){
      
      W <- max(x@phyDist) - x@phyDist	
    } else {
      W <- weights 
    }
    diag(W) <- 0
    
    N <- .poss(pvr)
    c1 <- 1
    c2 <- 1
    
    tmpRes <- matrix(nrow = Naxis, ncol = 2)
    selected <- matrix(nrow = nrow(pvr), ncol = Naxis)
    reg <- data.frame(Vectors = 0, MoransI = 0, p = 0)
    tmpTrait <- trait
    
    for(i in 1:N){
      
      pvr <- as.matrix(pvr)
      tmpLM <- lm(tmpTrait ~ pvr[ ,c1])
      MI <- Moran.I(tmpLM$residuals, weight = W, scaled = scaled, alternative = alternative)
      tmpRes[c1, 1] <- round(MI$p.value, 5)
      tmpRes[c1, 2] <- MI$observed
      
      if(c1 == ncol(pvr)){
        
        selAxis <- which(tmpRes[ ,2] == max(tmpRes[ ,2]))
        tmpLM <- lm(tmpTrait ~ pvr[ ,selAxis[1]])
        tmpTrait <- tmpLM$residuals
        selected[ ,c2] <- pvr[ ,selAxis[1]]
        reg[c2,1] <- colnames(pvr)[selAxis[1]]
        reg[c2,2] <- tmpRes[selAxis[1], 2]
        reg[c2,3] <- tmpRes[selAxis[1], 1]
        tmpColNames <- colnames(pvr)
        pvr <- pvr[ ,-selAxis[1]]
        tmpColNames <- tmpColNames[-selAxis[1]]
        pvr <- as.matrix(pvr)
        colnames(pvr) <- tmpColNames 
        if(significance){
          
          if(tmpRes[selAxis[1], 1] > sig.treshold | c2 == Naxis){
            
            break
          }	
        } else {
          
          if(tmpRes[selAxis[1], 2] > MI.treshold | c2 == Naxis){
            
            break
          }
        }
        
        
        c2 <- c2 + 1
        c1 <- 0
        
        tmpRes <- matrix(nrow = ncol(pvr), ncol = 2)
      }
      
      c1 <- c1 + 1
    }
    
    selected <- as.matrix(selected[ ,1:nrow(reg)])
    colnames(selected) <- reg$Vectors
    
    selection <- list(Method = "Moran's I", Vectors = selected, Id = reg)
    x@Selection <- selection
  } 
  
  if(method == "AIC" | method == "aic" | method == "stepwise"){
    
    pvr <- x@Eigen$vectors
    Naxis <- ncol(pvr)
    
    iniLM <- glm(trait ~ c1, data = as.data.frame(pvr))
    
    scopeForm <- trait ~ c1
    for(i in 2:ncol(x@Eigen$vectors)){
      
      charForm <- paste("~ . + ", colnames(x@Eigen$vectors)[i], sep ="") 
      upForm <- as.formula(charForm)
      scopeForm <- update.formula(scopeForm, upForm)
    }
    
    stepLM <- stepAIC(iniLM, scope = scopeForm, data = as.data.frame(pvr), ...)
    sels <- length(attr(terms(stepLM), "variables"))
    reg <- data.frame(Vectors = rep(0, (sels - 2)))
    reg$Vectors <- as.character(attr(terms(stepLM), "variables")[3:sels])
    
    selected <- data.frame(x = rep(0,nrow(pvr)))
    for(i in 1:(sels - 2)){
      
      selected[,i] <- pvr[ ,which(colnames(pvr) == reg$Vectors[i])]
    }
    colnames(selected) <- reg$Vectors
    
    selection <- list(Method = "stepwise AIC", Vectors = selected, Id = reg, AIC = stepLM$aic, Formula = stepLM$formula)
    x@Selection <- selection
  }
  
  if(method == "psr" | method == "PSR"){
    
    psr <- x@PSR
    pvr <- x@Eigen$vectors
    Naxis <- ncol(pvr)
    id <- "c1"
    reg <- data.frame(Vectors = 0, cumulR2 = psr$r.squared[1], deltaR2 = psr$r.squared[1])
    selected <- data.frame(x = rep(NA, nrow(pvr)))
    selected[ ,1] <- pvr[ ,Naxis]
    c1 <- 2
    for(i in 2:(Naxis-1)){
      
      deltaR2 <- psr$r.squared[i] - psr$r.squared[(i - 1)]
      if(deltaR2 >= psr.treshold){
        
        selected[ ,c1] <- pvr[ ,i]
        id <- c(id, colnames(pvr)[i])
        reg[c1, 2] <- psr$r.squared[i]
        reg[c1, 3] <- deltaR2 
        c1 <- c1 + 1
      }
    }
    
    reg[ ,1] <- id
    colnames(selected) <- reg$Vectors
    
    selection <- list(Method = "PSR", Vectors = selected, Id = reg)
    x@Selection <- selection
  }
  
  if(method == "sequential"){
    
    relAccVal <- numeric(Naxis)
    pvr <- x@Eigen
    relVal <- pvr$values/sum(pvr$values)
    for(i in 1:Naxis){
      relAccVal[i] <- sum(relVal[1:i])
    }
    pvr <- x@Eigen@vectors
    selected <- pvr[,which(relAccVal <= accvalue.treshold)]
    reg <- colnames(selected)
    selection <- list(Method = "Sequential", Vectors = selected, Id = reg)
    x@Selection <- selection
  }
  
  if(is.null(envVar)){
    pvrOLS <- lm(trait ~ selection$Vectors)
    x@PVR <- list(R2 = summary(pvrOLS)$r.squared, Residuals = pvrOLS$residuals)
  } else{
    
    pvrOLS <- lm(trait ~ selection$Vectors + envVar)
    x@PVR <- list(R2 = summary(pvrOLS)$r.squared, Residuals = pvrOLS$residuals, p = (summary(pvrOLS))$coefficient[2, 4])
  }
  if(!is.null(envVar)){
    
    ABC <- lm(trait ~ selection$Vectors + envVar)
    AB <- lm(trait ~ envVar)
    BC <- lm(trait ~ selection$Vectors)
    D <- 1 - summary(ABC)$r.squared
    B <- summary(AB)$r.squared + summary(BC)$r.squared - summary(ABC)$r.squared
    A <- summary(AB)$r.squared - B
    C <- summary(BC)$r.squared - B
    x@VarPart <- list(a = A, b = B, c = C, d = D)
  }
  
  return(x)
}