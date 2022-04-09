PLS.PCA.CLASS<-function(DATA,Y,VarName,prop,Method,Mes)
  {
  # Entradas
  # DATA : Matriz de atributos
  # Y    : vector de clases (para dos clases)
  # prop  : definir la cantidad de variables
  # Method : 1 variables.select PLS , 2 PLS.Regression, 3 PCA 
  DATA<-scale(DATA, center = TRUE, scale = TRUE)
  VarName<-colnames(DATA)
  if (Method==1){
    Aux <-pcaSVD(DATA,prop)
    for(i in 1:length(Y)){
      if(Y[i]==-1){Y[i]=2}
    }
    PLSvar<-variable.selection(DATA,Y,nvar=Aux$n)
    BestVar=VarName[PLSvar]
    print(PLSvar)
    X<-DATA[,PLSvar]
    colnames(X)<-VarName[PLSvar]
    YX <- data.frame(Y,X)
    StoreLDA<-lda(X,Y,tol = 1.0e-4,method="moment") 
    if(Mes==1){
      Accuracy.rate.VarSelct <- function(formula, data, indices){
        Modlda <- lda(formula, data, CV = TRUE)
        tab <- table(Modlda$class, data$Y)
        num <- sum(diag(tab))
        denom <- sum(tab)
        error <- signif(num/denom, 4)
        return(error)
      }
    }else{
      Accuracy.rate.VarSelct <- function(formula, data, indices){
        Modlda<-lda(formula, data, CV = TRUE)
        Class<-table(data$Y)
        ClassTest<-table(Modlda$class,data$Y)
        error<-signif(mean(diag(ClassTest)/Class), 4)  
        return(error)
      } # MisClass. 
    }
   boot.out <- boot(data = YX, statistic =  Accuracy.rate.VarSelct  , R = 200, formula = Y~X)
   Error <- mean(boot.out$t)
 }
  if (Method==2){
    Aux <-pcaSVD(DATA,prop) 
    nc<-length(DATA[1,])
    pls.out <- pls.regression(Xtrain = DATA, Ytrain = Y, ncomp = Aux$n)
    T <- pls.out$T
    pls.all <- pls.regression(Xtrain = DATA, Ytrain = Y, ncomp = nc)
    StoreLDA<-lda(T,Y,tol = 1.0e-4,method="moment") 
    YX <- data.frame(Y, T)
    X<-pls.all$T
    W<-pls.all$R
    RdYTm<-rep(0,nc)
    RdYTh<-rep(0,nc)
    VIP<-rep(0,nc)
    for (i in 1:nc){
      for(j in 1:nc){
        RdYTm[j]<-(cor(Y,X[,j]))^2
        RdYTh[j]<-RdYTm[j]*(W[i,j])^2
      }
      VIP[i]<-sqrt((nc/sum(RdYTm))*(sum(RdYTh)))
    }
    #VIP<-VIP/sqrt(sum(VIP^2)) 
   
    if(Mes==1){
      Accuracy.rate <- function(formula, data, indices){
        Modlda <- lda(formula, data, CV = TRUE)
        tab <- table(Modlda$class, data$Y)
        num <- sum(diag(tab))
        denom <- sum(tab)
        error <- signif(num/denom, 4)
        return(error)
      }
    }else{
      Accuracy.rate <- function(formula, data, indices){
        Modlda<-lda(formula, data, CV = TRUE)
        Class<-table(data$Y)
        ClassTest<-table(Modlda$class,data$Y)
        error<-signif(mean(diag(ClassTest)/Class), 4)  
        return(error)
      } # MisClass. 
    }
    boot.out <- boot(data = YX, statistic = Accuracy.rate, R = 200, formula = Y~T)
    Error <- mean(boot.out$t)
    BestVar=NULL
  }
 
  if (Method==3){
    Aux <-pcaSVD(DATA,prop) 
    X<-Aux$u[,1:Aux$n]
    StoreLDA<-lda(X,Y,tol = 1.0e-4,method="moment") 
    YX <- data.frame(Y, X)
    if(Mes==1){
      Accuracy.rate.PCA<- function(formula, data, indices){
        Modlda <- lda(formula, data, CV = TRUE)
        tab <- table(Modlda$class, data$Y)
        num <- sum(diag(tab))
        denom <- sum(tab)
        error <- signif(num/denom, 4)
        return(error)
      }
    }else{
      Accuracy.rate.PCA <- function(formula, data, indices){
        Modlda<-lda(formula, data, CV = TRUE)
        Class<-table(data$Y)
        ClassTest<-table(Modlda$class,data$Y)
        error<-signif(mean(diag(ClassTest)/Class), 4)  
        return(error)
      } # MisClass. 
    }
    boot.out <- boot(data = YX, statistic = Accuracy.rate.PCA, R = 200, formula = Y~X)
    Error <- mean(boot.out$t)
    BestVar=NULL
  }
  if (Method!=2){
    PLS.VIP<-NULL
  }else{
    VarRank<-sort.int(VIP,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix
    PLS.VIP<-data.frame(VAR=VarName[VarRank],VIP=VIP[VarRank])
  }
  return(list(LDA=StoreLDA,BestVar=BestVar,Accuracy=Error,VIP=PLS.VIP))
}
pcaSVD <-function(DATA,prop) 
{
  # Entradas
  # DATA    : matriz de datos
  # Cmin, Cmax
  # prop : proporción de varianza acumulda
  X<-as.matrix(DATA)                                   # Transformar la variable DATA  a Tipo de dato matriz
  XScale <- scale(X, center = TRUE, scale = TRUE)      # Escalar y Centrar la  matriz de datos
  xSvd <- svd(XScale)                                  # Descomposición del valor singular de la matriz de datos XScale       
  D <- xSvd$d                                          # Matriz de cargas  de las componentes principales
  U <- xSvd$u                                          # Matriz de scores  de las componentes principales
  V <- xSvd$v                                          # Matriz diagonal que contiene las desvia
  propVar <- D^2/sum(D^2)                              # Método empírico para detectar sepeccionar nímero de componentes  
  propAcum <- cumsum(propVar)
  n <- length(propAcum[propAcum <= prop])              # cantidad de componentes a seleccionar   
  return(list(xSvd, d = D, u = U, v = V, n = n,XScale=XScale))
  # Salidas
  # d : desviaciones de las componentes principales (sqrt(varianza))
  # u : matriz de vectores singulares izquierdos
  # v : matriz de vectores singulares
  # n :numero de componentes seleccionadas que resumen la proporcion de varianza deseada
}   