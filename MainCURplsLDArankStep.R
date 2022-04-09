Main.CUR.PLS.LDA.rank<-function(DATA,Y,VarName,k,w,IterMax,met.cur,met.class,mes.class){
  # HASVC : algoritmo para la selección de variables en la clasificación
  # 
  # ENTRADAS
  # DATA   : matriz de atributos
  # Y       : Vector de clases
  # VarName : Nombre de la variables
  # prop    : proporción de varianza acumulada por las componentes
  # w1,w2   : pesos para ponderar la frecuencia y el error al factor importancia 
  # IterMax : Número máximo de iteraciones
  # met.cur : método usado para usadas para seleccionar columnas en CUR: met.cur<-1 para met.pls=2
  # met.pls : si 1 aplicar PLS, si 2 aplicar var.selection (prop=[val1,val2] donde val=0.95, val2=[0.6:0.8]) y met.cur=1, si 3 no plica pls 
  # met.rank: metodo para la selección de variables mediante el ranking (1 umbral, 2 número de variables)
  library("rCUR")
  library("boot")
  library("plsgenomics")
  library("klaR")
  library("e1071")
  #Inicialización
  met.rank<-2
  m<-dim(DATA)[1]
  n<-dim(DATA)[2]
  # VarName<-colnames(DATA[,n-2])
  Hfvar<-rep(0,n)                                        # almacena la frecuencia  histórica de selección  de cada variable (Ratio) 
  HEbcv<-rep(0,n)                                        # almacena el errorBCV histórico de los modelos en a ha estado presente una variable 
  HIvar<-rep(0,n)                                        # almacena la importancia de la variable dentro de los modelos en los cuales ha estado presente
  Hfimp<-rep(0,n)                                        # almacena el factor de permanencia para el filtrado
  colnames(DATA)<-VarName
  Iter<-1; Stop0<-0; Stop1<-0;Stop2<-0                   # condiciones de parada (inicilización)
  pca<-pcaSVD(DATA,n-2)                                  # Descomposición del valor singular
  VarSelectAnt<-c(0)#rep(0,k)
  DATA<-pca$XScale
  A<-BuildModelCUR(DATA,Y,pca$xSvd,k,met.cur)
  VarSelect<-A$SVar                             # seleción inicial de variables                                      # ratios historicos
  for(i in 1:length(VarSelect)){
    j<-VarSelect[i]
    Hfvar[j]<-Hfvar[j]+1                                 # cáculo de la frecuencia de aparición del ratio 
    HEbcv[j]<-HEbcv[j]+A$FitnessBCV                      # cálculo del factor de importancia del ratio de acuerdo al modelo                                
    HIvar[j]<-(HIvar[j]+A$VIP[i])
  }
  while(Stop0<1){
    Iter<-Iter+1
    print(Iter)
    A<-BuildModelCUR(DATA,Y,pca$xSvd,k,met.cur)
    VarSelect<-A$SVar                                                          # seleción inicial de variables                                      
    for(i in 1:length(VarSelect)){
      j<-VarSelect[i]
      Hfvar[j]<-Hfvar[j]+1                                                      # cáculo de la frecuencia de aparición del ratio 
      HEbcv[j]<-HEbcv[j]+A$FitnessBCV                                           # cálculo del factor de importancia del ratio de acuerdo al modelo                                
      HIvar[j]<-(HIvar[j]+A$VIP[i])
    }
    # condición de parada
    if(Iter==IterMax){Stop0<-1}# Stop0
  } # while
  HH<-matrix(0,nrow=n,ncol=4)
  HH[,1]<-Hfvar
  for (y in 1:length(HH[,1])){
    if(Hfvar[y]!=0){
      HH[y,2]<-HEbcv[y]/Hfvar[y]
      HH[y,3]<-HIvar[y]/Hfvar[y]
      Hfimp[y]<-w*((Hfvar[y]/IterMax)+HEbcv[y]/Hfvar[y])+(1-w)*((HIvar[y]/Hfvar[y]))
      HH[y,4]<-Hfimp[y]
    }
  }
  VarRank<-sort.int(Hfimp,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables según su balor de factor de importancia
  rownames(HH)<-VarName[VarRank];HH[,1]<-HH[VarRank,1];HH[,2]<-HH[VarRank,2];HH[,3]<-HH[VarRank,3];HH[,4]<-HH[VarRank,4]
  colnames(HH)<-c("Frec","M.EBCV","VIP","Rank")
  # Selección de variables mediante ranking
  met.rank<-1
  if(met.rank==1){
    # metodo de umbral
    P<-Hfimp
    VarRank<-(1:length(P)) 
    a<- quantile(P,names=FALSE,na.rm=TRUE)
    U<-median(P)
    VarSelect<-VarRank[P>U]                       
  }
  if(met.rank==2){
    # metodo de número de variables fijo
    P<- HH[y,4]
    VarRank<-sort.int(P,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables según su balor de factor de importancia

    VarSelect<-VarRank[1:k] 
  } # if else (met.rank)
  
  # Métodos de clasificación
  if(met.class==1){method="lda"};if(met.class==2){method="qda"};if(met.class==3){method="rda"}
  if(met.class==4){method="sknn"};if(met.class==5){method="svm"};if(met.class==6){method="NaiveBayes"}
  
  # Criterio de evaluación
  if(mes.class==1){criterion="CR"};if(mes.class==2){criterion="AC"}
  if(mes.class==3){criterion="AS"};if(mes.class==4){criterion="CF"}
  Var.Selection<-stepclass(x=DATA[,VarSelect],grouping=Y,method=method,improvement = 0.0005,maxvar = k, criterion = criterion)
  VarSelect<- Var.Selection$model$name
  k<-length(VarSelect)
  SVarCUR<-rep(0,k)
  for(i in 1:k){SVarCUR[i]<-which(VarName==VarSelect[i])}

  Aux<-PLS.class.error<-PLS.LDA.errorBCV(DATA,Y,SVarCUR,k)
  
  return(list(RANKING=HH,Umbral=U,Summary=Var.Selection,PLS.class.error.boot=Aux$Fitness))
  # Salidas
  # BestVar: variables del mejor modelo Aparición histórica de los ratios en los modelos generados
  # ZModel : Modelo descriminante
  # EBCV   : Presicion  mejor modelo  
} # function
#----------------------------------------------------------------------------------------------------------------
# Descompocisión del Valor Singular
pcaSVD <-function(DATA,NVar) 
{
  # Entradas
  # DATA    : matriz de datos
  # Cmin, Cmax
  # prop : proporción de varianza acumulda
  X<-as.matrix(DATA[,1:NVar])                       # Transformar la variable DATA  a Tipo de dato matriz
  XScale <- scale(X, center = TRUE, scale = TRUE)      # Escalar y Centrar la  matriz de datos
  xSvd <- svd(XScale)                                  # Descomposición del valor singular de la matriz de datos XScale       
  D <- xSvd$d                                          # Matriz de cargas  de las componentes principales
  U <- xSvd$u                                          # Matriz de scores  de las componentes principales
  V <- xSvd$v                                          # Matriz diagonal que contiene las desvia
  return(list(xSvd = xSvd, d = D, u = U, v = V,XScale = XScale))
  # Salidas
  # d : desviaciones de las componentes principales (sqrt(varianza))
  # u : matriz de vectores singulares izquierdos
  # v : matriz de vectores singulares
  # n :numero de componentes seleccionadas que resumen la proporcion de varianza deseada
}
#----------------------------------------------------------------------------------------------------------------
BuildModelCUR<-function(DATA,Y,xSvd,k,met.cur)
{
  # DATA    : Matriz de datos (Emp X Ratios)
  # Nvar    : Número de variables a tener en cuenta
  # Y       : Vector de clases
  # met.pls : método de seleccion de columnas usado en algoritmo ColumnSelect (CUR)
  # met.pls : Si se aplica PLS par construción del modelo
  # xSvd    : descompocision del valor singular
  # NC      : numero de columnas a Selecionar
  Stop<-0
  Dmatriz <- as.matrix(DATA)    # Transformar la variable DATA  a Tipo de dato matriz                                         
  # método de selección
  if(met.cur==1){ 
    method<-"random"
    alpha = 1
  } # Met 1
  if(met.cur==2){
    method<-"exact.num.random" 
    alpha = 1
  } 
  if(met.cur==3){
    method<-"top.scores" 
    alpha = 1
  } 
  #descomposición matricial CUR 
  while(Stop==0){
    CURResults <- CUR(Dmatriz, c= k, r= dim(Dmatriz)[1], k =k, sv = xSvd ,method,alpha, weighted = FALSE, beta = 4,matrix.return = TRUE) 
    SVarCUR <- CURResults@C.index              # variables Seleccionados                  
    if(length(SVarCUR)>1){
      Stop<-1
    }
  }
  Aux<-PLS.LDA.errorBCV(Dmatriz,Y,SVarCUR,k)
  return(list(SVar=SVarCUR,Dmatriz=Dmatriz,FitnessBCV=Aux$Fitness,VIP=Aux$VIP))                  
  # Salidas
  # SRatios: Ratios seleccionados para el modelo
  # CoefModel:coeficientes del modelo descriminante
  # FitnessBCV: error de boostrap clasificación del modelo
}
#----------------------------------------------------------------------------------------------------------------

PLS.LDA.errorBCV<-function(Dmatriz,Y,SVarCUR,k){
  # Entrada
  # X: Subconjunto de la matriz  de datos seleccionado mediante la descomposición CUR
  # Y: variable dependiente (contiene las clases)
  
  # met.pls=1 Calcula las componetes PLS de SVarCUR para calcular el modelo descriminante 
    XX<-Dmatriz[,SVarCUR]
    pls.out <- pls.regression(Xtrain =XX, Ytrain = Y, ncomp = k)
    X <- pls.out$T
    W <- pls.out$R
    colnames(X)<-colnames(Dmatriz[,SVarCUR])
    YX <- data.frame(Y, X)
    RdYTm<-rep(0,k)
    RdYTh<-rep(0,k)
    VIP<-rep(0,k)
    for (i in 1:k){
      for(j in 1:k){
        RdYTm[j]<-(cor(Y,X[,j]))^2
        RdYTh[j]<-RdYTm[j]*(W[i,j])^2
      }
      VIP[i]<-sqrt((k/sum(RdYTm))*(sum(RdYTh)))
    }
    VIP<-VIP/sqrt(sum(VIP^2))

  # met.pls=2 no  Calcula las componentes PLS de SVarCUR se calcula el modelo descriminante directamente SVarCUR
 
  Mes<-1
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
      aux1<-diag(ClassTest)/Class
      aux2<-mean(aux1)
      error<-signif(aux2, 4)  
      return(error)
    } # MisClass.
  }
  boot.out <- boot(data = YX, statistic = Accuracy.rate , R = 200, formula = Y~X)
  Error <- mean(boot.out$t)
  return(list(Fitness=Error,VIP=VIP))
  # Salidas
  # Fitness: precisión del modelo descriminate
  # Coef   : modelo descriminate 
}
