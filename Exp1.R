Exp1<-function(DataExp){
  #-------------------------------------#
  # Paquetes Requeridos
  #-------------------------------------#
  library(boot)
  library(lattice)
  library(MASS)
  library(Matrix)
  library(rCUR)
  #-------------------------------------#
  # Par?metros del Experimento
  #-------------------------------------#
  k <-c(4,6,8)                          # cantidad de variables en la solucion k-columnas
  N=25                                  # Canttidad de corridas para una configuraci?n
  w<-0.5                                # n?mero de modelo a generar
  #-------------------------------------#
  for (DS in 1:3){                      # Para cada uno de los conjuntos de datos        
  
    if(DS==1){ #DataSet 1
       #-------------------------------------#
       print("Pietruszkiewicz Dataset")
       #-------------------------------------#
       VAR<-DataExp$VAR1[1:22]
       Col<-1:length(VAR)
       DATA<-DataExp$DataSet1[,1:22]
       Y<-DataExp$DataSet1[,23]
       NV<-length(VAR)                       # cantidad de variables           
       Nmet  = 2                             # metodo de para ColumnSelect (ver funcion column select)                                    
       IterMax=200
       #---------------------------------------------------------------------------#
       #Variables auxiliares
       HH1<-matrix(0,nrow=NV,ncol=3); HH2<-matrix(0,nrow=NV,ncol=3); HH3<-matrix(0,nrow=NV,ncol=3)          #almacena el ranking hist?tico de las variables de una corrida
       # Almacena la frecuencia de las variables en la soluci?n final
       StoreEBCV<-matrix(0,nrow=N+2,ncol=length(k));    
       colnames(StoreEBCV)<-c("k4","k6","k8")   
       # Almacena la frecuencia de las variables en la soluci?n final
       BestVARfrec<-matrix(0,nrow=NV,ncol=length(k));  
       rownames(BestVARfrec)<-VAR;colnames(BestVARfrec)<-c("BVF1","BVF2","BVF3")
       HBestVAR<-data.frame(cbind(k4=0,k6=0,k8=rep(0,max(k+1))))
       RANK1<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       RANK2<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       RANK3<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       #EXPERIMENTACION#
       #----------------------------------------------------------------------------#
       for(i in 1:length(k)){
         BestEBCV<-0
         for(l in 1:N){
           print("Ejecuci?n")
           print(c(DS,i,l))
           FunctionData<-Main.CUR.LDA(DATA,Y,VAR,k[i],w,IterMax=IterMax,2)
           StoreEBCV[l,i]<-FunctionData$EBCV                                  # mejor precisi?n del mejor modelo de la ejecuci?n l para la cantidad de variables k[i]
           AuxVar<-rownames(FunctionData$ZModel)                             # nombre de variables del modelo enla corrida l
           if(FunctionData$EBCV<BestEBCV){                                   # actualizar el  mejor modelo  hasta el momento                                                                                                                         
             HBestVAR[2:length(AuxVar),i]<-rownames(FunctionData$ZModel)
             HBestVAR[1,i]<-FunctionData$EBCV
             BestEBCV<-FunctionData$EBCV
           }
          
           for(q in 1:length(AuxVar)){
             j<-which(VAR==AuxVar[q])
             BestVARfrec[j,i]<-BestVARfrec[j,i]+1                              # c?culo de la frecuencia de aparici?n de la variable en el mejor modelo de una ejecuci?n                                                                     # c?lculo del factor de permanencia (FALTA ->  w3*HIvar[j])
           }
           if(i==1){
             HH1<-HH1+FunctionData$HH
             if(l==N){
               Frec<-HH1[,1]
               M.Error<-HH1[,2]
               Rank<-HH1[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK1[,1]<-VAR[VarRank]
               RANK1[,2]<-Frec[VarRank]
               RANK1[,3]<-(M.Error[VarRank])/N
               RANK1[,4]<-(Rank[VarRank])/N
               write.matrix(RANK1, file="RANK11.csv",sep=";")
             }
           }
           if(i==2){
             HH2<-HH2+FunctionData$HH 
             if(l==N){
               Frec<-HH2[,1]
               M.Error<-HH2[,2]
               Rank<-HH2[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK2[,1]<-VAR[VarRank]
               RANK2[,2]<-Frec[VarRank]
               RANK2[,3]<-(M.Error[VarRank])/N
               RANK2[,4]<-(Rank[VarRank])/N
               write.matrix(RANK2, file="RANK12.csv",sep=";")
             }
           }
           if(i==3){
             HH3<-HH3+FunctionData$HH
             if(l==N){
               Frec<-HH3[,1]
               M.Error<-HH3[,2]
               
               Rank<-HH3[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK3[,1]<-VAR[VarRank]
               RANK3[,2]<-Frec[VarRank]
               RANK3[,3]<-(M.Error[VarRank])/N
               RANK3[,4]<-(Rank[VarRank])/N
               write.matrix(RANK3, file="RANK13.csv",sep=";")
             }
           }
           write.matrix(StoreEBCV, file="StoreEBCV1.csv",sep=";")
           write.matrix(BestVARfrec, file="BestVARfrec1.csv",sep=";")
          
         } # for cantidad de corridas N
         StoreEBCV[N+1,i]<-mean(StoreEBCV[1:N,i]) 
         StoreEBCV[N+2,i]<-sd(StoreEBCV[1:N,i]) 
         Store<-list(RANK1=RANK1,RANK2=RANK2,RANK3=RANK3,StoreEBCV=StoreEBCV,BestVARfrec=BestVARfrec)
         save(Store,file="Store1.RData") # Salva resultados para dataset1
         } # for (k[i] cantidad de variables)
     } # if DataSet1
  ############################################################################################################################################
     if(DS==2){
       #-------------------------------------#
       print("Du Jardin Dataset")
       #-------------------------------------#
       VAR<-DataExp$VAR2[1:39]
       Col<-1:length(VAR)
       DATA<-DataExp$DataSet2[,1:39]
       Y<-DataExp$DataSet2[,40]
       NV<-length(VAR)                       # cantidad de variables           
       Nmet  = 2                             # m?todo de para ColumnSelect                                       
       IterMax=300
       #-------------------------------------#
       HH1<-matrix(0,nrow=NV,ncol=3); HH2<-matrix(0,nrow=NV,ncol=3); HH3<-matrix(0,nrow=NV,ncol=3)          #almacena el ranking hist?tico de las variables de una corrida
       # Almacena la frecuencia de las variables en la soluci?n final
       StoreEBCV<-matrix(0,nrow=N+2,ncol=length(k));    
       colnames(StoreEBCV)<-c("k5","k7","k9")   
       # Almacena la frecuencia de las variables en la soluci?n final
       BestVARfrec<-matrix(0,nrow=NV,ncol=length(k));  
       rownames(BestVARfrec)<-VAR;colnames(BestVARfrec)<-c("BVF1","BVF2","BVF3")
       HBestVAR<-data.frame(cbind(k4=0,k6=0,k8=rep(0,max(k+1))))
       RANK1<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       RANK2<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       RANK3<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       #EXPERIMENTACION#
       #----------------------------------------------------------------------------
       for(i in 1:length(k)){
         BestEBCV<-0
         for(l in 1:N){
           print("Ejecuci?n")
           print(c(DS,i,l))
           FunctionData<-Main.CUR.LDA(DATA,Y,VAR,k[i],w,IterMax=IterMax,2)
           StoreEBCV[l,i]<-FunctionData$EBCV                                  # mejor precisi?n del mejor modelo de la ejecuci?n l para la cantidad de variables k[i]
           AuxVar<-FunctionData$BestVar                                       # nombre de variables del modelo enla corrida l
           if(FunctionData$EBCV<BestEBCV){                                   # actualizar el  mejor modelo  hasta el momento                                                                                                                         
             HBestVAR[2:length(AuxVar),i]<-rownames(FunctionData$ZModel)
             HBestVAR[1,i]<-FunctionData$EBCV
             BestEBCV<-FunctionData$EBCV
           }
           for(q in 1:length(AuxVar)){
             j<-which(VAR==AuxVar[q])
             BestVARfrec[j,i]<-BestVARfrec[j,i]+1                              # c?culo de la frecuencia de aparici?n de la variable en el mejor modelo de una ejecuci?n                                                                     # c?lculo del factor de permanencia (FALTA ->  w3*HIvar[j])
           }
           if(i==1){
             HH1<-HH1+FunctionData$HH
             if(l==N){
               Frec<-HH1[,1]
               M.Error<-HH1[,2]
               Rank<-HH1[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK1[,1]<-VAR[VarRank]
               RANK1[,2]<-Frec[VarRank]
               RANK1[,3]<-(M.Error[VarRank])/N
               RANK1[,4]<-(Rank[VarRank])/N
               write.matrix(RANK1, file="RANK21.csv",sep=";")
             }
           }
           if(i==2){
             HH2<-HH2+FunctionData$HH 
             if(l==N){
               Frec<-HH2[,1]
               M.Error<-HH2[,2]
               
               Rank<-HH2[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK2[,1]<-VAR[VarRank]
               RANK2[,2]<-Frec[VarRank]
               RANK2[,3]<-(M.Error[VarRank])/N
               RANK2[,4]<-(Rank[VarRank])/N
               write.matrix(RANK2, file="RANK22.csv",sep=";")
             }
           }
           if(i==3){
             HH3<-HH3+FunctionData$HH
             if(l==N){
               Frec<-HH3[,1]
               M.Error<-HH3[,2]
               Rank<-HH3[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK3[,1]<-VAR[VarRank]
               RANK3[,2]<-Frec[VarRank]
               RANK3[,3]<-(M.Error[VarRank])/N
               RANK3[,4]<-(Rank[VarRank])/N
               write.matrix(RANK3, file="RANK23.csv",sep=";")
             }
           }                                                                               
           write.matrix(StoreEBCV, file="StoreEBCV2.csv",sep=";")
           write.matrix(BestVARfrec, file="BestVARfrec2.csv",sep=";")
           Store<-list(RANK1=RANK1,RANK2=RANK2,RANK3=RANK3,StoreEBCV=StoreEBCV,BestVARfrec=BestVARfrec)
           save(Store,file="Store2.RData")
         } # for cantidad de comidas.
         StoreEBCV[N+1,i]<-mean(StoreEBCV[1:N,i]) 
         StoreEBCV[N+2,i]<-sd(StoreEBCV[1:N,i])           
         } # for (k[i] cantidad de variables)
       } #  DataSet 2
     ####################################################################################################
     if(DS==3){
       #-------------------------------------#
       print(" Atilya DataSet")
       #-------------------------------------#
       VAR<-DataExp$VAR3[1:60]   
       DATA<-DataExp$DataSet3[,1:60]
       Y<-DataExp$DataSet3[,61]              # cantidad de variables           
       NV=length(VAR)
       Nmet  = 2                             # m?todo de para ColumnSelect                                     
       IterMax=500
       #-------------------------------------#
       HH1<-matrix(0,nrow=NV,ncol=3); HH2<-matrix(0,nrow=NV,ncol=3); HH3<-matrix(0,nrow=NV,ncol=3)          #almacena el ranking hist?tico de las variables de una corrida
       # Almacena la frecuencia de las variables en la soluci?n final
       StoreEBCV<-matrix(0,nrow=N+2,ncol=length(k));    
       colnames(StoreEBCV)<-c("k6","k9","k15")   
       # Almacena la frecuencia de las variables en la soluci?n final
       BestVARfrec<-matrix(0,nrow=NV,ncol=length(k));  
       rownames(BestVARfrec)<-VAR;colnames(BestVARfrec)<-c("BVF1","BVF2","BVF3")
       HBestVAR<-data.frame(cbind(k4=0,k6=0,k8=rep(0,max(k+1))))
       RANK<-data.frame(cbind(VarName=0,Frec=0,ME=0,Rank=rep(0,NV)))
       #EXPERIMENTACION#
       #----------------------------------------------------------------------------
       for(i in 1:length(k)){
         for(l in 1:N){
           print("Ejecuci?n")
           print(c(DS,i,l))
           FunctionData<-Main.CUR.LDA(DATA,Y,VAR,k[i],w,IterMax=IterMax,2)
           StoreEBCV[l,i]<-FunctionData$EBCV                                  # mejor precisi?n del mejor modelo de la ejecuci?n l para la cantidad de variables k[i]
           AuxVar<-FunctionData$BestVar                                       # nombre de variables del modelo enla corrida l
         
           if(FunctionData$EBCV<HBestEBCV){                                   # actualizar el  mejor modelo  hasta el momento                                                                                                                         
             HBestVAR[2:length(AuxVar),i]<-rownames(FunctionData$ZModel)
             HBestVAR[1,i]<-FunctionData$EBCV
           }
           for(q in 1:length(AuxVar)){
             j<-which(VAR==AuxVar[q])
             BestVARfrec[j,i]<-BestVARfrec[j,i]+1                              # c?culo de la frecuencia de aparici?n de la variable en el mejor modelo de una ejecuci?n                                                                     # c?lculo del factor de permanencia (FALTA ->  w3*HIvar[j])
           }
           if(i==1){
             HH1<-HH1+FunctionData$HIST
             if(l==N){
               Frec<-HH1[,1]
               M.Error<-HH1[,2]
               Rank<-HH1[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK1[,1]<-VAR[VarRank]
               RANK1[,2]<-Frec[VarRank]
               RANK1[,3]<-(M.Error[VarRank])/N
               RANK1[,4]<-(Rank[VarRank])/N
               write.matrix(RANK1, file="RANK31.csv",sep=";")
             }
           }
           if(i==2){
             HH2<-HH2+FunctionData$HIST 
             if(l==N){
               Frec<-HH2[,1]
               M.Error<-HH2[,2]
               Rank<-HH2[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK2[,1]<-VAR[VarRank]
               RANK2[,2]<-Frec[VarRank]
               RANK2[,3]<(-M.Error[VarRank])/N
               RANK2[,4]<-(Rank[VarRank])/N
               write.matrix(RANK2, file="RANK32.csv",sep=";")
             }
           }
           if(i==3){
             HH3<-HH3+FunctionData$HIST
             if(l==N){
               Frec<-HH3[,1]
               M.Error<-HH3[,2]
               Rank<-HH3[,3]
               VarRank<-sort.int(Rank,partial = NULL, na.last = NA, decreasing = TRUE, method = c("shell", "quick"), index.return = TRUE)$ix # ordenar variables seg?n su balor de factor de importancia
               RANK3[,1]<-VAR[VarRank]
               RANK3[,2]<-Frec[VarRank]
               RANK3[,3]<-(M.Error[VarRank])/N
               RANK3[,4]<-(Rank[VarRank])/N
               write.matrix(RANK3, file="RANK33.csv",sep=";")
             }
           }                                                                        
           write.matrix(StoreEBCV, file="StoreEBCV3.csv",sep=";")
           write.matrix(BestVARfrec, file="BestVARfrec3.csv",sep=";")
        }# for  la cantidad de corridas
         StoreEBCV[N+1,i]<-mean(StoreEBCV[1:N,i]) 
         StoreEBCV[N+2,i]<-sd(StoreEBCV[1:N,i]) 
         Store<-list(RANK1=RANK1,RANK2=RANK2,RANK3=RANK3,StoreEBCV=StoreEBCV,BestVARfrec=BestVARfrec)
         save(Store,file="Store3.RData")
      } # for (k[i] cantidad de variables)
    } # if  DataSet 3 
  } # Para cada uno de los conjuntos de datos
}# Function Exp1
###############################################################################################
Main.CUR.LDA<-function(DATA,Y,VarName,k,w,IterMax,met.cur){
  # HASVC-CUR.LDA : algoritmo para la selecci?n de k variables en la clasificaci?n
  # 
  # ENTRADAS
  # DATA    : matriz de atributos
  # Y       : Vector de clases
  # VarName : Nombre de la variables
  # k       : cantidad de variables
  # w       : peso para ponderar la frecuencia(w) y el error (1-w) en el factor importancia para crear ranking
  # IterMax : N?mero m?ximo de iteraciones
  # met.cur : m?todo usado para usadas para seleccionar columnas en CUR: met.cur<-1 :Random,  met.cu<-2 
  #library("rCUR");library("boot");library("rCUR");library("plsgenomics")
  #########################################################################################################################
  #Inicializaci?n
  m<-dim(DATA)[1]
  n<-dim(DATA)[2]
  Hfvar<-rep(0,n)                                        # almacena la frecuencia  hist?rica de selecci?n  de cada variable (Ratio) 
  HEbcv<-rep(0,n)                                        # almacena el errorBCV hist?rico de los modelos en a ha estado presente una variable 
  HIvar<-rep(0,n)                                        # almacena la importancia de la variable dentro de los modelos en los cuales ha estado presente
  Hfimp<-rep(0,n)                                        # almacena el factor de permanencia para el filtrado
  colnames(DATA)<-VarName
  Iter<-1; Stop0<-0; Stop1<-0;Stop2<-0                    # condiciones de parada (inicilizaci?n)
  pca<-pcaSVD(DATA,n)                                     # Descomposici?n del valor singular
  DATA<-pca$XScale
  A<-BuildModelCUR(DATA,Y,pca$xSvd,k,met.cur)                       # seleci?n inicial de variables                                      # ratios historicos
  VarSelect<-rownames(A$CoefModel)                                            # seleci?n inicial de variables    
  for(i in 1:length(VarSelect)){
    j<-which(VarName==VarSelect[i])
    Hfvar[j]<-Hfvar[j]+1                                                      # c?culo de la frecuencia de aparici?n del ratio 
    HEbcv[j]<-HEbcv[j]+A$FitnessBCV                                           # c?lculo del factor de importancia del ratio de acuerdo al modelo                                
    Hfimp[j]<-Hfimp[j]+(w*((Hfvar[j])/IterMax)+(1-w)*((HEbcv[j])/(Hfvar[j]))) # c?lculo del factor de permanencia (FALTA ->  w3*HIvar[j])
  }
  BestVar <-rownames(A$CoefModel)                                         # mejores ratios hasta el momento
  BestCoef<-A$CoefModel                                                     # mejores Coeficientes
  BestFitnessBCV<-A$FitnessBCV                                              # Mejor Error
  #############
  while(Stop0<1){
    Iter<-Iter+1 
    A<-BuildModelCUR(DATA,Y,pca$xSvd,k,met.cur)
    if(A$FitnessBCV>BestFitnessBCV){                                            # actualizar el  mejor modelo  hasta el momento                                                                                                                         
      BestCoef<-A$CoefModel
      BestVar<-rownames(BestCoef)
      BestFitnessBCV<-A$FitnessBCV
    }                                                                           # ratios historicos
    VarSelect<-rownames(A$CoefModel)                                            # seleci?n inicial de variables    
    for(i in 1:length(VarSelect)){
      j<-which(VarName==VarSelect[i])
      Hfvar[j]<-Hfvar[j]+1                                                      # c?culo de la frecuencia de aparici?n del ratio 
      HEbcv[j]<-HEbcv[j]+A$FitnessBCV                                           # c?lculo del factor de importancia del ratio de acuerdo al modelo                                
      Hfimp[j]<-Hfimp[j]+(w*((Hfvar[j])/IterMax)+(1-w)*((HEbcv[j])/(Hfvar[j]))) # c?lculo del factor de permanencia (FALTA ->  w3*HIvar[j])
    }
    # condici?n de parada
    if(Iter==IterMax){
      Stop1<-1
    } # Stop1
    if(BestFitnessBCV>0.99){
      Stop2<-1
    } # Stop2
    Stop0<-Stop1+Stop2
  } # while
  HH<-matrix(rep(0,3*length(VarName)),ncol=3)
  HH[,1]<-Hfvar
  for (y in 1:length(HH[,1])){
    if(Hfvar[y]!=0){
      HH[y,2]<-HEbcv[y]/Hfvar[y]
      HH[y,3]<-Hfimp[y]/Hfvar[y]
    }
  }
  rownames(HH)<-VarName
  colnames(HH)<-c("Frec","M.EBCV","Rank")
  return(list(BestVar=BestVar,ZModel=BestCoef,EBCV=BestFitnessBCV,HH=HH))
  # Salidas
  # BestVar: variables del mejor modelo Aparici?n hist?rica de los ratios en los modelos generados
  # ZModel : Modelo descriminante
  # EBCV   : Presicion  mejor modelo  
} # function
#----------------------------------------------------------------------------------------------------------------
# Descompocisi?n del Valor Singular
pcaSVD <-function(DATA,NVar) 
{
  # Entradas
  # DATA    : matriz de datos
  # Cmin, Cmax
  # prop : proporci?n de varianza acumulda
  X<-as.matrix(DATA[,1:NVar])                       # Transformar la variable DATA  a Tipo de dato matriz
  XScale <- scale(X, center = TRUE, scale = TRUE)      # Escalar y Centrar la  matriz de datos
  xSvd <- svd(XScale)                                  # Descomposici?n del valor singular de la matriz de datos XScale       
  D <- xSvd$d                                          # Matriz de cargas  de las componentes principales
  U <- xSvd$u                                          # Matriz de scores  de las componentes principales
  V <- xSvd$v                                          # Matriz diagonal que contiene las desvia
  # cantidad de componentes a seleccionar   
  return(list(xSvd = xSvd, d = D, u = U, v = V,XScale = XScale))
  # Salidas
  # d : desviaciones de las componentes principales (sqrt(varianza))
  # u : matriz de vectores singulares izquierdos
  # v : matriz de vectores singulares
}
#----------------------------------------------------------------------------------------------------------------
BuildModelCUR<-function(DATA,Y,xSvd,k,met.cur)
{
  # DATA    : Matriz de datos (Emp X Ratios)
  # Nvar    : N?mero de variables a tener en cuenta
  # Y       : Vector de clases
  # met.pls : m?todo de seleccion de columnas usado en algoritmo ColumnSelect (CUR)
  # met.pls : Si se aplica PLS par construci?n del modelo
  # xSvd    : descompocision del valor singular
  # NC      : numero de columnas a Selecionar
  Stop<-0
  Dmatriz <- as.matrix(DATA)    # Transformar la variable DATA  a Tipo de dato matriz                                         
  # m?todo de selecci?n
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
  #descomposici?n matricial CUR 
  while(Stop==0){
    CURResults <- CUR(Dmatriz, c= k, r= dim(Dmatriz)[1], k = k, sv = xSvd ,method,alpha, weighted = FALSE, beta = 4,matrix.return = TRUE) 
    SVarCUR <- CURResults@C.index              # variables Seleccionados                  
    if(length(SVarCUR)>1){
      Stop<-1
    }
  }
  Aux<-PLS.LDA.errorBCV(Dmatriz,Y,SVarCUR,k)
  return(list(SVar=SVarCUR,Dmatriz=Dmatriz,CoefModel=Aux$Coef,FitnessBCV=Aux$Fitness,VIP=Aux$VIP))                  
  # Salidas
  # SRatios: Ratios seleccionados para el modelo
  # CoefModel:coeficientes del modelo descriminante
  # FitnessBCV: error de boostrap clasificaci?n del modelo
}
#----------------------------------------------------------------------------------------------------------------
PLS.LDA.errorBCV<-function(Dmatriz,Y,SVarCUR,k){
  # Entrada
  # X: Subconjunto de la matriz  de datos seleccionado mediante la descomposici?n CUR
  # Y: variable dependiente (contiene las clases)
  X<-Dmatriz[,SVarCUR]
  nc<-length(SVarCUR)
  StoreLDA<-lda(X,Y,tol = 1.0e-4,method="moment")
  YX <- data.frame(Y, X)
  Accuracy.rate <- function(formula, data, indices){
    Modlda<-lda(formula, data, CV = TRUE)
    Class<-table(data$Y)
    ClassTest<-table(Modlda$class,data$Y)
    aux1<-diag(ClassTest)/Class
    aux2<-mean(aux1)
    error<-signif(aux2, 4)  
    return(error)
  } # MisClass.
  boot.out <- boot(data = YX, statistic = Accuracy.rate , R = 100, formula = Y~X)
  Error <- mean(boot.out$t)
  return(list(Fitness=Error,Coef=StoreLDA$scaling))
  # Salidas
  # Fitness: precisi?n del modelo descriminate
  # Coef   : modelo descriminate 
}