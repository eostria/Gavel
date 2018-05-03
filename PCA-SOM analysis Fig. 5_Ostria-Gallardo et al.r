### PCA-SOM analysis Fig. 5_Ostria-Gallardo et al.###


###############PCA###############
#################################

## Read in data
> data<-read.csv("Gavel.light.ORF.ave.csv", sep=",",quote="",header=T)

## Explore and write a new csv file from data
> head(data)
> colnames(data)<-c("SeqID","DS","SH","SS","SU","sd","average","cv")
> write.csv(data,"Gavel.light.2017.csv")
> data<-read.csv("Gavel.light.2017.csv",sep=",",quote="",header=T)
> head(data)


## Choose columns with count No, average, standar deviation (sd), and coefficient of variation (cv)
> data.cv<-subset(data[2:9])

## Check data.cv
> head(data.cv)
> data.cv[is.na(data.cv)]<-0

## Let's give better column names
> colnames(data.cv)<-c("SeqID","DS","SH","SS","SU","sd","average","cv")

## Select only those genes upper 75% CV
> quantile(data.cv$cv)
> DataCV<-subset(data.cv,cv>0.424470783)
> head(DataCV)
> length(DataCV$SeqID) #[1] 11101

## Create a matrix of the data to perform a PCA on and scale it
> library (ggplot2)
> m.DataCV<-as.matrix(DataCV[2:5])
> sc.DataCV<-t(scale(t(m.DataCV)))
> tisdata<-as.matrix(sc.DataCV,dimnames=list(rownames(X)))
> tispca<-prcomp(tisdata,scale=T)
> summary(tispca)
> tis.pca.scores<-data.frame(tispca$x)

## Write out a master data file with origib¿nal data, scaled data, and PCA results
> allData<-cbind(DataCV,sc.DataCV,tis.pca.scores)

## Let's give better column names "sc" denotes scaled gene expression rather than the original value 
> colnames(allData)<-c("SeqID","DS","SH","SS","SU","sd","average","cv","sc.DS","sc.SH","sc.SS","sc.SU","PC1","PC2","PC3","PC4")
> head(allData)
#           SeqID         DS         SH          SS         SU         sd
#1  """m.54690"""  0.6275563  0.8352066   2.9973896  64.259159 31.4045077
#3  """m.61325"""  0.6884803  2.7874280   0.6826512   4.112466  1.6851985
#22 """m.42375"""  0.7575304  0.5475163   1.7028975   1.698551  0.6112249
#24 """m.59733"""  5.9022665 34.5429136  55.2533270  51.996075 22.5916591
#30  """m.8359""" 46.8945175 53.0970590 145.6338833 113.164306 47.7886232
#38 """m.32224"""  0.9904402  3.1360189   8.0172086   3.551932  2.9506767
#     average        cv      sc.DS      sc.SH      sc.SS      sc.SU         PC1
#1  17.179828 1.8279873 -0.5270667 -0.5204546 -0.4516052  1.4991265 -1.40256513
#3   2.067756 0.8149890 -0.8184650  0.4270546 -0.8219240  1.2133344 -0.76688543
#22  1.176624 0.5194735 -0.6856616 -1.0292570  0.8610150  0.8539036 -0.65771340
#24 36.923645 0.6118480 -1.3731342 -0.1053810  0.8113473  0.6671679 -0.02426466
#30 89.697441 0.5327758 -0.8956718 -0.7658807  1.1704970  0.4910555 -0.09015822
#38  3.923900 0.7519755 -0.9941651 -0.2670170  1.3872441 -0.1260620  0.80377897
#           PC2        PC3           PC4
#1  -0.05376640 -0.1635600  0.000000e+00
#3   0.07217832 -1.4554605  8.326673e-16
#22  0.94906095  1.5646790 -2.386980e-15
#24  1.75762788  0.4344667 -3.691492e-15
#30  1.39558251  1.5823864 -3.226586e-15
#38  1.64748728  1.3714973 -3.330669e-15

## Write out the PC scores and other data.
> write.table(allData,file="pca_scores_Gavel2017.txt")
> write.table(tispca$rotation, "loadings_Gavel2017.txt")

## Let's visualize the PCA space by plotting the different PCs against each other
> library(ggplot2)
> p<-ggplot(allData,aes(PC1,PC2))
> p+geom_point(alpha=0.1)+theme_bw()
> p<-ggplot(allData,aes(PC2,PC3))
> p+geom_point(alpha=0.1)+theme_bw()
> p<-ggplot(allData,aes(PC1,PC3))
> p+geom_point(alpha=0.1)+theme_bw()
> p<-ggplot(allData,aes(PC2,PC4))
> p+geom_point(alpha=0.1)+theme_bw()
> p<-ggplot(allData,aes(PC1,PC4))
> p<-ggplot(allData,aes(PC1,PC4))
> p+geom_point(alpha=0.1)+theme_bw()
> p<-ggplot(allData,aes(PC3,PC4))
> p+geom_point(alpha=0.1)+theme_bw()

> ###############SOM###############
> #################################

##Let´s do a self-organizing map (SOM). Read in the Kohonen package
> library(kohonen)

## Create a matrix using scaled gene expression values, for only the columns you want to perform a SOM
> sc.sig<-as.matrix(allData[9:12])

## Set a random seed. SOMs use a random seed, but setting it beforehand will make your results reproducible.
> set.seed(2)

## Perform a SOM. Please think and check the topology and number of nodes. Spend some time trying different combinations and visualize the results to get an appropriate SOM node number  
> ssom<-som(sc.sig,somgrid(2,3,"hexagonal"))
> summary(ssom)

## Look at the training of the SOM ("changes"), the types of patterns represented in each node ("codes"), the total number of assigned genes in each node ("counts"), and the average distance of each gene from the node pattern ("quality")
> plot(ssom,type="changes")
> plot(ssom,type="codes")
> plot(ssom,type="counts")
> plot(ssom,type="quality")


## Get the node "codes"
> codes<-ssom$codes
> write.table(codes,file="Gave2012_codes.txt")


## Combine the PCA information with the SOM node assignments and the distance of each gene from its respective SOM node
> data.som<-cbind(allData,ssom$unit.classif,ssom$distances)

## Rename the clumn heading better
> colnames(data.som)<-c("SeqID","DS","SH","SS","SU","sd","average","cv","sc.DS","sc.SH","sc.SS","sc.SU","PC1","PC2","PC3","PC4","som.node","distance")
> head(data.som)
#           SeqID         DS         SH          SS         SU         sd
#1  """m.54690"""  0.6275563  0.8352066   2.9973896  64.259159 31.4045077
#3  """m.61325"""  0.6884803  2.7874280   0.6826512   4.112466  1.6851985
#22 """m.42375"""  0.7575304  0.5475163   1.7028975   1.698551  0.6112249
#24 """m.59733"""  5.9022665 34.5429136  55.2533270  51.996075 22.5916591
#30  """m.8359""" 46.8945175 53.0970590 145.6338833 113.164306 47.7886232
#38 """m.32224"""  0.9904402  3.1360189   8.0172086   3.551932  2.9506767
#     average        cv      sc.DS      sc.SH      sc.SS      sc.SU         PC1
#1  17.179828 1.8279873 -0.5270667 -0.5204546 -0.4516052  1.4991265 -1.40256513
#3   2.067756 0.8149890 -0.8184650  0.4270546 -0.8219240  1.2133344 -0.76688543
#22  1.176624 0.5194735 -0.6856616 -1.0292570  0.8610150  0.8539036 -0.65771340
#24 36.923645 0.6118480 -1.3731342 -0.1053810  0.8113473  0.6671679 -0.02426466
#30 89.697441 0.5327758 -0.8956718 -0.7658807  1.1704970  0.4910555 -0.09015822
#38  3.923900 0.7519755 -0.9941651 -0.2670170  1.3872441 -0.1260620  0.80377897
#           PC2        PC3           PC4 som.node    distance
#1  -0.05376640 -0.1635600  0.000000e+00        6 0.007985779
#3   0.07217832 -1.4554605  8.326673e-16        6 1.076692739
#22  0.94906095  1.5646790 -2.386980e-15        4 0.386934792
#24  1.75762788  0.4344667 -3.691492e-15        4 0.423236729
#30  1.39558251  1.5823864 -3.226586e-15        4 0.519489755
#38  1.64748728  1.3714973 -3.330669e-15        3 0.473594194

## Write out the SOM data
> write.table(data.som,file="Gave2017_som.data.txt")

## Visualize the SOM node membership on the PCA space.
> p<-ggplot(data.som,aes(PC1,PC2,colour=factor(som.node)))
> p+geom_point(alpha=0.2)+theme_bw()+scale_colour_manual(values=c("red","blue","green","gray25","purple","gold"))

> p<-ggplot(data.som,aes(PC2,PC3,colour=factor(som.node)))
> p+geom_point(alpha=0.2)+theme_bw()+scale_colour_manual(values=c("red","blue","green","gray25","purple","gold"))

> p<-ggplot(data.som,aes(PC1,PC3,colour=factor(som.node)))
> p+geom_point(alpha=0.2)+theme_bw()+scale_colour_manual(values=c("red","blue","green","gray25","purple","gold"))

## Generate box plots of each node
> library(reshape)


## Isolate just the seqIDs, the nodes, and the scaled gene expression
> box.data<-data.som[c(1,9:12,17)]


## Melt the data to reformat it
> melt.data<-melt(box.data,id=c("SeqID","som.node"))

## Finally, let´s visualize the expression patterns of node members. Use the subset function to set the node you want to look
> node<-subset(melt.data,som.node==4)
> p<-ggplot(node,aes(x=variable,y=value))
> p+geom_point(colour="black",position="jitter",alpha=0.1)+geom_boxplot(outlier.size=0,size=2,alpha=0.6)+theme_bw()

> node<-subset(melt.data,som.node==1)
> p<-ggplot(node,aes(x=variable,y=value))
> p+geom_point(colour="black",position="jitter",alpha=0.1)+geom_boxplot(outlier.size=0,size=2,alpha=0.6)+theme_bw()

> node<-subset(melt.data,som.node==2)
> p<-ggplot(node,aes(x=variable,y=value))
> p+geom_point(colour="black",position="jitter",alpha=0.1)+geom_boxplot(outlier.size=0,size=2,alpha=0.6)+theme_bw()

> node<-subset(melt.data,som.node==3)
> p<-ggplot(node,aes(x=variable,y=value))
> p+geom_point(colour="black",position="jitter",alpha=0.1)+geom_boxplot(outlier.size=0,size=2,alpha=0.6)+theme_bw()

> node<-subset(melt.data,som.node==5)
> p<-ggplot(node,aes(x=variable,y=value))
> p+geom_point(colour="black",position="jitter",alpha=0.1)+geom_boxplot(outlier.size=0,size=2,alpha=0.6)+theme_bw()

> node<-subset(melt.data,som.node==6)
> p<-ggplot(node,aes(x=variable,y=value))
> p+geom_point(colour="black",position="jitter",alpha=0.1)+geom_boxplot(outlier.size=0,size=2,alpha=0.6)+theme_bw()

################################################################################################################################## 
