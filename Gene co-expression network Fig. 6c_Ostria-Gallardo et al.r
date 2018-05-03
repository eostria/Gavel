## Gene co-expression network (Fig. 6c in Ostria-Gallardo et al)##

> library(WGCNA)

> options(stringsAsFactors=F)
> enableWGCNAThreads()
## Allowing parallel execution with up to 3 working processes.
> library(igraph)

## Read in data
> data<-read.csv("Gave.interestORF.light.csv",sep=",",quote="")
> row.names(data)<-make.names(data[,1],T)## Assigning column 1 as row names
> dataORF=as.data.frame(data[,-c(1)])## Eliminates column 1 from data frame
> head(dataORF)# Explore the first 7 rows of data frame
> dim(dataORF) #[1] 436  24
> id=rownames(dataORF)

## Select a soft threshold to get a scale-free topology fit index
> counts.t<-t(dataORF)
> powers=c(c(1:30),seq(from=12,to=30,by=2))

## Call the network topology analysis function
> sft=pickSoftThreshold(counts.t,powerVector=powers,RsquaredCut=0.85,verbose=5)

## Create TOM
> adjacency=adjacency(counts.t,power=17)
> k=as.vector(apply(adjacency,2,sum,na.rm=T))
> sizeGrWindow(10,5)
> par(mfrow=c(1,2))
> hist(k)
> scaleFreePlot(k,main="Scale free topology\n") ##  scaleFreeRsquared slope 1              0.72 -0.72
> dev.off()

## Calculation of the topological overlap matrix, and the corresponding dissimilarity, from a given adjacency matrix.
> TOM=TOMsimilarity(adjacency)
> colnames(TOM)=id
> rownames(TOM)=id

## Get hub genes
> hub=rowSums(TOM)

## igraph adjacency matrix
> net=graph.adjacency(TOM,mode="undirected",weighted=T,diag=F)
> summary(net) ## IGRAPH UNW- 436 94830 --  + attr: name (v/c), weight (e/n)

## Adjacency matrix for visualization
> subTOM=(TOM>0.05)*TOM
> subnet=graph.adjacency(subTOM,mode="undirected",weighted=T,diag=F)
> summary(subnet)

## Finding community structure by greedy optimization of modularity
> community.fastgreedy=fastgreedy.community(subnet)
> table(community.fastgreedy$membership)


## Visualization of the network
> V(subnet)$color<-"lightgray"
> V(subnet)[community.fastgreedy$membership==1]$color<-"greenyellow"
> V(subnet)[community.fastgreedy$membership==4]$color<-"gold1"
> V(subnet)[community.fastgreedy$membership==15]$color<-"firebrick1"
> V(subnet)[community.fastgreedy$membership==6]$color<-"darkorchid1"
> V(subnet)[community.fastgreedy$membership==7]$color<-"dodgerblue1"
> V(subnet)[hub>30]$color<-"red"
> v.label=rep("",length(V(subnet)))
> v.label=V(subnet)$name
> v.size=rep(3,length(V(subnet)))
> V(subnet)$shape<-"circle"
> pdf("GaveNetwork.pdf",useDingbats=F)
> plot(subnet,layout=layout.auto,vertex.size=v.size,vertex.frame.color=NA,vertex.label=v.label,vertex.label.cex=0.05,edge.color="snow3",edge.width=E(subnet)$weight)
> dev.off()
