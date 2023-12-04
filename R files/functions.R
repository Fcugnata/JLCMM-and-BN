#############################################################################
####plotbn.R
####Plot a Bayesian network as a graph whose nodes are barplots representing 
####the marginal distributions of the corresponding variables.
#############################################################################

plotbn<-function(data,x,v_class,cex_p=2){

fit<-bn.fit(x, data)

qd<-gRain::querygrain(as.grain(fit), nodes = colnames(data), 
        type = "marginal")

# convert arc to edge 
mat<-matrix(0,length(bnlearn::nodes(x)),length(bnlearn::nodes(x)))
rownames(mat) <- bnlearn::nodes(x)
colnames(mat) <- bnlearn::nodes(x)
for (i in 1:dim(arcs(x))[1]){
mat[bnlearn::nodes(x)==arcs(x)[i,1],bnlearn::nodes(x)==arcs(x)[i,2]]<-1
}
# create the graphAM object from the bn object
g1 <- graphAM(adjMat=mat,edgemode="directed")
g1layout <- agopen(g1, name="foo")
plot(g1layout)

dd<-list()
c2<-colnames(data)
for (i in c2){
dd[[i]]<-makeNodeDrawFunction(qd[[i]],pp=0.8)
if (i%in%v_class){
dd[[i]]<-makeNodeDrawFunction_class(qd[[i]],pp=0.8)
}
}

attrs <- list(node=list(height=2,width=3.8,shape="box", fixedsize=TRUE),edge=list(arrowsize=0.5))
g1layout <- agopen(g1, name="foo",attrs=attrs)
plot(g1layout, drawNode=dd)
}


##########################################################################
####whatIFplot.R
####Produce a BN plot with highlighted the evidence and of consequent 
####conditional probabilities of the target nodes obtained with
####the querygrain function in the package gRain.
##########################################################################


whatIFplot<-function(data,x,evidenceName, states,cex_p=2){
if (!is.list(states)){
states <- as.list(states)}
fit<-bn.fit(x, data)
jtree = compile(as.grain(fit))
evidence<-which(colnames(data)%in%evidenceName)
kk<-colnames(data)[!colnames(data)%in%evidenceName]

jprop<-setEvidence(jtree, nodes = evidenceName, states = states)
qd<-querygrain(jprop)

# convert arc to edge 
mat<-matrix(0,length(bnlearn::nodes(x)),length(bnlearn::nodes(x)))
rownames(mat) <- bnlearn::nodes(x)
colnames(mat) <- bnlearn::nodes(x)
for (i in 1:dim(arcs(x))[1]){
mat[bnlearn::nodes(x)==arcs(x)[i,1],bnlearn::nodes(x)==arcs(x)[i,2]]<-1
}
# create the graphAM object from the bn object
g1 <- graphAM(adjMat=mat,edgemode="directed")
g1layout <- agopen(g1, name="foo")
plot(g1layout)

dd<-list()
c2<-which(colnames(data)%in%kk)
for (i in c2){
dd[[i]]<-makeNodeDrawFunction(qd[[colnames(data)[i]]],pp=0.8)
if (i%in%c(13,14)){
dd[[i]]<-makeNodeDrawFunction_class(qd[[colnames(data)[i]]],pp=0.8)
}

}


j_e=0
for (i in evidenceName){
j_e=j_e+1
i_ev<-which(colnames(data)==i)
jj<-length(states[[j_e]])
if (jj!=1){
l<-levels(data[,i_ev])
count<-states[[j_e]]
names(count)<-l
}
else {
l<-levels(data[,i_ev])
count<-rep(0,length(l))
names(count)<-l
count[states[[j_e]]]<-100
}
dd[[i_ev]]<-makeNodeDrawFunction_ev(count,pp=0.8)
}

attrs <- list(node=list(height=2,width=3.8,shape="box", fixedsize=TRUE),edge=list(arrowsize=0.5))
g1layout <- agopen(g1, name="foo",attrs=attrs)
plot(g1layout, drawNode=dd)
}


##########################################################################
####Other functions
##########################################################################

barGlyph<-function (x, xposc, yposc, labels = names(x), main = NULL, height = 10, width = 16,
    density = NULL, angle = 45, col=rep("lightgoldenrod1",5),col2="black",col3="white", cex1=2,cex2=2,cexT=2,border = NULL, lty = NULL, ...) 
{
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("pie: `x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(1:length(x))
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
xpos<-xposc-width/2
ypos<-yposc-height/2
vp<-dx*100
vp<-round(vp, digits = 1)
vp<-paste(vp," ",sep = "")
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("lightblue", "mistyrose", "lightcyan", "lavender", 
                "cornsilk", "white")
        else par("fg")
    if (!is.null(col)) 
        col <- rep(col, length.out = nx) 
    if (!is.null(border)) 
        border <- rep(border, length.out = nx)
    if (!is.null(lty))
        lty <- rep(lty, length.out = nx)
    if (!is.null(angle)) 
        angle <- rep(angle, length.out = nx)
    if (!is.null(density)) 
        density <- rep(density, length.out = nx)

    p1<-width/10
    h1<-height/4
    k<-(height-h1)/nx
    s<-k/6 
    kk<-(k/5*4)-(s/nx)
    rect(xpos,ypos,xpos+width,ypos+height-h1,border=col2,col=col3)
    xc <- xpos+(p1*8)
    segments(xc, ypos, x1 = xc, y1 = ypos+height-h1)
      for (i in 1:nx) {
        #n <- max(2, floor(edges * dx[i]))
        #t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
        yc <- ypos+((i-1)*kk)+s*i
        ll<-dx[i]*(width-(p1*8))  
        rect(xc, yc, xc+ll, yc+kk,density = density[i], col = col[i])
        #t2p <- 2 * pi * mean(x[i + 0:1])
        xt <-xpos-p1/2
        yt <- (yc+yc+k)/2
        text(xt, yt-p1/5, labels[i],  cex = 1,pos=4) 
        xt2<-xpos+(p1*9.1) 
        text(xt2, yt-p1/5, vp[i],  cex = 1.2,pos=2)
        xtit<-xpos+width/2
        ytit<-ypos+3.5*h1
        text(xtit, ytit, main,  cex = 1.5)
    }
    invisible(NULL)
}

makeNodeDrawFunction <- function(x,pp=2) {
force(x)
function(node, ur, attrs, radConv) {
nc <- getNodeCenter(node)
na <- name(node)
hh<-getNodeHeight(node)
ww<- getNodeLW(node) + getNodeRW(node)
barGlyph(x,xpos=getX(nc),ypos=getY(nc),col=rep("lightgoldenrod1",5),height=150, width=270,main=na,cexT=pp)
}
}

makeNodeDrawFunction_class <- function(x,pp=2) {
force(x)
function(node, ur, attrs, radConv) {
nc <- getNodeCenter(node)
na <- name(node)
hh<-getNodeHeight(node)
ww<- getNodeLW(node) + getNodeRW(node)
barGlyph(x,xpos=getX(nc),ypos=getY(nc),col=c("darkgreen","blue","darkred"),height=150, width=270,main=na,cexT=pp)
}
}

makeNodeDrawFunction_ev <- function(x,pp=2) {
force(x)
function(node, ur, attrs, radConv) {
nc <- getNodeCenter(node)
na <- name(node)
hh<-getNodeHeight(node)
ww<- getNodeLW(node) + getNodeRW(node)
barGlyph(x,xpos=getX(nc),ypos=getY(nc),col=c("lightgoldenrod1"),height=150, width=270,main=na,col3="gray87",cexT=pp)
}
}



