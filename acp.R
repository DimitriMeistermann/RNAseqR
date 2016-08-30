###############################################################################
# FONCTIONS POUR L'ACP                                                        #
###############################################################################

# ACP
# Paramètre :	Table R des données
# Sortie :	un objet de type ACP
ACP<-function(d,transpose=T,scale=F,center=T) {
	if(transpose) d<-t(d);
	resacp <-prcomp(x = d,retx = T,center = center,scale = scale);
	resacp$n.obs<-dim(d)[1];
	resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
	return(resacp);
}

	
# Fonctions présentes dans R pour interpréter une ACP :
# soient "acp" un objet de type acp et T la table des données
#
# coefficients de corrélation entre les variables cor(T)
# Inertie des axes, axe i : acpn$sdev^2 , acpn$sdev[i]^2
# Nouvelles coordonnées des individus : acpn$scores
# De l'individu i sur CP j : acpn$scores[i,j]
# Graphique des inerties (valeurs propres) : plot(acp)
# Plan principal : biplot(acp)
# Plan i x j : biplot(acp,c(c1=i, c2=j))


#Tableau des inerties et des pourcentages cumulés
# Paramètre : résultat ACP
# Sortie : tableau des inerties, pourcentages et pourcentages cumulés
VP <- function(resacp) {
	N <- length(resacp$sdev);
	tab <- matrix(nrow=N,ncol=3);
	s <- sum(resacp$sdev^2);
	s1 <- 0;
	for (i in 1:N) {
		tab[i,1] <- resacp$sdev[i]^2;
		tab[i,2] <- resacp$sdev[i]^2/s;
		s1 <- s1+resacp$sdev[i]^2;
		tab[i,3] <- 100*s1/s;
	};
	return(tab)}


# Corrélations entre les axes et les variables initiales
# Paramètres :	table R des données
#		résultat ACP (produit par princomp)
# Sortie :	la matrice des corrélations
AXEVAR <- function(resacp) {return(resacp$rotation)}


# Corrélations entre les k premiers axes et les variables initiales
# Paramètres :	table R des données
#		résultat ACP (produit par princomp)
#		nombre axes
# Sortie :	la matrice des corrélations
AXEVARk <- function(d,resacp,k) {
	return(resacp$rotation[,1:k])
}

# Contribution de la ligne i à l'inertie de l'axe j
# Paramètres :	résultat ACP
#		numéro ligne
#		numéro axe
# Sortie :	pourcentage de la contribution
CTRij <- function(resacp,i,j) {
	x <- resacp$rotation[i,j]^2/(resacp$n.obs * resacp$sdev[j]^2);
        x <- 100*x;
	return(x)}


# Tableau des contribution des lignes aux axes
# Paramètres :	résultat ACP
#		nombre axes
# Sortie :	tableau des pourcentages des contributions
CTR <- function(resacp, nbax) {
      	matrice <- matrix(nrow=resacp$n.obs,ncol=nbax);
        row.names(matrice) <- row.names(resacp$x);
        for (j in 1:nbax) 
		for (i in 1:resacp$n.obs) matrice[i,j] <- CTRij(resacp,i,j);
     
       return(matrice)}

# Fonction utilitaire
SOMME2 <- function(resacp) {
	N <- resacp$n.obs ; M <- ncol(resacp$x);
	s2 <- vector (mode = "numeric", N);
	for (i in 1:N)
		for (j in 1:M) s2[i] <- s2[i] + resacp$x[i,j]^2;
	return(s2)
}

# Cosinus ** 2 des angles de projection
# Paramètres :	résultat ACP
#		nombre axes
# Sortie :	tableau des cos2 des angles de projection
COS2TETA <- function(resacp, nbax) {
	N <- resacp$n.obs ; 
	c2teta <- matrix(nrow=N,ncol=nbax);
	row.names(c2teta) <- row.names(resacp$x);
	s2 <- SOMME2(resacp);
	for (i in 1:N)
		for (j in 1:nbax) c2teta[i,j] <- resacp$x[i,j]^2 / s2[i];
	return(c2teta)
}

############################################################
#PLOTS													   #
############################################################
# Raccourci pour faire afficher un plan de projection
# Paramètres :	résultat ACP
#		premier axe choisi
#		deuxième axe choisi	
PLAN <- function(resacp,i,j) {biplot(resacp,c(c1=i,c2=j))}


# visualisation 3d ggplot

acp3d<-function(pca, comp=1:3, group=rep(1,pca$n.obs), plotVars = FALSE, pointSize=2, plotText=FALSE){
	if(!require("rgl")) stop("You must install rgl");
	if(length(comp)!=3) stop("You must give a vector of 3 integer for comp parameter")
	if(!plotVars){
		x<-pca$x
	}else{
		x<-pca$rotation
	}
	if(is.null(levels(group))){ colors="black"}
	else{
		hashCol<-rainbow(nlevels(group))
		names(hashCol)<-levels(group)
		colors<-hashCol[group]
	}

	
	percentVar <- pca$percentVar
	plot3d(x[,comp[1]],x[,comp[2]],x[,comp[3]],
		xlab=paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance"), 
		ylab=paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance"), 
		zlab=paste0("PC",comp[3],": ",round(percentVar[comp[3]] * 100),"% variance"),
	col=colors,size=pointSize,type=ifelse(plotText,"n","p"))
	
	legend3d("topright", legend = names(hashCol), pch = 16, col = hashCol, cex=1, inset=c(0.02))
	
	if(plotText) text3d(x[,comp[1]],x[,comp[2]],x[,comp[3]],texts=rownames(x),cex=pointSize,col=colors)
	if(plotVars) spheres3d(x=0,y=0,z=0, radius = 1,alpha=0.5,color="white")
	spheres3d(x=0,y=0,z=0, radius = 0.005,alpha=1,color="red")
}

# visualisation 2d ggplot
acp2d<-function(pca, comp=1:2,group=rep(1,pca$n.obs), plotVars = FALSE, pointSize=2, plotText=FALSE){
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(comp)!=2) stop("You must give a vector of 2 integer for comp parameter");
	
	percentVar <- pca$percentVar
    
	if(plotText){
		functPlot=geom_text
	}else{
		functPlot=geom_point
	}
	
	if(!plotVars){
		d <- data.frame(PC1=pca$x[,comp[1]], PC2=pca$x[,comp[2]], group=group);
		
		ggplot(data=d, mapping = aes(x=PC1, y=PC2,colour=group, label = rownames(d))) + 
		functPlot(size=pointSize)+
		xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
		ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance")) +
		coord_fixed()
	}else{
		if(pca$n.obs==length(group)) group=rep(1,nrow(pca$rotation))
		d <- data.frame(PC1=pca$rotation[,comp[1]], PC2=pca$rotation[,comp[2]], group=group);
		
		ggplot(data=d, mapping = aes(x=PC1, y=PC2,colour=group, label = rownames(d))) + 
		functPlot(size=pointSize)+
		xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
		ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance")) +
		coord_fixed()
	}
}

###############################################################################
# Projections d'une partition (obtenue ici avec hclust) sur un plan factoriel #
# (ici : projection sur le plan principal)                                    #
###############################################################################

# Paramètres : résultat de l'acp (renvoyé par princomp),
#              résultat de la cah (renvoyé par hclust),
#              nombre de classes de la cah

# Sortie : les éléments des classes projetés sur le plan principal
#          une couleur par classes

CAHsurACP<-function(acp,cah,k){
    n<-acp$obs
    classe<-cutree(cah,k)
    couleur<-vector("numeric",n)
    # pour avoir la liste des couleurs, on peut raper sous R : colors()
    liste_coul<-c("blue","red", "green", "grey", "orange", "turquoise", "yellow")
    for(i in 1:n){
        couleur[i]<-liste_coul[classe[i]]
    }
    plot (acp$scores[,1],acp$scores[,2], col=couleur)
}
