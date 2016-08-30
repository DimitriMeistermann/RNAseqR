# Calculs des matrices de distances par analyse de rang

revDist<-function(vectA,vectB){
	if(length(vectB)>length(vectA)){
		vectT<-vectA;
		vectA<-vectB;
		vectB<-vectT;
		rm(vectT);
	}
	compVectB<-1:length(vectB);
	names(compVectB)<-vectB;
	scoreT<-0;
	
	for(i in 1:length(vectA)){
		score<-abs(compVectB[vectA[i]]-1);
		scoreT<-scoreT+score;
	}
	return(scoreT);
}

compareVectReplace<-function(data){
	#return a dist matrix
	#@param data : dataframe or matrx (rows = genes, cols = samples)
	# tri alph
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix("",nrows,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		listOfgenesBySample[,i]<-names(sort(col,decreasing = T,method="shell"));
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	tot<-length(as.dist(distMat))
	c<-1;
	comptT<-0
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-revDist(listOfgenesBySample[,i],listOfgenesBySample[,j]);
			comptT<-comptT+1;
			print(paste0(comptT,"/",tot));
		}
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}