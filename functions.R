#' Preprocess 10X
#'
#' 
#' This function takes a data-frame (genes x cells), 
#' creates a Seurat object with it and filters the object for default tags such as  
#' Min and max nFeature_RNA and % of MT
#' 
#'
#' @param data data-frame
#' @param percent_mt integer [0-100]
#' @param max_features integer [0-Inf]
#' @param min_features integer [0-Inf]
#' @return Preprocessed Seurat object
#' @export
preprocess_10X<-function(data,name='10X-project',percent_mt=20, max_features=5000, min_features=200){

	parameters = list(
			min = min_features,
			max = max_features,
			mt = percent_mt
		)
	seurat.obj = CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
	seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = '^mt-|^MT-')
	seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt<20)
	return(seurat.obj)
}
#' Load 10X
#'
#' 
#' This function takes the path to a 10X output folder and instanciates the Seurat object
#'
#' @param file string (path to file)
#' @param percent_mt integer [0-100]
#' @param max_features integer [0-Inf]
#' @param min_features integer [0-Inf]
#' @return Preprocessed Seurat object
#' @export
load_10X<-function(file,name='10X-project',percent_mt=20, max_features=5000, min_features=200){
	data = Read10X(data.dir = file)
	seurat.obj = preprocess_10X(data,percent_mt=percent_mt,max_features=max_features,min_features=min_features)
	return(seurat.obj)
}
#' Reduce Dimensions
#'
#' 
#' This function reduces the dimensions of the Normalized Seurat object
#' Runs PCA and then UMAP and then performs clustering 
#'
#' @param seurat.obj S4 instance
#' @param ndims integer [3-100]
#' @param res double [0-3]
#' @return Seurat object with reduction embeddings
#' @export
reduce_dim<-function(seurat.obj,ndims=15,res=.1){
	seurat.obj = SCTransform(seurat.obj)
	seurat.obj <- RunPCA(seurat.obj, verbose = TRUE)
	seurat.obj <- RunUMAP(seurat.obj, dims = 1:ndims)
	seurat.obj <- FindNeighbors(seurat.obj, dims = 1:ndims)
	seurat.obj <- FindClusters(seurat.obj, resolution = res)
	return(seurat.obj)	
}
#' Integrate batches
#'
#' 
#' Use the SCT approach of Seurat to integrate batches together 
#' 
#'
#' @param obj.list list of the S$ objects to be merged
#' @return Integrated Seurat object
#' @export
integrate_batches<-function(obj.list){

	features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 5000)
	obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features, verbose = TRUE)
	anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", 
	    anchor.features = features, verbose = TRUE)
	integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
	    verbose = TRUE)
	integrated<-reduce_dim(integrated)
	return(integrated)
	
}
#' Merge samples 
#'
#' 
#' This function allows to instantiate a merged Seurat object of several distinct samples
#' 
#'
#' @param files character vector of the files to merge 
#' @return Merged Seurat object 
#' @export
merge_samples<-function(files){

	to.merge=list()
	for(f in names(files)){
		to.merge[[f]]=load_10X(files[[f]],name=f)
	}	
	to.merge = merge(to.merge[[1]],to.merge[2:length(to.merge)])
	to.merge = SCTransform(to.merge)
	to.merge = reduce_dim(to.merge)
	return(to.merge)
}
#' Recluster 
#'
#' 
#' Takes a subset of a Seurat object and performs reclustering by instanciating a new object  
#' 
#'
#' @param sub.obj Seurat object 
#' @return Processed reclustering
#' @export
recluster<-function(sub.obj,name='10X-subcluster',batch=NA){

	if(is.na(batch)){
		new = preprocess_10X(sub.obj[['RNA']]@data)
		new = SCTransform(new)
		new = reduce_dim(new)
	}else{
		obj.list=list()
		for(name in names(sub.obj)){
			obj.list[[name]]=preprocess_10X(sub.obj[['RNA']]@data)
			obj.list[[name]]=SCTransform(obj.list[[name]])
		}
		new = integrate_batches(obj.list)
	}
	return(new)
	
}
#' Run PHATE 
#'
#' 
#' Run PHATE and MAGIC on a Seurat object. If the integrated flag if up, the data extracted 
#' is the normalized and integrated data from seurat and no further tranformation is performed
#' prior to running PHATE. Otherwise raw data is extracted and tranformation is performed
#' 
#' @param seurat.obj Seurat object 
#' @param integrated Boolean
#' @param feature Str (Feature; Gene to plot)
#' @param t Double - PHATE hyperparameter 
#' @param subset Str - Subset a specific group of cell from the Seurat Object 
#' @return list 
#' @export
run_phate<-function(seurat.obj, integrated=F, feature=NA, t=NA, subset=NA){
	
	if(!is.na(subset)){
		seurat.obj=subset(seurat.obj, idents=subset)
	}

	if(integrated){
		data.phate=t(data.frame(seurat.obj[['integrated']]@scale.data))
	}
	else{
		data.phate=t(data.frame(seurat.obj[['RNA']]@data))
		keep_cols <- colSums(data.phate > 0) > 10
		data.phate <- data.phate[,keep_cols]
		keep_rows <- rowSums(data.phate) > 1000
		data.phate <- data.phate[keep_rows,]
		data.phate <- library.size.normalize(data.phate)
		data.phate <- sqrt(data.phate)
	}
	if(is.na(t)){
		data.phate_PHATE <- phate(data.phate)
	}
	else{
		data.phate_PHATE <- phate(data.phate, t=t)
	}
	
	if(!integrated){
		data.magic= magic(data.phate, t=4)
	}
	if(is.na(feature)){
		p=ggplot(data.phate_PHATE) + geom_point(aes(PHATE1, PHATE2))
	}
	else{
		p=ggplot(data.phate_PHATE) + geom_point(aes(PHATE1, PHATE2, col=data.magic$result[,feature]))+scale_color_viridis(option='B')
	}
	return(list(plot=p, data.phate=data.phate_PHATE, data.magic=data.magic))
}
#' Run Monocle trajectory 
#'
#' 
#' This function allows one to perform a monocle-pseudotime analysis
#' on a Seurat Object
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param batch Str - Column name of the batch
#' @return list
#' @export
run_trajectory<-function(seurat.obj, batch=NA){
	#Extract data, phenotype data, and feature data from the SeuratObject
	data <- as(as.matrix(seurat.obj@assays$RNA@data), 'sparseMatrix')
	
	pd <- new('AnnotatedDataFrame', data = seurat.obj@meta.data)
	
	fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
	fd <- new('AnnotatedDataFrame', data = fData)
	
	#Construct monocle cds
	HSMM <- newCellDataSet(data,
	                         phenoData = pd,
	                         featureData = fd,
	                         lowerDetectionLimit = 0.5,
	                         expressionFamily = negbinomial.size())
	
	HSMM <- estimateSizeFactors(HSMM)
	HSMM <- estimateDispersions(HSMM)
	HSMM <- detectGenes(HSMM, min_expr = 0.1)
	if(!is.na(batch)){
		HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
		            reduction_method = 'tSNE',
		            residualModelFormulaStr = paste0("~",batch," + num_genes_expressed",
		            verbose = T)
	}
	disp_table <- dispersionTable(HSMM)
	ordering_genes <- subset(disp_table,
	                  mean_expression >= 0.5 &
	                  dispersion_empirical >= 0.5 * dispersion_fit)$gene_id
	HSMM <- setOrderingFilter(HSMM, ordering_genes)
	HSMM <- reduceDimension(HSMM, max_components = 2,
	    method = 'DDRTree')
	HSMM <- orderCells(HSMM)
	p=plot_cell_trajectory(HSMM, color_by = "State")
	return(list(monocle_object=HSMM, trajectory=p))
}
#' Single-cell heatmap
#'
#' 
#' Make a heatmap of average expression from a seurat object 
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param group_name Str Name of the metadata feature to group 
#' @param group_name Str Name of the metadata feature to aggregate
#' @return list
#' @export
sc_heatmap<-function(seurat.obj, features){

}
#' Proportion of clusters
#'
#' 
#' Aggregate clusters of celltypes and extract proportion of a group
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param group_name Str Name of the metadata feature to group 
#' @param group_name Str Name of the metadata feature to aggregate
#' @return list
#' @export
sc_aggregate<-function(seurat.obj,group_name='orig.ident',by='seurat_clusters',levs=NA){
	df=data.frame(matrix(ncol=length(unique(seurat.obj@meta.data[,group_name]))))
	colnames(df)=unique(seurat.obj@meta.data[,group_name])
	first=T
	clusters=c()
	for(c in unique(seurat.obj@meta.data[,by])){
		t=table(seurat.obj@meta.data[which(seurat.obj@meta.data[,by]==c),group_name])/sum(table(seurat.obj@meta.data[which(seurat.obj@meta.data[,by]==c),group_name]))
		clusters=c(clusters,c)
		absent=colnames(df)[!colnames(df) %in% names(t)]
		if(length(absent)>0){
			for(a in absent){
				t[a]=0
			}
		}
		df=rbind(df,t[colnames(df)])
	}
	df=na.omit(df)
	df$clusters=clusters
	m=melt(df, id.vars='clusters')
	m$clusters=paste0(m$clusters)
	id1=unique(m$variable)[1]
	m.sub=subset(m,variable==id1)
	order=m.sub[order(m.sub$value),'clusters']
	m$clusters=factor(m$clusters,levels=order)
	if(!is.na(levs[1])){
		m$variable=factor(m$variable,levels=levs)
	}
	g=ggplot(m,aes(clusters,value, fill=variable))+geom_bar(stat='identity',col='black')+theme(axis.text.x=element_text(angle = 60, hjust = 1))
	return(list(a=g,b=m,c=df))
}
#' Proportion of clusters
#'
#' 
#' Aggregate clusters of celltypes and extract proportion of a group
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param tcr_folder Str Path to the tcr folder 
#' @param tag Str tag the columns on the metadata
#' @return Seurat object with clonotype annotation
#' @export
add_clonotype <- function(seurat.obj,tcr_folder,tag='clono'){
	
	tcr <- read.csv(paste(tcr_folder,"/filtered_contig_annotations.csv", sep=""))

	# Remove the -1 at the end of each barcode.
	# Subsets so only the first line of each barcode is kept,
	# as each entry for given barcode will have same clonotype.
	tcr$barcode <- gsub("-1", "", tcr$barcode)
	tcr <- tcr[!duplicated(tcr$barcode), ]

	# Only keep the barcode and clonotype columns. 
	# We'll get additional clonotype info from the clonotype table.
	tcr <- tcr[,c("barcode", "raw_clonotype_id")]
	names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

	# Clonotype-centric info.
	clono <- read.csv(paste(tcr_folder,"/clonotypes.csv", sep=""))

	# Slap the AA sequences onto our original table by clonotype_id.
	tcr <- merge(tcr, clono[, c(paste0("clonotype_id"), paste0("cdr3s_aa"))])

	# Reorder so barcodes are first column and set them as rownames.
	tcr <- tcr[, c(2,1,3)]
	rownames(tcr) <- tcr[,1]
	tcr[,1] <- NULL
	colnames(tcr)[colnames(tcr)=='clonotype_id']=paste0(colnames(tcr)[colnames(tcr)=='clonotype_id'],'.',tag)
	colnames(tcr)[colnames(tcr)=='cdr3s_aa']=paste0(colnames(tcr)[colnames(tcr)=='cdr3s_aa'],'.',tag)

	# Add to the Seurat object's metadata.
	clono_seurat <- AddMetaData(object=seurat.obj, metadata=tcr)
	return(clono_seurat)
}
#' Cell-Cycle regression
#'
#' 
#' Annotate the cells with their estimated cell cycle 
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param cc.filename Str path to cell cycle genes
#' @param species Str h.sapens or m.musculus
#' @param regress Boolean 
#' @return Seurat object with cell-cycle annotation or alleviated of cell-cycle effect 
#' @export
cellCycle_analysis<-function(seurat.obj, cc.filename, species='h.sapiens'){

	cc = read.table(cc.filename, header=T)
	if(species == 'h.sapiens'){
		feat = 'Gene_Name'
	}else{
		feat = 'mGeneName'
	}
	g2m.genes = subset(cc, Stage=='G2-M')[,feat]
	s.genes = subset(cc, Stage=='G1-S')[,feat]
	seurat.obj <- CellCycleScoring(seurat.obj, 
		s.features = s.genes, 
		g2m.features = g2m.genes, 
		set.ident = TRUE)
	return(seurat.obj)
}
#' Cell-Cycle regression
#'
#' 
#' Annotate the cells with their estimated cell cycle 
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param cc.filename Str path to cell cycle genes
#' @param species Str h.sapens or m.musculus
#' @param regress Boolean 
#' @return Seurat object with cell-cycle annotation or alleviated of cell-cycle effect 
#' @export
regress_out<-function(seurat.obj, vars){
	seurat.obj = SCTransform(seurat.obj, vars.to.regress = vars)
	seurat.obj = reduce_dim(seurat.obj)
	return(seurat.obj)
}

#' Write Garnett
#'
#' 
#' Write a new classifier file. Provide the annotation as follows:
#'
#'			 example = list(
#'				'T cells' = list(
#'						'expressed' = c('Cd3d','Cd3g','Cd3e'),
#'						subtypes = NA),
#'				'CD4' = list(
#'						'expressed' = c('Cd3d','Cd3g','Cd3e'),
#'						subtypes='T cells'),
#'				'CD8' = list(
#'						'expressed' = c('Cd3d','Cd3g','Cd3e'),
#'						subtypes='T cells'),
#'				'B cells' = list(
#'						'expressed' = c('Cd19','Ms4a1'),
#'						subtypes = NA	
#'						)
#'				)
#' 
#' 
#' @param dict List Annotation 
#' @param filename Str Output filename
#' @export
write_garnett<-function(dict,  filename){

	vec = c()
	for(n in names(dict)){
		tmp = paste0('>', paste(n ,'expressed: ' , sep='\n'), paste(dict[[n]][['expressed']], collapse=', '))
		if(!is.na(dict[[n]][['subtypes']])){
			tmp=paste0(paste(tmp,'subtype of: ',sep='\n'),paste(dict[[n]][['subtypes']],'','\n\n'))
		}else{
			tmp=paste(tmp,'',sep='\n\n')
		}
		vec=c(vec,tmp)

	}
	write.table(vec,filename,quote=F,row.names=F,col.names=F)

}
#' Train Garnett
#'
#' 
#' Train a garnett classifier based on an annotation file and a seurat obj
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param marker.path Str path to clasifier file
#' @param species Str h.sapens or m.musculus
#' @param unknown Integer: Number of unknown allowed  
#' @return Trained classifier 
#' @export
train_garnett<-function(seurat.obj,marker.path,species='h.sapiens', unknown = 50){

	cds = as.CellDataSet(seurat.obj)
	cds <- estimateSizeFactors(cds)

	if(species == 'h.sapiens'){
		gene.db = org.Hs.eg.db
	}else{
		gene.db = org.Mm.eg.db
	}

	marker_check <- check_markers(cds, marker.path,
                              db=gene.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

	marker.plot = plot_markers(marker_check)
	

	classifier <- train_cell_classifier(cds = cds,
                                         marker_file = marker.path,
                                         db=gene.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = unknown,
                                         marker_file_gene_id_type = "SYMBOL")
	return(classifier)
}
#' Run Garnett
#'
#' 
#' Run a garnett classifier based on a seurat obj
#' 
#' 
#' @param seurat.obj Seurat object 
#' @param classifier Trained classifier
#' @param species Str h.sapens or m.musculus
#' @return Annotation of cells
#' @export
run_garnett<-function(seurat.obj,classifier,species='h.sapiens'){
	
	cds = as.CellDataSet(seurat.obj)
	cds <- estimateSizeFactors(cds)
	
	if(species == 'h.sapiens'){
		gene.db = org.Hs.eg.db
	}else{
		gene.db = org.Mm.eg.db
	}
	cds <- classify_cells(cds, classifier,
                           db = gene.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")
	return(pData(cds))

}



obj = readRDS('pbmc.nov2019.rds')
Idents(obj)='true.ident'
u = unique(as.character(Idents(obj)))
u = u[!u %in% c('AP-1104.fresh','AP-1104.2days','MM-066.fresh')]
sub = recluster(subset(obj, idents = u),name='pbmc_good')
saveRDS(sub, 'pbmc_good.rds')
Idents(sub)='seurat_clusters'

p = FeaturePlot(s, features=c('CD3D',
						'CD4',
						'CD8A',
						'NKG7',
						'MS4A1',
						'S100A8',
						'LYZ',
						'CD14',
						'FCGR3A',
						'CLEC10A',
						'CLEC9A',
						'PPBP',
						'CCR7',
						'S100A4','nFeature_RNA'),order=T,combine=F,cols=rev(brewer.pal(11,'RdYlBu')))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

pl=cowplot::plot_grid(plotlist = p)
ggplot2::ggsave(plot = pl, 'featurePlot.pbmc_good.png')
d = DimPlot(s, label=T, label.size=6, group.by='orig.ident')+NoLegend()+NoAxes()
ggplot2::ggsave(plot = d, 'dimplot.orig.pbmc_good.png')
d = DimPlot(s, label=T, label.size=6)+NoLegend()+NoAxes()
ggplot2::ggsave(plot = d, 'dimplot.pbmc_good.png')

sub = 
d = DimPlot(s, group.by='Phase')+NoAxes()
ggplot2::ggsave(plot = d, 'dimplot.cc.pbmc_good.png')

## Remove variation associated with cell-cycle and with percent.mt

sub@meta.data$celltype = NA
sub@meta.data[which(sub@meta.data$seurat_clusters %in% c(7,0)),'celltype']='CD14.Mono'
sub@meta.data[which(sub@meta.data$seurat_clusters %in% c(10)),'celltype']='CD16.Mono'
sub@meta.data[which(sub@meta.data$seurat_clusters %in% c(14)),'celltype']='DCs'




dat = sub@meta.data
dat$XIST=sub[['SCT']]@data['XIST',rownames(dat)]
dat$DDX3Y=sub[['SCT']]@data['DDX3Y',rownames(dat)]
dat$PPBP=sub[['SCT']]@data['PPBP',rownames(dat)]
dat$LYZ=sub[['SCT']]@data['LYZ',rownames(dat)]
dat$MS4A1=sub[['SCT']]@data['MS4A1',rownames(dat)]
dat$CD4=sub[['SCT']]@data['CD4',rownames(dat)]
dat$MBP=sub[['SCT']]@data['MBP',rownames(dat)]
dat$CNP=sub[['SCT']]@data['CNP',rownames(dat)]
dat$CLDN5=sub[['SCT']]@data['CLDN5',rownames(dat)]

genes=c('XIST','DDX3Y','PPBP','LYZ','MS4A1','CD4','MBP','CNP','CLDN5')
to_plot = reshape2::melt(dat[,c('orig.ident',genes)])
write.table(to_plot, 'plot_genes.dat',quote=F)


a = aggregate(dat$celltype, by=list(dat$sex,dat$condition),FUN=table)
write.table(a, 'cell_proportions.dat',quote=F)

## Differential expression
Idents(sub)='celltype'
first=T
for(i in unique(celltype)){
	s = subset(sub, idents=i)
	Idents(s)='condition'
	m = FindMarkers(s, ident.1='MS',test.use='MAST')
	m$cluster = i
	m$genes=rownames(m)
	if(first){
		first=F
		df=m
	}else{
		df=rbind(df, m)
	}
}
write.table(df, 'ms_de.dat',quote=F)

avg = AverageExpression(sub, )


for(i in unique(fr@meta.data$celltype)){
	sub = subset(fr, idents=i)
	Idents(sub)='condition'
	marks = FindMarkers(sub,ident.1='MS',test.use='MAST',logfc.threshold=.1)
	marks$cluster = i
	marks$genes = rownames(marks)
	if(first){
		first=F
		out=marks
	}else{
		out=rbind(out,marks)
	}
}
for(u in unique(out$cluster)){
	sub = subset(out, cluster==u)
	write.table(sub , paste0('DEGs/',gsub(' ','\\.',u),'.degs.txt'),quote=F)
}

f1 = FindMarkers(obj,ident.1='S1',ident.2='S2',test.use='MAST',logfc.threshold=0)
f2 = FindMarkers(obj,ident.1='S1',ident.2='S3',test.use='MAST',logfc.threshold=0)
f3 = FindMarkers(obj,ident.1='S1',ident.2='S4',test.use='MAST',logfc.threshold=0)
f4 = FindMarkers(obj,ident.1='S2',ident.2='S3',test.use='MAST',logfc.threshold=0)
f5 = FindMarkers(obj,ident.1='S3',ident.2='S4',test.use='MAST',logfc.threshold=0)
f6 = FindMarkers(obj,ident.1='S2',ident.2='S4',test.use='MAST',logfc.threshold=0)

d=cbind(f1$avg_logFC,f2[rownames(f1),'avg_logFC'],f3[rownames(f1),'avg_logFC'],f4[rownames(f1),'avg_logFC'],f5[rownames(f1),'avg_logFC'],f6[rownames(f1),'avg_logFC'])
colnames(d)=c('S1vS2','S1vS3','S1vS4','S2vS3','S3vS4','S2vS4')
rownames(d) = rownames(f1)