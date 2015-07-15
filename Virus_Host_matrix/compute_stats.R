source("lp_brim.R")
for (i in 1:99){
	cat('|i|',i,'|');
	current_matrix<-as.matrix(mat_2_whole_random$perm[[i]])
	colnames(current_matrix)<-colnames(mat_2_whole)
	rownames(current_matrix)<-rownames(mat_2_whole)
	current_matrix_modularity<-bBRIM(current_matrix)
}