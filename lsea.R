
#' @export lsea
lsea <- function(de_tbl, rnk_name, var_name = NULL, rownames = TRUE, nperm = 10000, minSize = 2){
	if(rownames) anno_df <- annotate.lipid.species(rownames(de_tbl))
	else if(!rownames & is.null(var_name)) message("Error: Provide column name of lipid species if not using rownames")
	else anno_df <- annotate.lipid.species(de_tbl[[var_name]])

	anno_df$Chain.new <- anno_df$Chain
	anno_df$Chain.new[anno_df$Longest.Tail %in% c(12, 14, 15, 16)] <- "12..16"
	anno_df$Chain.new[anno_df$Longest.Tail %in% c(17, 18, 20)] <- "17..20"
	anno_df$Chain.new[anno_df$Longest.Tail %in% c(22, 24, 26)] <- "22..26"

	tail_lengths <- split(anno_df$Species, anno_df$Longest.Tail)
	saturations <- split(anno_df$Species, anno_df$Saturation)
	classes <- split(anno_df$Species, anno_df$Class)
	chains <- split(anno_df$Species, anno_df$Chain.new)

	high_res <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, anno_df$Chain.new, sep = "-"))
	high_res1 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, sep = "-"))
	high_res2 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Chain.new, sep = "-"))
	high_res3 <- split(anno_df$Species, paste(anno_df$Saturation, anno_df$Chain.new, sep = "-"))
	high_res4 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, anno_df$Longest.Tail, sep = "-"))
	lipid_lists <- c(tail_lengths, saturations, classes, chains, high_res, high_res1, high_res2, high_res3, high_res4)
	rnk <- setNames(de_tbl[[rnk_name]], rownames(plot_df))
	rnk[is.na(rnk)] <- 0
	res <- fgsea::fgsea(stats = rnk, pathways = lipid_lists, nperm = nperm, minSize = minSize) %>% arrange(NES) %>% data.frame()
	return(res)
}

#' @export get.lipid.category
get.lipid.category <- function(species){
	category_list <- list(
		"Sterol" = c("CE"), 
		"Sphingolipid" = c("Cer", "LacCER", "HexCER", "LCER", "SM", "dhCer"), 
		"Glycerolipid" = c("DG", "TG", "DAG", "TAG"), 
		"Fatty.Acyl" = c("FA"),
		"Glycerophospholipid" = c("LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "PA"), 
		"Ether" = c("PE.O", "PE.P"),
		"Cholesterol" = c("Chol")
	)
	inv_list <- inverse.list(category_list)
	return(unlist(sapply(species, function(x){
		if(x %in% names(inv_list)) inv_list[[x]]
		else NA
	})))
}

#' @export get.chain.group
get.chain.group <- function(lengths){
	group_list <- list(
		"SCFA" = as.character(c(1:4)), 
		"MCFA" = as.character(c(5:12)), 
		"LCFA" = as.character(c(13:21)), 
		"VLCFA" = as.character(c(22:32))
	)
	inv_list <- inverse.list(group_list)
	return(unlist(sapply(lengths, function(x){
		if(x %in% names(inv_list)) inv_list[[as.character(x)]]
		else NA
	})))
}

#' @export annotate.lipid.species
annotate.lipid.species <- function(input_names){
	# set up annotation vectors and lists
	two_chain <- c("PI", "PC", "PE", "PG", "PS", "PG", "DG", "DAG", "LacCER", "SM", "HexCER", "Cer", "dhCer", "PG", "PI", "PS")
	category_list <- list(
		"Sterol" = c("CE"), 
		"Sphingolipid" = c("Cer", "LacCER", "HexCER", "LCER", "SM", "dhCer"), 
		"Glycerolipid" = c("DG", "TG", "DAG", "TAG"), 
		"Fatty.Acyl" = c("FA"),
		"Glycerophospholipid" = c("LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "PA"), 
		"Ether" = c("PE.O", "PE.P"),
		"Cholesterol" = c("Chol")
	)
	# clean up lipid name format
	lipid_names <- gsub(" |-|/|\\\\|:|;|_|~", "\\.", input_names)
	# split the details in the species names by periods
	temp <- strsplit(lipid_names, "\\.")
	# remove numbers from the class names
	class_name <- unlist(lapply(temp, function(x){
		return(gsub("[0-9]*","", x[[1]]))
	}))
	if(!all(class_name %in% unlist(category_list))){
		missing <- class_name[!class_name %in% unlist(category_list)]
		w.missing <- which(!class_name %in% unlist(category_list))
		prnt <- paste(w.missing, missing, sep = ": ", collapse = "\n")
		stop(paste("Unknown lipid classes...\n", prnt))
	} 
	
	# pre-define data.frame and fill with entries using for loop
	structure_anno <- data.frame(matrix(nrow = length(class_name), ncol = 5))
	colnames(structure_anno) <- c("Class","Total.Carbons", "Longest.Tail", "Total.DBs", "Saturation")
	structure_anno[,1] <- class_name
	# This script works based on specific class names. If new input data has different or new 
	# class names than the script or names must be edited.
	for(i in 1:length(class_name)){
		if(class_name[i] %in% c("TG", "TAG")){
			# For TAG class with three tails
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- gsub("[A-z]*", "", temp[[i]][4])
			structure_anno[i,4] <- temp[[i]][3]
			fa_dbs <- as.numeric(temp[[i]][5])
			if(fa_dbs > 1 | (as.numeric(structure_anno[i,4]) - fa_dbs) > 2){
				structure_anno[i,5] <- "PUFA"
			} else if(fa_dbs == 1 & (as.numeric(structure_anno[i,4]) - fa_dbs) < 2){
				structure_anno[i,5] <- "MUFA"
			} else if(fa_dbs == 0 & as.numeric(structure_anno[i,4]) == 1){
				structure_anno[i,5] <- "MUFA"
			}
			 else if(fa_dbs == 0 & as.numeric(structure_anno[i,4]) == 0){
				structure_anno[i,5] <- "SFA"
			} else {
				structure_anno[i,5] <- "UFA"
			}
		} else if(class_name[i] == "PE" & temp[[i]][2] %in% c("P","O")){
			# For PE classes with two chains and extra character for PE-P and PE-O
			structure_anno[i,1] <- paste("PE", temp[[i]][2], sep = ".")
			structure_anno[i,2] <- sum(as.numeric(temp[[i]][3]), as.numeric(temp[[i]][5]))
			structure_anno[i,3] <- max(as.numeric(temp[[i]][3]), as.numeric(temp[[i]][5]))
			structure_anno[i,4] <- as.numeric(temp[[i]][4]) + as.numeric(temp[[i]][6])
			if(as.numeric(temp[[i]][4]) > 1 | as.numeric(temp[[i]][6]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][4]) == 1 | as.numeric(temp[[i]][6]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		} else if(class_name[i] %in% two_chain){
			# For other classes with two chains
			structure_anno[i,2] <- sum(as.numeric(gsub("d", "", temp[[i]][2])), as.numeric(temp[[i]][4]))
			structure_anno[i,3] <- max(as.numeric(gsub("d", "", temp[[i]][2])), as.numeric(temp[[i]][4]))
			structure_anno[i,4] <- as.numeric(temp[[i]][3]) + as.numeric(temp[[i]][5])
			if(structure_anno[i,1] == "Cer" & as.numeric(temp[[i]][3]) == 0) structure_anno[i,1] <- "dhCer"
			if(as.numeric(temp[[i]][3]) > 1 | as.numeric(temp[[i]][5]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1 | as.numeric(temp[[i]][5]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		} else if(class_name[i] == "Chol"){
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][2]
			structure_anno[i,4] <- NA
			structure_anno[i,5] <- NA
		} else {
			# For classes with one chain
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][2]
			structure_anno[i,4] <- temp[[i]][3]
			if(as.numeric(temp[[i]][3]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		}
	}
	#rownames(structure_anno) <- lipid_names
	structure_anno$Species <- input_names
	structure_anno[,2:4] <- apply(structure_anno[,2:4], 2, as.numeric)
	structure_anno$Category <- get.lipid.category(structure_anno$Class)
	structure_anno$Chain <- get.chain.group(structure_anno$Longest.Tail)
	structure_anno[structure_anno$Class == "PA","Chain"] <- NA
	return(structure_anno[,c("Species", "Class", "Category", "Total.Carbons", "Longest.Tail", "Total.DBs", "Saturation", "Chain")])
}
