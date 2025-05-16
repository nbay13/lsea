
#' @export clr.transform
clr.transform <- function(mat){
	# optimal according to https://doi.org/10.1016/j.chemolab.2021.104248
	min_const <- min(mat[mat != 0]) * 0.65
	temp_mat <- mat
	temp_mat[temp_mat == 0] <- min_const
	mat_clr <- compositions::clr(temp_mat)
}

#' @export two.group.row.test
two.group.row.test <- function(data, labels, test = c("t", "w"), var_equal = FALSE, paired = FALSE, adjust_method = "fdr"){
	test <- match.arg(test)
	if(!is.factor(labels)) {
	 	message("labels are not factors, converting to factors")
	 	labels <- factor(labels)
	} 
	if(test == "t") ev <- paste("t.test(data[i,] ~ factor(labels), var.equal = ", var_equal, ", paired = ", paired,  ")", sep = "")
	else if(test == "w") ev <- paste("wilcox.test(data[i,] ~ factor(labels), paired = ", paired, ")", sep = "")
	df <- data.frame(matrix(nrow = nrow(data), ncol = 5))
	colnames(df) <- c("stat", "mean1", "mean2", "dm", "pvalue")
	rownames(df) <- rownames(data)
	inds <- which(apply(data, 1, var) > 0)
	for(i in inds){
		t_res <- eval(parse(text = ev))
		mean1 <- mean(data[i, labels == levels(labels)[1]])
		mean2 <- mean(data[i, labels == levels(labels)[2]])
		df[i,] <- c(t_res$statistic, mean1, mean2, mean1 - mean2, t_res$p.value)
	}
	df$padj <- p.adjust(df$pvalue, method = adjust_method)
	return(df)
}

#' @export lsea.old
lsea.old <- function(de_tbl, rnk_name, var_name = "Species", rownames = TRUE, nperm = 10000, minSize = 2, reformat = TRUE){
	if(rownames) anno_df <- annotate.lipid.species.old(rownames(de_tbl))
	else if(!rownames & !(var_name %in% colnames(de_tbl))) message("Error: Provide valid column name of lipid species names")
	else anno_df <- annotate.lipid.species.old(de_tbl[[var_name]])

	anno_df$Chain.new <- anno_df$Chain
	anno_df$Chain.new[anno_df$Longest.Tail %in% c(12, 14, 15, 16)] <- "12-16"
	anno_df$Chain.new[anno_df$Longest.Tail %in% c(17, 18, 20)] <- "17-20"
	anno_df$Chain.new[anno_df$Longest.Tail %in% c(22, 24, 26)] <- "22-26"

	tail_lengths <- split(anno_df$Species, anno_df$Longest.Tail)
	saturations <- split(anno_df$Species, anno_df$Saturation)
	classes <- split(anno_df$Species, anno_df$Class)
	chains <- split(anno_df$Species, anno_df$Chain.new)

	high_res <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, anno_df$Chain.new, sep = "_"))
	high_res1 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, sep = "_"))
	high_res2 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Chain.new, sep = "_"))
	high_res3 <- split(anno_df$Species, paste(anno_df$Saturation, anno_df$Chain.new, sep = "_"))
	high_res4 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, anno_df$Longest.Tail, sep = "_"))
	lipid_lists <- c(tail_lengths, saturations, classes, chains, high_res, high_res1, high_res2, high_res3, high_res4)
	rnk <- setNames(de_tbl[[rnk_name]], anno_df$Species)
	rnk[is.na(rnk)] <- 0
	res <- fgsea::fgseaSimple(stats = rnk, pathways = lipid_lists, nperm = nperm, minSize = minSize)
	res <- res %>% dplyr::arrange(desc(NES))
	new_edge <- apply(res, 1, function(x){
		paste(unlist(x$leadingEdge), collapse = ",")
	})
	if(reformat){
		res$leadingEdge <- new_edge
	}
	return(data.frame(res))
}

#' @export lsea
lsea <- function(de_tbl, rnk_name, var_name = "Species", rownames = TRUE, nperm = 10000, minSize = 2, reformat = TRUE){
	if(rownames){
		anno_df <- get.acyl.tails(rownames(de_tbl))
		lipid_names <- gsub(" |-|/|\\\\|:|;|_|~", "\\.", rownames(de_tbl))
	} 
	else if(!rownames & !(var_name %in% colnames(de_tbl))) message("Error: Provide valid column name of lipid species names")
	else {
		anno_df <- get.acyl.tails(de_tbl[[var_name]])
		lipid_names <- gsub(" |-|/|\\\\|:|;|_|~", "\\.", de_tbl[[var_name]])
	}
	pa_df <- anno_df[anno_df$Class == "PA",]
	anno_df <- anno_df[!anno_df$Class %in% c("PA", "Cholesterol"),]

	classes <- split(anno_df$Species, anno_df$Class)
	classes[["PA"]] <- pa_df$Species
	tail_lengths <- split(anno_df$Species, anno_df$Total.Carbons)
	saturations <- split(anno_df$Species, anno_df$Saturation)
	chains <- split(anno_df$Species, anno_df$Chain.Group)

	high_res <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, anno_df$Chain.Group, sep = "_"))
	high_res1 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, sep = "_"))
	high_res2 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Chain.Group, sep = "_"))
	high_res3 <- split(anno_df$Species, paste(anno_df$Saturation, anno_df$Chain.Group, sep = "_"))
	high_res4 <- split(anno_df$Species, paste(anno_df$Class, anno_df$Saturation, anno_df$Total.Carbons, sep = "_"))

	lipid_lists <- c(tail_lengths, saturations, classes, chains, high_res, high_res1, high_res2, high_res3, high_res4)
	rnk <- setNames(de_tbl[[rnk_name]], lipid_names)
	rnk[is.na(rnk)] <- 0
	res <- fgsea::fgseaSimple(stats = rnk, pathways = lipid_lists, nperm = nperm, minSize = minSize)
	res <- res %>% dplyr::arrange(desc(NES))
	new_edge <- apply(res, 1, function(x){
		paste(unlist(x$leadingEdge), collapse = ",")
	})
	if(reformat){
		res$leadingEdge <- new_edge
	}
	return(data.frame(res))
}

# from https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @export gg.colors
gg.colors <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 ggplot aes position_dodge
#' @export structure.enrichment.plot
structure.enrichment.plot <- function(de_tbl, anno_tbl, group_names, class = NULL, p_thresh = 0.05, color_pal = NULL, size_range = c(1, 4), facet_rows = 3){
	if(is.null(color_pal)) color_pal <- gg.colors(2)
	merge_df <- data.frame(cbind(de_tbl, anno_tbl[rownames(de_tbl),]))
	temp <- merge_df %>% dplyr::filter(padj < p_thresh) %>% dplyr::group_by(Class, Longest.Tail, Total.DBs, dm > 0) %>% dplyr::count() %>% data.frame()
	colnames(temp)[4] <- "sign"
	temp$label <- temp$n
	temp$label[temp$n < 2] <- ""
	temp <- temp %>% tidyr::complete(Longest.Tail, Total.DBs, sign) %>% data.frame()
	temp$Longest.Tail <- factor(temp$Longest.Tail)
	if(is.null(class)){
		gg <- ggplot(temp %>% na.omit(), aes(x = Longest.Tail, y = Total.DBs, size = n, color = sign, group = sign)) + 
		ggplot2::geom_point(position = position_dodge(width = 0.8)) +
		ggplot2::geom_text(aes(label= label), position = position_dodge(width = 0.8), color = "black") +
		ggplot2::theme_classic() +
		ggplot2::scale_color_manual(values = color_pal, limits = c(TRUE, FALSE), labels = group_names) +
		ggplot2::facet_wrap(~ Class, axes = "all_x", nrow = facet_rows) +
		ggplot2::theme(legend.position = "top") +
		ggplot2::scale_size(range = size_range) + 
		ggplot2::labs(size = "# of lipid species", color = "Increase in")		
	} else {
		gg <- ggplot(temp %>% dplyr::filter(Class == class) %>% na.omit(), aes(x = Longest.Tail, y = Total.DBs, size = n, color = sign, group = sign)) + 
		ggplot2::geom_point(position = position_dodge(width = 0.8)) +
		ggplot2::geom_text(aes(label= label), position = position_dodge(width = 0.8), color = "black") +
		ggplot2::theme_classic() +
		ggplot2::scale_color_manual(values = color_pal, limits = c(TRUE, FALSE), labels = group_names) +
		ggplot2::theme(legend.position = "top") +
		ggplot2::scale_size(range = size_range) + 
		ggplot2::facet_wrap(~ Class) +
		ggplot2::labs(size = "# of lipid species", color = "Increase in")
	}
	return(gg)
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

#' @export get.chain.group
get.chain.group.2 <- function(lengths){
	group_list <- list(
		"12-16" = as.character(c(12:16)), 
		"17-21" = as.character(c(17:21)), 
		"22-26" = as.character(c(22:26))
		)
	inv_list <- inverse.list(group_list)
	return(unlist(sapply(lengths, function(x){
		if(x %in% names(inv_list)) inv_list[[as.character(x)]]
		else NA
	})))
}

# based on inverseList from BioCor
#' @export inverse.list
inverse.list <- function(x){
	stopifnot(length(names(x)) == length(x))
    stopifnot(all(sapply(x, function(x) {
        is.character(x) || is.na(x)
    })))
	values <- unlist(x, use.names = FALSE)
    names <- rep(names(x), lengths(x))
    split(names, values)
}

#' @export annotate.lipid.species.old
annotate.lipid.species.old <- function(input_names){
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
	rownames(structure_anno) <- input_names
	structure_anno$Species <- lipid_names
	structure_anno[,2:4] <- apply(structure_anno[,2:4], 2, as.numeric)
	structure_anno$Category <- get.lipid.category(structure_anno$Class)
	structure_anno$Chain <- get.chain.group(structure_anno$Longest.Tail)
	structure_anno[structure_anno$Class == "PA","Chain"] <- NA
	return(structure_anno[,c("Species", "Class", "Category", "Total.Carbons", "Longest.Tail", "Total.DBs", "Saturation", "Chain")])
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
	tail_length_list <- vector("list", length(10:28))
	names(tail_length_list) <- as.character(10:28)
	dbs_list <- vector("list", length(0:12))
	names(dbs_list) <- as.character(0:12)
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
			tail_length_list[[gsub("[A-z]*", "", temp[[i]][4])]] <- c(tail_length_list[[gsub("[A-z]*", "", temp[[i]][4])]], lipid_names[i])
			dbs_list[[temp[[i]][5]]] <- c(dbs_list[[temp[[i]][5]]], lipid_names[i])
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
			tail_length_list[[temp[[i]][3]]] <- c(tail_length_list[[temp[[i]][3]]], lipid_names[i])
			tail_length_list[[temp[[i]][5]]] <- c(tail_length_list[[temp[[i]][5]]], lipid_names[i])
			dbs_list[[temp[[i]][4]]] <- c(dbs_list[[temp[[i]][4]]], lipid_names[i])
			dbs_list[[temp[[i]][6]]] <- c(dbs_list[[temp[[i]][6]]], lipid_names[i])
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
			tail_length_list[[gsub("d", "", temp[[i]][2])]] <- c(tail_length_list[[gsub("d", "", temp[[i]][2])]], lipid_names[i])
			tail_length_list[[temp[[i]][4]]] <- c(tail_length_list[[temp[[i]][4]]], lipid_names[i])
			dbs_list[[temp[[i]][3]]] <- c(dbs_list[[temp[[i]][3]]], lipid_names[i])
			dbs_list[[temp[[i]][5]]] <- c(dbs_list[[temp[[i]][5]]], lipid_names[i])
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
			tail_length_list[[temp[[i]][2]]] <- c(tail_length_list[[temp[[i]][2]]], lipid_names[i])
			dbs_list[[temp[[i]][3]]] <- c(dbs_list[[temp[[i]][3]]], lipid_names[i])
			if(as.numeric(temp[[i]][3]) > 1){
				structure_anno[i,5] <- "PUFA"
			} else if(as.numeric(temp[[i]][3]) == 1){
				structure_anno[i,5] <- "MUFA"
			} else {
				structure_anno[i,5] <- "SFA"
			}
		}
	}
	rownames(structure_anno) <- input_names
	structure_anno$Species <- lipid_names
	structure_anno[,2:4] <- apply(structure_anno[,2:4], 2, as.numeric)
	structure_anno$Category <- get.lipid.category(structure_anno$Class)
	structure_anno$Chain <- get.chain.group(structure_anno$Longest.Tail)
	structure_anno[structure_anno$Class == "PA","Chain"] <- NA


	tail_length_list <- tail_length_list[unlist(lapply(tail_length_list, function(x) length(x) > 0))]
	dbs_list <- dbs_list[unlist(lapply(dbs_list, function(x) length(x) > 0))]

	names(dbs_list) <- paste("contains.DB", names(dbs_list), sep = ".")
	names(tail_length_list) <- paste("contains.length", names(tail_length_list), sep = ".")


	for(i in 1:length(dbs_list)){
		temp <- names(dbs_list)[i]
		structure_anno[[temp]] <- FALSE
		structure_anno[structure_anno$Species %in% dbs_list[[temp]],temp] <- TRUE
	}
	for(i in 1:length(tail_length_list)){
		temp <- names(tail_length_list)[i]
		structure_anno[[temp]] <- FALSE
		structure_anno[structure_anno$Species %in% tail_length_list[[temp]],temp] <- TRUE
	}

	structure_anno[["contains.SFA"]] <- structure_anno[,grep("contains\\.DB\\.[0]", colnames(structure_anno))]
	structure_anno[["contains.MUFA"]] <- structure_anno[,grep("contains\\.DB\\.[1]", colnames(structure_anno))]
	structure_anno[["contains.PUFA"]] <- apply(structure_anno[,grep("contains\\.DB\\.[2-9]", colnames(structure_anno))], 1, any)

	structure_anno[["contains.12-16"]] <- apply(structure_anno[,grep("contains\\.length\\.1[2-6]", colnames(structure_anno))], 1, any)
	structure_anno[["contains.17-20"]] <- apply(structure_anno[,grep("contains\\.length\\.(1[7-9]|20)", colnames(structure_anno))], 1, any)
	structure_anno[["contains.21-26"]] <- apply(structure_anno[,grep("contains\\.length\\.2[1-6]", colnames(structure_anno))], 1, any)

	structure_anno <- structure_anno[,!colnames(structure_anno) %in% c("Chain", "Saturation")]
	structure_anno <- structure_anno[,c("Species", setdiff(colnames(structure_anno), "Species"))]
	return(structure_anno)
}

#' @export get.tail.saturation
get.tail.saturation <- function(n_db){
	res <- n_db
	res[n_db == 0] <- "SFA"
	res[n_db == 1] <- "MUFA"
	res[n_db > 1] <- "PUFA"
	return(res)
}


#' @export get.acyl.tails
get.acyl.tails <- function(input_names){
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
	structure_anno <- data.frame(matrix(nrow = length(class_name), ncol = 3))
	extras <- list()
	colnames(structure_anno) <- c("Class","Total.Carbons", "Total.DBs")
	structure_anno[,1] <- class_name
	# This script works based on specific class names. If new input data has different or new 
	# class names than the script or names must be edited.
	for(i in 1:length(class_name)){
		if(class_name[i] %in% c("TG", "TAG")){
			# For TAG class with three tails
			structure_anno[i,2] <- gsub("[A-z]*", "", temp[[i]][4])
			fa_dbs <- as.numeric(temp[[i]][5])
			structure_anno[i,3] <- fa_dbs
		} else if(class_name[i] == "PE" & temp[[i]][2] %in% c("P","O")){
			# For PE classes with two chains and extra character for PE-P and PE-O
			structure_anno[i,1] <- paste("PE", temp[[i]][2], sep = ".")
			structure_anno[i,2] <- as.numeric(temp[[i]][3]) 
			extras[[as.character(i)]] <- c(paste("PE", temp[[i]][2], sep = "."), as.numeric(temp[[i]][5]), as.numeric(temp[[i]][6]))
			structure_anno[i,3] <- as.numeric(temp[[i]][4])
		} else if(class_name[i] %in% two_chain){
			# For other classes with two chains
			structure_anno[i,2] <- as.numeric(gsub("d", "", temp[[i]][2]))
			if(structure_anno[i,1] == "Cer" & as.numeric(temp[[i]][3]) == 0) structure_anno[i,1] <- "dhCer"
			extras[[as.character(i)]] <- c(structure_anno[i,1], as.numeric(temp[[i]][4]), as.numeric(temp[[i]][5]))
			structure_anno[i,3] <- as.numeric(temp[[i]][3])
		} else if(class_name[i] == "Chol"){
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- NA
		} else {
			# For classes with one chain
			structure_anno[i,2] <- temp[[i]][2]
			structure_anno[i,3] <- temp[[i]][3]
		}
	}
	#rownames(structure_anno) <- lipid_names
	extra_df <- data.frame(t(as.data.frame(extras)))
	colnames(extra_df) <- colnames(structure_anno)
	extra_df$Species <- lipid_names[as.numeric(names(extras))]
	rownames(extra_df) <- input_names[as.numeric(names(extras))]
	structure_anno$Species <- lipid_names
	rownames(structure_anno) <- input_names
	structure_anno <- data.frame(rbind(structure_anno, extra_df))
	structure_anno[,2:3] <- apply(structure_anno[,2:3], 2, as.numeric)
	structure_anno$Category <- get.lipid.category(structure_anno$Class)
	structure_anno$Chain <- get.chain.group(structure_anno$Total.Carbons)
	structure_anno$Chain.Group <- get.chain.group.2(structure_anno$Total.Carbons)
	structure_anno$Saturation <- get.tail.saturation(structure_anno$Total.DBs)
	return(structure_anno[order(structure_anno$Species),c("Species", "Class", "Category", "Total.Carbons", "Chain", "Chain.Group", "Total.DBs", "Saturation")])
}