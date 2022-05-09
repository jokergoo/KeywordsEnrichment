
library(GetoptLong)
library(ComplexHeatmap)
library(cola)
library(simplifyEnrichment)

anno_GO_keywords = function(split, genes, id_mapping = NULL, org_db = "org.Hs.eg.db",
	ora_fun = function(x) clusterProfiler::enrichGO(x, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 1), 
	padj_cutoff = 0.05, go_id_column = NULL, padj_column = NULL, verbose = TRUE, ...) {

	if(!is.null(id_mapping)) {
		genes = id_mapping[genes]
	}

	id_mapping = cola:::guess_id_mapping(genes, org_db = org_db)
	if(is.null(id_mapping)) {

	} else if(is.function(id_mapping)) {
		genes = id_mapping(genes)
	} else {
		genes = id_mapping[genes]
	}

	gene_list = split(genes, split)
	gene_list = lapply(gene_list, function(x) {
		x[!is.na(x)]
	})

	lt = lapply(names(gene_list), function(nm) {
		x = gene_list[[nm]]
		qqcat("Perform over-representation analysis on group '@{nm}', @{length(x)} genes...\n")
		as.data.frame(ora_fun(x))
	})
	names(lt) = names(gene_list)

	if(is.null(go_id_column)) {
		go_id_column = which(sapply(lt[[1]], function(x) all(grepl("^GO:\\d+$", x))))[1]
		if(length(go_id_column) == 0) {
			if(!is.null(rownames(lt[[1]]))) {
				go_id_column = rownames
				if(is.null(rownames(lt[[1]]))) {
					stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
				}
				if(verbose) {
					qqcat("Use row names of the data frame as `go_id_column`.\n")
				}
			} else {
				stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
			}
		} else {
			if(verbose) {
				qqcat("Use column '@{colnames(lt[[1]])[go_id_column]}' as `go_id_column`.\n")
			}
		}
	}

	test_padj_column = function(cn) {
		test_cn = c("p.adjust", "p_adjust", "padjust", "padj", "fdr", "FDR", "BH", "p.value", "p-value", "pvalue", "p_value")
		for(x in test_cn) {
			ind = which(cn %in% x)
			if(length(ind)) {
				return(ind[1])
			}
		}
		return(NULL)
	}

	if(is.null(padj_column)) {
		cn = colnames(lt[[1]])
		ind = test_padj_column(cn)
		if(length(ind)) {
			padj_column = ind
			if(verbose) {
				qqcat("Use column '@{colnames(lt[[1]])[padj_column]}' as `padj_column`.\n")
			}
		} else {
			stop_wrap("Cannot find the column the contains adjusted p-values in the data frames. Please explicitly set argument `padj_column`.")
		}
	}

	if(is.function(go_id_column)) {
		go_list = lapply(lt, function(tb) {
			go_id_column(tb)[tb[, padj_column] <= padj_cutoff]
		})
	} else {
		go_list = lapply(lt, function(tb) {
			tb[tb[, padj_column] <= padj_cutoff, go_id_column]
		})
	}

	anno_word_cloud_from_GO(split, go_list, ...)
}




