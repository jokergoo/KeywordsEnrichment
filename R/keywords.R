



n_grams = function(x, k) {
    n = length(x)
    if(n < k) {
    	return(character(0))
    } else {
    	si = seq(1, n - k + 1)
    	ei = seq(k, n)
    	sapply(seq_along(si), function(i) {
    		paste(x[seq(si[i], ei[i])], collapse = " ")
    	})
    }
}


tokenizer = function(n_gram = 1) {
	if(n_gram == 1) {
		function(x) {
			w = x$content
			w = strsplit(w, "[[:punct:]]|[[:space:]]+", perl = TRUE)[[1]]
			w = w[!grepl("^\\d*$", w)]
			unique(w)
		}
	} else if(n_gram == 2) {
		function(x) {
			w = x$content
			w = strsplit(w, "[[:punct:]]|[[:space:]]+", perl = TRUE)[[1]]
			l = !grepl("^\\d*$", w)
			rl = rle(l)
			i1 = cumsum(rl$length)
			i2 = c(0, i1) + 1; i2 = i2[1:length(i1)]
			i1 = i1[rl$value]
			i2 = i2[rl$value]
			unique(unlist(lapply(1:length(i1), function(i) {
				n_grams(w[seq(i2[i], i1[i])], 2)
			})))
		}
	} else if(n_gram == 3) {
		function(x) {
			w = x$content
			w = strsplit(w, "[[:punct:]]|[[:space:]]+", perl = TRUE)[[1]]
			l = !grepl("^\\d*$", w)
			rl = rle(l)
			i1 = cumsum(rl$length)
			i2 = c(0, i1) + 1; i2 = i2[1:length(i1)]
			i1 = i1[rl$value]
			i2 = i2[rl$value]
			unique(unlist(lapply(1:length(i1), function(i) {
				n_grams(w[seq(i2[i], i1[i])], 3)
			})))
		}
	}
}

make_term_document_matrix = function(term, n_gram = 1) {


	docs = VCorpus(VectorSource(term))

	sp = stopwords()
	sp2 = sapply(sp, function(x) substr(x, 1, 1) = toupper(substr(x, 1, 1)))
	docs = tm_map(docs, removeWords, c(sp, sp2))

	tdm = TermDocumentMatrix(
		docs,
		control = list(
			tokenize = tokenizer(n_gram),
			wordLengths = c(1, Inf),
			tolower = FALSE
		)
	)
	return(tdm)
}

prepare_keywords_tdm_from_GO = function(use_desc = FALSE, ontology = c("BP", "CC", "MF"), ...) {

	all_go = as.list(GO.db::GOTERM)

	onto = sapply(all_go, slot, "Ontology")
	term = sapply(all_go, slot, "Term")
	l = onto %in% ontology
	term = term[l]
	term_id = names(term)

	if(use_desc) {
		suppressMessages(term <- select(GO.db::GO.db, keys = term_id, columns = "DEFINITION")$DEFINITION)
		tdm = make_term_document_matrix(term, ...)
	} else {
		tdm = make_term_document_matrix(term, ...)
	}

	colnames(tdm) = term_id
	attr(tdm, "GO_dbInfo") = GO.db::GO_dbInfo()
	
	tdm
}


keyword_enrichment = function(term_id, tdm, min_bg = 5, min_term = 2) {
	tdm2 = tdm[slam::row_sums(tdm) >= min_bg, , drop = FALSE]

	l = colnames(tdm2) %in% term_id

	l2 = row_sums(tdm2[, l]) >= min_term
	tdm2 = tdm2[l2, , drop = FALSE]

	n = nrow(tdm2)
	n_term = numeric(n)
	n_bg = numeric(n)
	p = numeric(n)
	for(i in seq_len(n)) {
			if(i %% 100 == 0 || i == n) {
				message(strrep("\r", 100), appendLF = FALSE)
				message(qq("keyword enrichment, @{i}/@{n}..."), appendLF = FALSE)
			}
		
		v = as.vector(tdm2[i, ])
		s11 = sum(v & l)
		if(s11 < min_term) {
			next
		}
		s12 = sum(!v & l)
		s21 = sum(v & !l)
		s22 = sum(!v & !l)

		n_term[i] = s11
		n_bg[i] = s11 + s21

		p[i] = 1 - phyper(s11-1, s11+s21, s12+s22, s11+s21)

	}

		message("")

	df = data.frame(keyword = rownames(tdm2), n_term = n_term, n_bg = n_bg, p = p)
	df = df[df$n_term >= min_term, , drop = FALSE]
	df$padj = p.adjust(df$p)
	df[order(df$padj, df$p), , drop = FALSE]
}


# go_id = random_GO(100)
# keyword_enrichment_from_GO(go_id)
keyword_enrichment_from_GO = function(go_id, min_bg = 5, min_term = 2) {

	if(is.null(env$tdm_GO)) {
		# env$tdm_GO = readRDS("~/project/development/simplifyEnrichment/inst/extdata/tdm_GO.rds")
		env$tdm_GO = readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
	}

	df = keyword_enrichment(go_id, env$tdm_GO, min_bg, min_term)
	rownames(df) = NULL
	df
}
