setwd("~/project/development/KeywordsEnrichment/inst/extdata")

library(udpipe)

if(!exists("m_eng_ewt_loaded")) {
	m_eng_ewt = udpipe_download_model(language = "english")
	m_eng_ewt_path = m_eng_ewt$file_model
	m_eng_ewt_loaded = udpipe_load_model(file = m_eng_ewt_path)
}

library(pluralize)
library(textstem)
library(GetoptLong)
library(IRanges)
library(Matrix)

SINGLE_WORDS = read.table("single_word.csv", header = TRUE, sep = ",")
SINGLE_WORDS = SINGLE_WORDS[SINGLE_WORDS$action != "x", ]
lt = strsplit(SINGLE_WORDS$token, ",")
SINGLE_WORDS_MAP = structure(rep(ifelse(SINGLE_WORDS$action == "", SINGLE_WORDS$keyword, SINGLE_WORDS$action), times = sapply(lt, length)), names = unlist(lt))

PHRASES = read.table("phrase_list.csv", sep = "\t", quote = "")
l = PHRASES[, 1] != ""
PHRASES_MAP = structure(PHRASES[l, 1], names = PHRASES[l, 2])
PHRASES_WORDS = strsplit(PHRASES[, 2], " ")

PHRASES_N_WORDS = sapply(PHRASES_WORDS, length)

uncaptitalize = function(str) {
	substr(str, 0, 1) = tolower(substr(str, 0, 1))
	str
}

docs_to_keywords = function(docs, doc_id, udpipe_model) {

	docs[is.na(docs)] = ""

	# A-B is parsed as three tokens: A, -, B, so we removed "-" or "_"
	docs = gsub("_|-", " ", docs)

	message("- tokenization...")
	text_annotated = udpipe_annotate(udpipe_model, x = docs, doc_id = doc_id, parser = "none")
	text_annotated = as.data.frame(text_annotated)

	token = text_annotated$token
	
	# remove punctuation
	l = grepl("^[[:punct:]][[:alnum:]]{2,}", token)
	token[l] = gsub("^[[:punct:]]", "", token[l])
	l = grepl("[[:alnum:]]{2,}[[:punct:]]$", token)
	token[l] = gsub("[[:punct:]]$", "", token[l])


	message("- uncaptitalize words...")
	# a word contains vowel and no numbers
	l = grepl("^[A-Z][a-z]+$", text_annotated$token) & 
    grepl("[oeuiaOEUIA]", text_annotated$token)
	sw = uncaptitalize(token[l])
	l2 = sw %in% words::words[, 1]
	token[which(l)[l2]] = sw[l2]


	# manually validated capitalized words
	cp = readLines("capitalized_words.txt")
	l = text_annotated$token %in% cp
	token[l] = tolower(token[l])

	message("- singularization...")
	# change plural to singulars
	l = grepl("s$", token) & 
	    !grepl("sis$", token) & 
	    !grepl("'s$", token) & 
	    !grepl("tis$", text_annotated$token) & 
	    text_annotated$upos %in% c("NOUN", "PROPN", "INTJ")
	exclude_singular = scan("singular_exclude.txt", what = "character")
	l = l & !token %in% exclude_singular


	w = unique(token[l])
	singular_map = structure(singularize(w), names = w)

	token[l] = singular_map[token[l]]
	text_annotated$token_normalized = token
	text_annotated$lemma2 = lemmatize_words(token)


	lemma_map = read.table("adjust_words.txt")
	lemma_map = structure(lemma_map[, 2], names = lemma_map[, 1])
	l = text_annotated$lemma2 %in% names(lemma_map)
	text_annotated$lemma2[l] = lemma_map[ text_annotated$lemma2[l] ]


	## words may still be wrongly singularized by `lemmatize_words()`
	l = grepl("sis$", token) | 
	    grepl("'s$", token) | 
	    grepl("tis$", text_annotated$token)
	l = l | token %in% exclude_singular
	l = l & text_annotated$upos %in% c("NOUN", "PROPN", "INTJ")
	text_annotated$lemma2[l] = uncaptitalize(text_annotated$token[l])


	l = text_annotated$lemma2 == ""
	text_annotated$lemma2[l] = text_annotated$token[l]

	text_annotated$token_id = as.numeric(text_annotated$token_id)

	lt = split(text_annotated, text_annotated$doc_id)
	lt2 = lapply(lt, function(term) {

		qqcat(strrep("\b", 100))
		qqcat("doc: @{term$doc_id[1]}...")

		phrases = tapply(seq_len(nrow(term)), paste(term$paragraph_id, term$sentence_id, sep = "_"), function(ind) {
			token_id = term$token_id[ind]
			x = Rle(cumsum(c(1, diff(token_id))) - cumsum(rep(1, length(token_id))))
			ir = IRanges(start = ind[start(x)], end = ind[end(x)])
			kw2 = character(0)
			if(length(ir) > 0) {
				s = start(ir)
				e = end(ir)
				for(k in seq_along(s)) {
					kw2 = c(kw2, get_phrase(term$lemma2[ seq(s[k], e[k]) ]))
				}
			}

			unique(kw2)
		})
		phrases = unname(unique(unlist(phrases)))
		ll = phrases %in% PHRASES_MAP
		if(any(ll)) {
			phrases[ll] = PHRASES_MAP[ phrases[ll] ]
		}
		kw = c(phrases, SINGLE_WORDS_MAP[term$token_normalized])
		unique(kw[!is.na(kw)])
	})

	lt2
}

lt_to_tdm = function(lt2) {
	message("- save to the dgTMatrix object...")
	all_docs = unique(names(lt2))
	all_terms = unique(unlist(lt2))
	all_docs_map = structure(seq_along(all_docs), names = all_docs)
	all_terms_map = structure(seq_along(all_terms), names = all_terms)

	i = all_terms_map[unlist(lt2)]
	j = all_docs_map[rep(names(lt2), times = sapply(lt2, length))]

	tdm = sparseMatrix(i, j, x = 1, dims = c(length(all_terms),length(all_docs)))
	rownames(tdm) = all_terms
	colnames(tdm) = all_docs
	
	as(tdm, "dgTMatrix")
}

.get_phrase = function(words) {
	l = rep(TRUE, length(PHRASES_WORDS))
	kw = integer(0)
	for(i in seq_along(words)) {
		ind = PHRASES_N_WORDS >= i
		if(length(ind) == 0 ) {
			break
		}
		l[ind] = l[ind] & sapply(PHRASES_WORDS[ind], function(x) x[i] == words[i])
		if(i > 1) {
			kw = c(kw, which(PHRASES_N_WORDS == i & l))
		}
		l[PHRASES_N_WORDS == i] = FALSE
		if(!any(l)) {
			break
		}
	}
	if(length(kw)) {
		sapply(PHRASES_WORDS[kw], paste, collapse = " ")
	} else {
		character(0)
	}
}

get_phrase = function(words) {

	if(length(words) <= 1) {
		return(character(0))
	}

	nw = length(words)
	kw = character(0)
	for(i in seq_along(words)[-nw]) {
		kw = c(kw, .get_phrase(words[i:nw]))
	}
	kw
}


library(GO.db)
go_text = Term(GOTERM)
lt = docs_to_keywords(go_text, names(go_text), udpipe_model)
saveRDS(lt, file = "keywords_from_go_name.rds")

go_def = Definition(GOTERM)
lt = docs_to_keywords(go_def, names(go_def), udpipe_model)
saveRDS(lt, file = "keywords_from_go_definition.rds")


dag = ontology_kw()
docs = mcols(dag)[, "description"]
names(docs) = dag_all_terms(dag)
lt = docs_to_keywords(docs, names(docs), udpipe_model)
saveRDS(lt, file = "keywords_from_uniport_kw.rds")


library(GeneSummary)
tb = loadGeneSummary(organism = 9606)
docs = tapply(tb$Gene_summary, tb$Gene_ID, function(x) x[1])
gene_ids = names(docs)
docs = as.vector(docs)
names(docs) = gene_ids
lt = docs_to_keywords(docs, names(docs), udpipe_model)
saveRDS(lt, file = "keywords_from_refseq_genes.rds")


dag = ontology_pw()
docs = mcols(dag)[, "definition"]
names(docs) = dag_all_terms(dag)
lt = docs_to_keywords(docs, dag_all_terms(dag), udpipe_model)
saveRDS(lt, file = "keywords_from_pathway_ontology.rds")


dag = ontology_rdo()
docs = mcols(dag)[, "definition"]
names(docs) = dag_all_terms(dag)
lt = docs_to_keywords(docs, dag_all_terms(dag), udpipe_model)
saveRDS(lt, file = "keywords_from_disease_ontology.rds")


dag = ontology_vt()
docs = mcols(dag)[, "definition"]
names(docs) = dag_all_terms(dag)
lt = docs_to_keywords(docs, dag_all_terms(dag), udpipe_model)
saveRDS(lt, file = "keywords_from_vertebrate_trait_ontology.rds")


dag = ontology_hp()
docs = mcols(dag)[, "definition"]
names(docs) = dag_all_terms(dag)
lt = docs_to_keywords(docs, dag_all_terms(dag), udpipe_model)
saveRDS(lt, file = "keywords_from_human_phenotype_ontology.rds")




