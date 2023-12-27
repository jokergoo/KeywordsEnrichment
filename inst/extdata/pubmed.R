library(xml2)
library(parallel)
setwd("/Volumes/One Touch/pubmed")

xml_files = list.files(pattern = "gz$")

abstracts = list()
for(i in seq_along(xml_files)) {
		
	f = xml_files[i]
	xml = read_xml(f)
	docs <- xml_find_all(xml, "PubmedArticle/MedlineCitation/Article/Abstract/AbstractText") |> xml_text()
	
	docs = gsub("_|-", " ", docs)

	text_annotated = udpipe_annotate(m_eng_ewt_loaded, x = docs, doc_id = seq_along(docs), parser = "none")
	text_annotated = as.data.frame(text_annotated)

	# remove NA document
	text_annotated = text_annotated[text_annotated$sentence != "NA", ]

	# we use lemma normalization from another package
	library(textstem)
	text_annotated$lemma2 = lemmatize_words(text_annotated$token)

	l = text_annotated$lemma2 %in% names(lemma_map)
	text_annotated$lemma2[l] = lemma_map[ text_annotated$lemma2[l] ]


	t2 = text_annotated[text_annotated$upos %in% c("NOUN", "PROPN", "VERB", "SYM"), ]
	df_single = data.frame(
		action = "",
		lemma = tapply(t2$lemma2, t2$lemma2, unique),
		token = tapply(t2$token, t2$lemma2, function(x) paste(unique(x), collapse = ",")),
		upos = tapply(t2$upos, t2$lemma2, function(x) paste(unique(x), collapse = ",")),
		xpos = tapply(t2$xpos, t2$lemma2, function(x) paste(unique(x), collapse = ",")),
		n_docs = tapply(t2$doc_id, t2$lemma2, length)
	)

	abstracts[[i]] = df_single
}
names(abstracts) = xml_files


save(abstracts, file = "~/project/development/KeywordsEnrichment/pubmed_abstract.RData")

for(i in which(sapply(abstracts, is.null))) {
	f = xml_files[i]
	cat(f, "\n")
	
	oe = try({
		xml = read_xml(f, options = "RECOVER")
		ab <- xml_find_all(xml, "PubmedArticle/MedlineCitation/Article/Abstract/AbstractText") |> xml_text()
	})
	if(inherits(oe, "try-error")) {
		abstracts[[f]] = oe
	} else {
		abstracts[[f]] = ab
	}
}