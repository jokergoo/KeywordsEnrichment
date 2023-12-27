
setwd("~/project/development/KeywordsEnrichment/inst/extdata")

######################################################
## biomedical texts are from the following sources:
## - GO name
## - GO definition
## - RefSeq gene summary
## - UniProt keywords
## - pathway ontology
## - disease ontology
## - vertebrate trait ontology
## - human phenotype ontology
docs = NULL

library(simona)
dag = create_ontology_DAG_from_GO_db("BP")
docs = c(docs, mcols(dag)[, "definition"])
docs = c(docs, mcols(dag)[, "name"])

dag = create_ontology_DAG_from_GO_db("CC")
docs = c(docs, mcols(dag)[, "definition"])
docs = c(docs, mcols(dag)[, "name"])

dag = create_ontology_DAG_from_GO_db("MF")
docs = c(docs, mcols(dag)[, "definition"])
docs = c(docs, mcols(dag)[, "name"])

library(GeneSummary)
tb = loadGeneSummary(organism = 9606)

docs = c(docs, tapply(tb$Gene_summary, tb$Gene_ID, function(x) x[1]))


dag = ontology_kw()
docs = c(docs, mcols(dag)[, "description"])

dag = ontology_pw()
docs = c(docs, mcols(dag)[, "definition"])

dag = ontology_rdo()
docs = c(docs, mcols(dag)[, "definition"])

dag = ontology_vt()
docs = c(docs, mcols(dag)[, "definition"])

dag = ontology_hp()
docs = c(docs, mcols(dag)[, "definition"])


docs = docs[!is.na(docs)]
docs = docs[docs != ""]

length(docs)

# A-B is parsed as three tokens: A, -, B, so we removed "-" or "_"
docs = gsub("_|-|/", " ", docs)



##################################
## package udpipe is used to tokenize and lemmatize words in the documents

library(udpipe)

if(!exists("m_eng_ewt_loaded")) {
	m_eng_ewt = udpipe_download_model(language = "english")
	m_eng_ewt_path = m_eng_ewt$file_model
	m_eng_ewt_loaded = udpipe_load_model(file = m_eng_ewt_path)
}

text_annotated = udpipe_annotate(m_eng_ewt_loaded, x = docs, doc_id = seq_along(docs), parser = "none")
text_annotated = as.data.frame(text_annotated)

# remove NA document
text_annotated = text_annotated[text_annotated$sentence != "NA", ]

################################
## we will not use the lemmatized words from udpipe because many are wrongly lemmatized
## we will first normalize the original tokens by several rules
token = text_annotated$token

## step 1. remove punctuation as the first character or the last character
l = grepl("^[[:punct:]][[:alnum:]]{2,}", token)
token[l] = gsub("^[[:punct:]]", "", token[l])
l = grepl("[[:alnum:]]{2,}[[:punct:]]$", token)
token[l] = gsub("[[:punct:]]$", "", token[l])


## step 2. identify captitalized words
uncaptitalize = function(str) {
	substr(str, 0, 1) = tolower(substr(str, 0, 1))
	str
}

# a word contains vowel and no numbers
l = grepl("^[A-Z][a-z]+$", text_annotated$token) & 
    grepl("[oeuiaOEUIA]", text_annotated$token)
sw = uncaptitalize(token[l])

# and the word must be a "normal" english words (words::words contains a list of english words)
l2 = sw %in% words::words[, 1]
token[which(l)[l2]] = sw[l2]

# manually validated capitalized words
cp = readLines("capitalized_words.txt")
l = text_annotated$token %in% cp
token[l] = tolower(token[l])

# step 3. change plural to singulars
library(pluralize)
# possible plurals
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

# singular_rev_map = structure(w, names = lemmatize_words(singularize(w)))

# step 4. we use another package textstem for lemmatization because it seems it generates less wrong words than udpipe
text_annotated$lemma2 = textstem::lemmatize_words(token)

# # manually fix failed lemmatizations, most are upper-cased
# l = grepl("^[A-Z]+$", text_annotated$token) & 
#     grepl("[OEUIA]", text_annotated$token) & 
#     nchar(text_annotated$token) > 4
# tb = table(text_annotated$token[l])
# tb = tb[tb > 5]

# library(org.Hs.eg.db)
# gene_symbols = c(keys(org.Hs.eg.db, keytype = "SYMBOL"), keys(org.Hs.eg.db, keytype = "ALIAS"))

# tb = tb[!names(tb) %in% gene_symbols]
# cat(names(tb), sep = "\n")
# # then use google translate to validate whether they are true words


# step 5. manually curated words
lemma_map = read.table("adjust_words.txt")
lemma_map = structure(lemma_map[, 2], names = lemma_map[, 1])
l = text_annotated$lemma2 %in% names(lemma_map)
text_annotated$lemma2[l] = lemma_map[ text_annotated$lemma2[l] ]


# step 6. words may still be wrongly singularized by `lemmatize_words()`
l = grepl("sis$", token) | 
    grepl("'s$", token) | 
    grepl("tis$", text_annotated$token)
l = l | token %in% exclude_singular
l = l & text_annotated$upos %in% c("NOUN", "PROPN", "INTJ")
text_annotated$lemma2[l] = uncaptitalize(text_annotated$token[l])


# now text_annotated$lemma2 contains "normalized" and "lemmatized" words
l = text_annotated$lemma2 == ""
text_annotated$lemma2[l] = text_annotated$token[l]
saveRDS(text_annotated, file = "text_annotated.rds")

##########################################
## next we generate two lists
## - single words
## - phrases that contain multiple words

#### single words ####
t2 = text_annotated
t2 = t2[!t2$lemma2 %in% tm::stopwords(), ]
df_single = data.frame(
	action = "",
	keyword = tapply(t2$lemma2, t2$lemma2, unique),
	token = tapply(t2$token, t2$lemma2, function(x) paste(unique(x), collapse = ",")),
	upos = tapply(t2$upos, t2$lemma2, function(x) paste(unique(x), collapse = ",")),
	xpos = tapply(t2$xpos, t2$lemma2, function(x) paste(unique(x), collapse = ",")),
	n_docs = tapply(t2$doc_id, t2$lemma2, length)
)

df_single = df_single[df_single$n_docs >= 20, ]
df_single = df_single[nchar(df_single$keyword) > 1, ]

df_single = df_single[grep("[a-zA-Z]", df_single$keyword), ]



update_df_single = function(df_single) {
	if(file.exists("single_word.csv")) {
		df = read.csv("single_word.csv", header = TRUE)
		rownames(df) = df$keyword
		rownames(df_single) = df_single$keyword

		new_rn = setdiff(rownames(df_single), rownames(df))
		if(length(new_rn)) {
			df = rbind(df, df_single[new_rn, , drop = FALSE])

			write.csv(df, file = "single_word.csv", row.names = FALSE)
		}
	} else {
		write.csv(df_single, file = "single_word.csv", row.names = FALSE)
	}
}

update_df_single(df_single)

## then we manualy go through single_word.csv and do two things:
# 1. remove words that provide no useful information
# 2. merge words with the same stem. e.g. transcribe, transcription, transcriptional all to transcription


####### phrase #######

l = rep(TRUE, nrow(text_annotated))
l[text_annotated$lemma2 %in% tm::stopwords()] = FALSE

# the words excluded from the phrase
phrase_exclude = scan("phrase_exclude.txt", what = "character")
l = l & !text_annotated$lemma2 %in% c(phrase_exclude, '"')

# A phrase can not start with A or An
l = l & !(text_annotated$token_id == "1" & text_annotated$lemma2 %in% c("A", "An"))
t1 = text_annotated[l, ]

# phrases are looked for in a list of continous tokens,
# each row in `ir` is a range of indices of `t1`. E.g. a range of `[4, 6]` means token in row 4, 5, 6 are contious tokens in a sentence
t1$token_id = as.numeric(t1$token_id)
library(IRanges)
ir = tapply(seq_len(nrow(t1)), paste(t1$doc_id, t1$paragraph_id, t1$sentence_id, sep = "_"), function(ind) {
	cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", t1$doc_id[ind[1]])
	token_id = t1$token_id[ind]
	x = Rle(cumsum(c(1, diff(token_id))) - cumsum(rep(1, length(token_id))))
	IRanges(start = ind[start(x)], end = ind[end(x)])
})

names(ir) = dimnames(ir)
attributes(ir) = NULL
class(ir) = "list"

ir = unlist(as(ir, "IRangesList"))
ir = ir[width(ir) > 1]


## to identify high frequency phrases, we create a matrix where rows are preceding words and columns are succeeding words
## for example, in the phrase "cell cycle", cell is the preceding word and cycle is the succeeding word. We basically
## count the frequecy of "cell cycle" in the document set.
s = start(ir)
e = end(ir)
words = unique(unlist(lapply(seq_len(length(ir)), function(i) t1$lemma2[seq(s[i], e[i])])))

# the matrix is huge and sparse, here we use the TsparseMatrix class
library(Matrix)
pm1 = Matrix(data = 0, nrow = length(words), ncol = length(words), sparse = TRUE)
pm1 = as(pm1, "TsparseMatrix")
words_to_ind = structure(seq_along(words), names = words)

lemma_to_ind = words_to_ind[t1$lemma2]

# start(ir), end(ir), ind
library(Rcpp)
sourceCpp(code = "
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
sp_mat calc_freq(sp_mat pm, IntegerVector s, IntegerVector e, IntegerVector ind) {
	int n = s.size();
	for(int i = 0; i < n; i ++) {
		for(int j = s[i]; j < e[i]; j ++) {
			pm(ind[j-1]-1, ind[j]-1) ++;
		}
	}
	return(pm);
}

// [[Rcpp::export]]
IntegerVector leading_word_freq(int n, IntegerVector s, IntegerVector ind) {
	IntegerVector freq(n);
	for(int i = 0; i < s.size(); i ++) {
		freq[ ind[s[i]-1]-1 ] ++;
	}
	return(freq);
}
")


pm1 = calc_freq(pm1, start(ir), end(ir), unname(lemma_to_ind))
pm1 = as(pm1, "dgTMatrix")

# ldf = leading_word_freq(nrow(pm1), start(ir), unname(lemma_to_ind))
# names(ldf) = words

# saveRDS(pm1, file = "pm1.rds")

library(slam)
rs = row_sums(pm1)
names(rs) = names(words_to_ind)
cs = col_sums(pm1)
names(cs) = names(words_to_ind)


df = data.frame(
	action = "",
	i1 = pm1@i + 1, 
	i2 = pm1@j + 1,
	preceding = words[pm1@i+1], 
	succeeding = words[pm1@j+1],
	freq = pm1@x
)

nn = dim(pm1)[1]
df[["P(suc|pre)"]] = df$freq/rs[ (df$i1-1) %% nn + 1 ]
df[["P(pre|suc)"]] = df$freq/cs[ (df$i2-1) %% nn + 1 ]


upos = structure(text_annotated$upos, names = text_annotated$lemma2)
df$upos1 = upos[df$preceding]
df$upos2 = upos[df$succeeding]

# only the high frequent preceding-succeeding pairs
df = df[df$freq >= 20 & (df[["P(suc|pre)"]] > 0.2 | df[["P(pre|suc)"]] > 0.2), ]
df = df[order(-df$freq), ]

df_phrase = df

# update_df_phrase = function(df_phrase) {
# 	if(file.exists("two_words.csv")) {
# 		df = read.csv("two_words.csv", header = TRUE)
# 		rownames(df) = paste(df$preceding, df$succeeding, sep = "|")
# 		rownames(df_phrase) = paste(df_phrase$preceding, df_phrase$succeeding, sep = "|")

# 		new_rn = setdiff(rownames(df_phrase), rownames(df))
# 		if(length(new_rn)) {
# 			df = rbind(df, df_phrase[new_rn, , drop = FALSE])
# 			write.csv(df, file = "two_words.csv", row.names = FALSE)
# 		}
# 	} else {
# 		write.csv(df_phrase, file = "two_words.csv", row.names = FALSE)
# 	}
# }

# update_df_phrase(df_phrase)


## since `df_phrase` only contains relations of two words, to extend to phrases with more words,
## we go back to `ir`, for every set of continuous words (every row in `ir`), we obtain phrases
## based on `df_phrase` (remember `df_phrase` has already been filtered and only high frequent
## phrases are kept.
get_phrase = function(words) {

	if(length(words) <= 1) {
		return(character(0))
	}

	df = df_phrase[df_phrase$preceding %in% words & df_phrase$succeeding %in% words, , drop = FALSE]

	i = 1
	new_kw = character(0)
	kw_list = list()
	while(1) {
		if(words[i] %in% df$preceding) {
			if(length(new_kw) == 0) {
				if(words[i+1] %in% df$succeeding) {
					new_kw = paste(words[i], words[i+1], sep = " ")
				} 
			} else {
				if(words[i+1] %in% df$succeeding) {
					new_kw = paste(new_kw, words[i+1], sep = " ")
				} 
			}
		} else {
			if(length(new_kw) > 0) {
				kw_list[[length(kw_list) + 1]] = new_kw
				new_kw = character(0)
			}
		}
		i = i + 1
		if(i == length(words)) {
			if(length(new_kw) > 0) {
				kw_list[[length(kw_list) + 1]] = new_kw
				new_kw = character(0)
			}
			break
		}
	}

	unique(unlist(kw_list))
}

# get candicate phrase
s = start(ir)
e = end(ir)
phrase = list()
for(k in seq_along(s)) {
	cat(strrep("\b", 100))
	cat(k, "/", length(s))
	phrase[[k]] = get_phrase(t1$lemma2[ seq(s[k], e[k]) ])
}
cat("\n")



## only keep phrase
foo = sort(table(unlist(phrase)))
foo = foo[foo >= 20]
# the next two lines are optional. They basically reorder the phrases by the last words
ltt = strsplit(names(foo), " ")
foo = foo[order(sapply(ltt, function(x) x[length(x)]), sapply(ltt, function(x) x[length(x)-1]))]

write.csv(data.frame(phrase = names(foo), n = as.vector(foo)), file = "phrase_list.csv", row.names = FALSE)


## then we also manualy go through phrase_list.csv and do two things:
# 1. remove phrases that provide no useful information
# 2. reformat some phrases, e.g. "3 ' UTR" -> "3' UTR"






# ###

# df_single = df_single[df_single$action != "x", ]
# df_phrase = df_phrase[df_phrase$action != "x", ]

# l = !df_single$action %in% c("")
# map = structure(df_single$action[l], names = df_single$lemma[l])
# map = map[!duplicated(names(map))]
# new_nm = setdiff(names(map), names(lemma_map))
# lemma_map = c(lemma_map, map[new_nm])


# df_single$lemma2 = ifelse(df_single$action == "", df_single$lemma, df_single$action)


# # global objects:
# # - df_single
# # - df_phrase
# # - lemma_map
# get_keywords = function(docs, doc_id) {
# 	docs = gsub("_|-", " ", docs)

# 	docs[is.na(docs)] = ""

# 	## tokenize and tag
# 	m_eng_ewt = udpipe_download_model(language = "english-ewt")
# 	m_eng_ewt_path = m_eng_ewt$file_model
# 	m_eng_ewt_loaded = udpipe_load_model(file = m_eng_ewt_path)

# 	anno = udpipe_annotate(m_eng_ewt_loaded, x = docs, doc_id = doc_id, parser = "none")
# 	anno = as.data.frame(anno)
# 	anno = anno[, c("doc_id", "paragraph_id", "sentence_id", "token_id", "token", "lemma", "upos", "xpos")]

# 	xpos = c("CD", "FW", "NFP", "UH", "AFX", "LS", "JJ", "JJR", "JJS", "NN", "NNS", "NNP", "NNPS", "SYM", "VB", "VBD", "VBG", "VBN", "VBP", "VBZ")
# 	l = anno$xpos %in% xpos
# 	l[anno$lemma %in% anno$lemma[l]] = TRUE
# 	anno = anno[l, ]

# 	anno$token_id = as.numeric(anno$token_id)
# 	anno$lemma = lemmatize_words(anno$token)
# 	l = anno$lemma %in% names(lemma_map)

# 	anno$lemma[l] = lemma_map[ anno$lemma[l] ]

# 	lt = split(anno, anno$doc_id)
# 	lt2 = lapply(lt, function(term) {

# 		qqcat(strrep("\b", 100))
# 		qqcat("doc: @{term$doc_id[1]}...")

# 		v = tapply(seq_len(nrow(term)), paste(term$paragraph_id, term$sentence_id, sep = "_"), function(ind) {
# 			token_id = term$token_id[ind]
# 			x = Rle(cumsum(c(1, diff(token_id))) - cumsum(rep(1, length(token_id))))
# 			ir = IRanges(start = ind[start(x)], end = ind[end(x)])
# 			kw2 = character(0)
# 			if(length(ir) > 0) {
# 				s = start(ir)
# 				e = end(ir)
# 				for(k in seq_along(s)) {
# 					kw2 = c(kw2, get_phrase(term$lemma[ seq(s[k], e[k]) ]))
# 				}
# 			}

# 			unique(kw2)
# 		})
# 		unname(unique(unlist(v)))
# 	})
# 	cat("\n")
# 	all_docs = names(lt2)
# 	all_terms = unique(unlist(lt2))

# 	tdm = matrix(0L, nrow = length(all_terms), ncol = length(all_docs))
# 	rownames(tdm) = all_terms
# 	colnames(tdm) = all_docs
# 	for(i in seq_along(lt2)) {
# 		if(length(lt2[[i]])) {
# 			tdm[lt2[[i]], i] = 1L
# 		}
# 	}

# 	as(tdm, "dgTMatrix")
# }


# get_phrase = function(words) {

# 	if(length(words) == 1) {
# 		if(words %in% df_single$lemma2) {
# 			return(words)
# 		} else {
# 			return(character(0))
# 		}
# 	}

# 	df = df_phrase[df_phrase$preceding %in% words & df_phrase$succeeding %in% words, , drop = FALSE]

# 	i = 1
# 	new_kw = character(0)
# 	kw_list = list()
# 	while(1) {
# 		if(words[i] %in% df$preceding) {
# 			if(length(new_kw) == 0) {
# 				if(words[i+1] %in% df$succeeding) {
# 					new_kw = paste(words[i], words[i+1], sep = " ")
# 				} else {
# 					if(words[i] %in% df_single$lemma2) {
# 						kw_list[[length(kw_list) + 1]] = words[i]
# 					}
# 				}
# 			} else {
# 				if(words[i+1] %in% df$succeeding) {
# 					new_kw = paste(new_kw, words[i+1], sep = " ")
# 				} else {
# 					kw_list[[length(kw_list) + 1]] = new_kw
# 					new_kw = character(0)

# 					if(words[i] %in% df_single$lemma2) {
# 						kw_list[[length(kw_list) + 1]] = words[i]
# 					}
# 				}
# 			}
# 		} else {
# 			if(length(new_kw) > 0) {
# 				kw_list[[length(kw_list) + 1]] = new_kw
# 				new_kw = character(0)
# 			}
# 		}
# 		i = i + 1
# 		if(i == length(words)) {
# 			if(length(new_kw) > 0) {
# 				kw_list[[length(kw_list) + 1]] = new_kw
# 				new_kw = character(0)
# 			} else if(words[i] %in% df_single$lemma2) {
# 				kw_list[[length(kw_list) + 1]] = words[i]
# 			}
# 			break
# 		}
# 	}

# 	unique(unlist(kw_list))
# }


