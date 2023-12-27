
library(simplifyEnrichment)

set.seed(123)
go_id = random_GO(500)
mat = GO_similarity(go_id)

cl = binary_cut(mat)



go = rownames(mat)[cl == 5]

tdm1 = readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
tb1 = simplifyEnrichment:::keyword_enrichment(go, tdm1, )



library(GO.db)
go_text = Term(GOTERM)
go_def = Definition(GOTERM)

tdm2 = docs_to_keywords(go_def, names(go_def), udpipe_model)
tb2 = simplifyEnrichment:::keyword_enrichment(go, tdm2)


dag = ontology_kw()
docs = mcols(dag)[, "description"]
tdm2 = docs_to_keywords(docs, dag@terms, udpipe_model)
