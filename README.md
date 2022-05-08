## Example

```r
load(url("https://github.com/jokergoo/KeywordsEnrichment/raw/master/sig_mat.RData"))
expr = t(scale(t(sig_mat)))
km = kmeans(expr, centers = 3)$cluster

source("https://raw.githubusercontent.com/jokergoo/KeywordsEnrichment/master/anno.R")
Heatmap(expr, row_split = km, 
    show_row_names = FALSE, show_row_dend = FALSE, 
    show_column_names = FALSE) + 
rowAnnotation(
    keywords = anno_GO_keywords(split = km, genes = rownames(expr), max_words = 30)
)
```

<img width="1063" alt="image" src="https://user-images.githubusercontent.com/449218/167316913-42c51641-b40e-4694-962a-b9c4054dabd8.png">
