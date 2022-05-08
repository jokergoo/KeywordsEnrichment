## Example

```{r}
load("sig_mat.RData")
expr = t(scale(t(sig_mat)))
km = kmeans(expr, centers = 3)$cluster

source("")
Heatmap(expr, row_split = km, 
	show_row_names = FALSE, show_row_dend = FALSE, 
	show_column_names = FALSE) + 
rowAnnotation(
	keywords = anno_GO_keywords(split = km, genes = rownames(expr), max_words = 30)
)
```
