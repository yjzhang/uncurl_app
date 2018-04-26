User Guide
==========

Uncurl-app is a tool for exploratory analysis of single-cell RNA-seq datasets.

## Loading data

The app can take in sparse or dense matrix formats. The formats that we accept are the outputs of, for example, `scipy.io.mmwrite` and `numpy.savetxt`, where the rows represent genes and the columns represent cells. Files should have one of the extensions `.mtx`, `.mtx.gz`, `.txt`, `.txt.gz`.

It also requires a list of gene names, where gene names are one per line.

## Options

After uploading the data, you will be redirected to the data preview screen. There are a number of options, most of which are set to some reasonable defaults.

The most important options are the number of cell types (k), and the visualization options. TSNE-based visualizations might take longer.

Other options include:
- Whether or not to exclude cells with low or high read counts. By default, this is set to exclude the lowest and highest 5% of cells by read count.
- How many cells to include in the visualization. More cells will take longer and might slow down the visualization. This is set to be a maximum of 1500 cells.

After selecting the options, click on Submit in order to start the preprocessing.

## Interaction

This is divided into three views: the scatterplot, the barplot, and the gene set queries.

### Scatterplot view

There are three different scatterplot views showing either the cells or identified cell types, selected using the radio buttons. "Cluster means" shows the cell archetypes. "Processed cells" shows the uncurl-based visualization, and "Unprocessed cells" shows a visualization that does not depend on uncurl. 

### Barplot view

Click on a cell on the scatterplot view to show the corresponding top genes. From the dropdown, it's possible to select whether to view the top genes by c-score (ratio of expression to highest expression in another cluster), or the p-value calculated from this score.

### Enrichr

Click the "Submit Enrichr query" button to query the currently selected gene set in Enrichr.

Sometimes, it will show that Enrichr timed out. If that happens, wait a few minutes before trying again.

### Merge/split clusters

To merge or split clusters, use the box select tool, and select cells belonging to the desired clusters. The selected clusters will be shown below the scatterplot. Click the Merge Clusters or Split Clusters button to run the indicated operation. It might take several minutes before the operation completes.
