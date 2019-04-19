User Guide
==========

Uncurl-app is a tool for exploratory analysis of single-cell RNA-seq datasets.

## Integration with split-seq

Uncurl-app can directly work with the output of <https://github.com/yjzhang/split-seq-pipeline>. To do so after running the split-seq pipeline, go to the folder ending with `DGE_filtered`, and run the command `uncurl_app_split_seq \`pwd\``. Then go to http://127.0.0.1:5000/ in a web browser. This will automatically load the data for analysis.


## Loading data

The app can take in sparse or dense matrix formats. The formats that we accept are the outputs of, for example, `scipy.io.mmwrite` and `numpy.savetxt`, where the rows represent genes and the columns represent cells. Files should have one of the extensions `.mtx`, `.mtx.gz`, `.txt`, `.txt.gz`.

It also requires a list of gene symbols, where gene names are one per line.

## Options

After uploading the data, you will be redirected to the data preview screen. There are a number of options, most of which are set to some reasonable defaults.

The most important options are the number of cell types (k), and the visualization options. TSNE-based visualizations might take longer.

Other options include:
- Whether or not to exclude cells with low or high read counts. By default, this is set to exclude the lowest and highest 5% of cells by read count.
- How many cells to include in the visualization. More cells will take longer and might slow down the visualization. This is set to be a maximum of 10000 cells.

After selecting the options, click on Submit in order to start the preprocessing.

## Interaction

There are three primary views: the scatterplot, the barplot, and the database queries.

### Scatterplot view

This view shows a dimensionality-reduced view of the cells.

There are three different scatterplot views showing either the cells or identified cell types, selected using the radio buttons. "Cluster means" shows the cell archetypes. "Post-processed cells" shows the uncurl-based visualization, and "Pre-processed cells" shows a visualization that does not depend on uncurl. 

### Barplot view

Click on a cell on the scatterplot view to show the corresponding top genes for the cell's cluster. From the dropdown, it is possible to find the top genes in a number of different ways: 

select whether to view the top genes by c-score (ratio of expression to highest expression in another cluster), or the p-value calculated from this score.

### Database queries

There are three databases that can be used to identify cell types.

#### Enrichr

Click the "Submit Enrichr query" button to query the currently selected gene set in Enrichr.

#### CellMarker


#### CellMeSH

### Merge/split clusters

To merge or split clusters, use the box select tool, and select cells belonging to the desired clusters. The selected clusters will be shown below the scatterplot. Click the Merge Clusters or Split Clusters button to run the indicated operation. It might take several minutes before the operation completes.
