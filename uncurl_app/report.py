# generates a static report page
# TODO
# what should the report contain? A description of all the cell types...
# 1. scatterplot, both pre and post-processed, generated using plotly (this should include the js files directly)
# 2. list of top 50 genes per cluster using 1-vs-rest ratio
# 3. list of top 10 cell types per cluster using prob method on cellmesh
# 4. basic information about each cluster: num cells, mean/median read counts, mean/median gene counts
import json

from flask import render_template

from .utils import SimpleEncoder
from .interaction_views import interaction_views, get_sca, get_sca_gene_names, get_sca_top_1vr, scatterplot_data, update_cellmesh_result

@interaction_views.route('/user/<user_id>/report')
def generate_report(user_id):
    """
    Returns a rendered template for a report identifying the cell types
    for each cluster.
    """
    sca = get_sca(user_id)
    top_genes = get_sca_top_1vr(user_id)
    gene_names = get_sca_gene_names(user_id)
    cellmesh_results_clusters = {}
    cluster_top_genes = {}
    cluster_cell_counts = {}
    cluster_mean_reads = {}
    # scatterplot data
    scatterplot = scatterplot_data(sca.dim_red, sca.labels)
    # for each cluster
    for cluster_id in top_genes.keys():
        cluster_cells = (sca.labels == cluster_id)
        cluster_cell_counts[cluster_id] = cluster_cells.sum()
        cluster_mean_reads[cluster_id] = sca.read_counts[sca.cell_subset][sca.cell_sample][cluster_cells].mean()
        top_50_genes = top_genes[cluster_id][:50]
        cluster_top_genes[cluster_id] = top_50_genes
        selected_gene_names = [gene_names[x[0]] for x in top_50_genes]
        # do a cellmesh query
        cellmesh_results = update_cellmesh_result(user_id, selected_gene_names, 'prob', return_json=False)
        cellmesh_results_clusters[cluster_id] = cellmesh_results
    return render_template('report.html',
            user_id=user_id,
            results=cellmesh_results_clusters, #json.dumps(cellmesh_results_clusters, cls=SimpleEncoder),
            top_genes=cluster_top_genes,
            mean_reads=cluster_mean_reads,
            cell_counts=cluster_cell_counts,
            scatterplot=scatterplot,
            cluster_ids=sorted(top_genes.keys()))

