# db queries without running uncurl

import json

from flask import request, render_template, Blueprint

from .utils import SimpleEncoder

db_query = Blueprint('db_query', __name__,
        template_folder='templates')

def pmid_to_link(pmid):
    return '<a href="https://www.ncbi.nlm.nih.gov/pubmed/{0}">{0}</a>'.format(pmid)


@db_query.route('/db_query')
def db_query_index():
    return render_template('db_query_index.html')

@db_query.route('/db_query/submit', methods=['POST'])
def db_query_submit():
    top_genes = [x.strip().upper() for x in request.form['top_genes'].split('\n')]
    print('update_cellmarker:', top_genes)
    # gene_set is a string.
    db = request.form['database_select']
    if db == 'cellmarker':
        test_type = request.form['test_type']
        cells_or_tissues = request.form['cells_or_tissues']
        species = request.form['species']
        import cellmarker
        result = []
        if test_type == 'hypergeom':
            result = cellmarker.hypergeometric_test(top_genes, cells_or_tissues, return_header=True, species=species, return_cl=True)
        cell_types = [result[0]]
    elif db == 'cellmesh':
        test_type = request.form['mesh_test_type']
        import cellmesh
        result = []
        if test_type == 'hypergeom':
            result = cellmesh.hypergeometric_test(top_genes, return_header=True)
        elif test_type == 'norm_hypergeom':
            result = cellmesh.normed_hypergeometric_test(top_genes, return_header=True)
        elif test_type == 'prob':
            from cellmesh import prob_method
            result = prob_method.prob_test(top_genes, return_header=True)
        elif test_type == 'gsva':
            from cellmesh import gsva_ext_method
            result = gsva_ext_method.gsva_ext_test(top_genes, return_header=True)
        cell_types = [result[0]]
    for i in range(1, min(20, len(result))):
        ri = result[i]
        gene_pmids = []
        genes = ri[3]
        for g in genes:
            gene_pmids.append('{0}: {1}'.format(g, ', '.join(pmid_to_link(x) for x in ri[4][g])))
        cell_types.append((ri[0], ri[1], ri[2], ', '.join(ri[3]), ', '.join(gene_pmids)))
    return json.dumps(cell_types, cls=SimpleEncoder)

@db_query.route('/db_query/cell_info', methods=['POST'])
def get_cell_info():
    """
    Returns all genes/pmids associated with a cell type, ranked in order of hits, that pass a certain threshold.
    """
    # TODO
    pass

