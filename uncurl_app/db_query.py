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
    top_genes = [x.strip().upper() for x in request.form['top_genes'].split()]
    print('update_cellmarker:', top_genes)
    # gene_set is a string.
    db = request.form['database_select']
    try:
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
            species = request.form['cellmesh_species']
            import cellmesh
            result = []
            if test_type == 'hypergeom':
                result = cellmesh.hypergeometric_test(top_genes, species=species, return_header=True)
            elif test_type == 'norm_hypergeom':
                result = cellmesh.normed_hypergeometric_test(top_genes, return_header=True, species=species)
            elif test_type == 'prob':
                from cellmesh import prob_method
                result = prob_method.prob_test(top_genes, return_header=True, species=species)
            elif test_type == 'gsva':
                from cellmesh import gsva_ext_method
                result = gsva_ext_method.gsva_ext_test(top_genes, return_header=True, species=species)
            cell_types = [result[0]]
        elif db == 'cellmesh_anatomy':
            mesh_subset = request.form['anatomy_mesh_subset']
            if len(mesh_subset) > 1:
                mesh_subset = [x.strip() for x in mesh_subset.split(',')]
            else:
                mesh_subset = None
            species = request.form['anatomy_species']
            # TODO: validate mesh_subset
            import cellmesh
            result = cellmesh.hypergeometric_test(top_genes, return_header=True, db_dir=cellmesh.ANATOMY_DB_DIR,
                    cell_type_subset=mesh_subset, species=species)
            cell_types = [result[0]]
        elif db == 'go':
            from cellmesh import go_query
            top_genes = [x.capitalize() for x in top_genes]
            result = go_query.gene_set_query(top_genes, return_header=True)
            for r in result[1:]:
                r[3] = ', '.join(r[3])
            return json.dumps(result, cls=SimpleEncoder)
        for i in range(1, min(20, len(result))):
            ri = result[i]
            gene_pmids = []
            genes = ri[3]
            for g in genes:
                gene_pmids.append('{0}: {1}'.format(g, ', '.join(pmid_to_link(x) for x in ri[4][g])))
            cell_types.append((ri[0], ri[1], ri[2], ', '.join(ri[3]), ', '.join(gene_pmids)))
        return json.dumps(cell_types, cls=SimpleEncoder)
    except Exception as e:
        import traceback
        text = traceback.format_exc()
        print(text)
        return 'Error: ' + str(e)

@db_query.route('/db_query/cell_info', methods=['POST'])
def get_cell_info():
    """
    Returns all genes/pmids associated with a cell type, ranked in order of hits, that pass a certain threshold.
    """
    # TODO
    pass

@db_query.route('/db_query/get_mesh_tree')
def get_mesh_tree():
    """
    Returns a tree of all mesh terms.
    """
    # TODO
    import cellmesh
    tree, id_to_name = cellmesh.get_cellmesh_anatomy_tree()
    new_tree = {}
    def process_tree(key, children, new_parent):
        new_key = (key, id_to_name[key])
        new_node = {}
        new_parent[new_key] = new_node
        for k, v in children.items():
            process_tree(k, v, new_node)
    for k, v in tree.items():
        process_tree(k, v, new_tree)
    results_dict = {
            'tree': tree,
            'id_to_name': id_to_name
    }
    return json.dumps(results_dict, cls=SimpleEncoder)
