import json
from multiprocessing import Process
import os
import pickle
import uuid

import bokeh
from flask import Markup, render_template, request, redirect, send_from_directory, url_for
from werkzeug import secure_filename

import numpy as np
import scipy.io
from scipy import sparse
import uncurl

from app import app

from . import vis
from .generate_analysis import generate_uncurl_analysis
from .data_stats import Summary

def load_upload_data(path=None):
    if 'fileinput' in request.files:
        f = request.files['fileinput']
        output_filename = secure_filename(f.filename)
        if output_filename == '':
            return
    else:
        return
    init_f = None
    if 'startinput' in request.files:
        init_f = request.files['startinput']
    # allow for mtx input data
    input_type = request.form['inputtype']
    if input_type == 'dense':
        data_filename = 'data.txt'
        if output_filename.endswith('.gz'):
            data_filename = 'data.txt.gz'
        f.save(os.path.join(path, data_filename))
        data = np.loadtxt(os.path.join(path, data_filename))
    elif input_type == 'sparse':
        data_filename = 'data.mtx'
        if output_filename.endswith('.mtx.gz'):
            data_filename = 'data.mtx.gz'
        f.save(os.path.join(path, data_filename))
        data = scipy.io.mmread(os.path.join(path, data_filename))
    init = None
    if init_f is not None and init_f.filename != '':
        init = np.loadtxt(init_f)
        init_f.save(os.path.join(path, 'init.txt'))
    return data, output_filename, init


def load_gene_names(path=None):
    if 'genenames' in request.files:
        f = request.files['genenames']
        if path is not None:
            f.save(os.path.join(path, 'gene_names.txt'))
    else:
        return None
    gene_names = np.loadtxt(f, dtype=str)
    print(gene_names)
    return gene_names

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/cluster')
def cluster():
    return render_template('clustering.html')

@app.route('/cluster/input', methods=['POST'])
def cluster_input():
    try:
        data, output_filename, init = load_upload_data()
    except:
        return error('Error: no file found', 400)
    k = int(request.form['k'])
    user_id = str(uuid.uuid4())
    dist_type = request.form['disttype']
    P = Process(target=cluster_thread, args=(data, k, user_id, init, dist_type))
    P.start()
    return redirect(url_for('clustering_result', user_id=user_id))

def cluster_thread(data, k, user_id, init=None, dist_type='Poisson'):
    """
    Thread for performing the clustering operation - currently unused.
    """
    path = os.path.join('/tmp/uncurl/', user_id)
    try:
        os.mkdir(path)
        print path
    except:
        pass
    if dist_type=='Poisson':
        assignments, centers = uncurl.poisson_cluster(data, k, init)
        with open(os.path.join(path, 'centers.txt'), 'w') as output_file:
            np.savetxt(output_file, centers)
        vis.vis_clustering(data, assignments, user_id)
    elif dist_type=='Negative binomial':
        assignments, P, R = uncurl.nb_cluster(data, k)
        with open(os.path.join(path, 'P.txt'), 'w') as output_file:
            np.savetxt(output_file, P)
        with open(os.path.join(path, 'R.txt'), 'w') as output_file:
            np.savetxt(output_file, R)
        vis.vis_clustering(data, assignments, user_id)
    elif dist_type=='Zero-inflated Poisson':
        assignments, L, M = uncurl.zip_cluster(data, k)
        with open(os.path.join(path, 'L.txt'), 'w') as output_file:
            np.savetxt(output_file, L)
        with open(os.path.join(path, 'M.txt'), 'w') as output_file:
            np.savetxt(output_file, M)
        vis.vis_clustering(data, assignments, user_id)
    with open(os.path.join(path, 'assignments.txt'), 'w') as output_file:
        np.savetxt(output_file, assignments, fmt='%1.0f')

@app.route('/clustering/results/<user_id>')
def clustering_result(user_id):
    if os.path.exists(os.path.join('/tmp/uncurl/', user_id, 'assignments.txt')):
        try:
            visualization = open(os.path.join('/tmp/uncurl/', user_id, 'vis_clustering.html')).read()
        except:
            visualization = ''
        poisson = True
        if os.path.exists(os.path.join('/tmp/uncurl/', user_id, 'centers.txt')):
            pass
        else:
            poisson=False
        visualization = Markup(visualization)
        return render_template('clustering_user.html',
                user_id=user_id, has_result=True,
                visualization=visualization, poisson=poisson)
    else:
        return render_template('clustering_user.html',
                user_id=user_id, has_result=False)

@app.route('/state_estimation')
def state_estimation():
    return render_template('state_estimation.html')

@app.route('/state_estimation/input', methods=['POST'])
def state_estimation_input():
    user_id = str(uuid.uuid4())
    path = os.path.join('/tmp/uncurl/', user_id)
    os.makedirs(path)
    # save state estimation params? save request.form
    with open(os.path.join(path, 'inputs.json'), 'w') as f:
        f.write(json.dumps(request.form))
    # TODO: if file is large, start a new thread. otherwise just
    # run the thing
    P = Process(target=state_estimation_preproc, args=(user_id, path))
    P.start()
    #state_estimation_preproc(user_id, path)
    return redirect(url_for('state_estimation_result', user_id=user_id))

@app.route('/state_estimation/results/<user_id>/start', methods=['POST'])
def state_estimation_start(user_id):
    """
    Actually start the process of state estimation.

    This saves a file called 'params.json' in /tmp/uncurl/<user_id>
    containing all parameters used in state estimation.
    """
    path = os.path.join('/tmp/uncurl/', user_id)
    gene_names = os.path.join(path, 'gene_names.txt')
    if not os.path.exists(gene_names):
        gene_names = None
    init = os.path.join(path, 'init.txt')
    if not os.path.exists(init):
        init = None
    # load json params
    with open(os.path.join(path, 'preprocess.json')) as f:
        preprocess = json.load(f)
    for key in request.form.keys():
        preprocess[key] = request.form[key]
    # params.json contains all input parameters to the state estimation
    with open(os.path.join(path, 'params.json'), 'w') as f:
        json.dump(preprocess, f)
    P = Process(target=state_estimation_thread, args=(user_id, gene_names, init, path, preprocess))
    P.start()
    return redirect(url_for('state_estimation_result', user_id=user_id))


@app.route('/state_estimation/results/<user_id>/')
def state_estimation_result(user_id):
    path = os.path.join('/tmp/uncurl/', user_id)
    if os.path.exists(os.path.join(path, 'mds_data.txt')):
        #try:
        #    visualization = open(os.path.join('/tmp/uncurl/', user_id, 'vis_state_estimation.html')).read()
        #except:
        #    visualization = ''
        #visualization = Markup(visualization)
        return redirect(url_for('route_user', user_id=user_id))
    elif os.path.exists(os.path.join(path, 'preprocess.json')):
        uncurl_is_running = os.path.exists(os.path.join(path, 'submitted'))
        with open(os.path.join(path, 'preprocess.json')) as f:
            preprocess = json.load(f)
        with open(os.path.join(path, 'vis_summary.html')) as f:
            v = f.read()
        v = Markup(v)
        return render_template('state_estimation_user.html',
                user_id=user_id, has_preview=True,
                uncurl_is_running=uncurl_is_running,
                visualization=v,
                min_reads=preprocess['min_reads'],
                max_reads=preprocess['max_reads'],
                cells=preprocess['cells'],
                genes=preprocess['genes'])
    else:
        return render_template('state_estimation_user.html',
                user_id=user_id, uncurl_is_running=False,
                has_result=False)

@app.route('/<x>/results/<user_id>/<filename>')
def state_estimation_file(x, user_id, filename):
    if x!='test':
        path = os.path.join('/tmp/uncurl/', user_id)
    else:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                'test_data', user_id)
    print(path)
    return send_from_directory(path, filename)

@app.route('/<x>/results/<user_id>/data_download')
def data_download(x, user_id):
    if x!='test':
        path = os.path.join('/tmp/uncurl/', user_id)
    else:
        path = os.path.join('test_data', user_id)
    files = os.listdir(path)
    files.sort()
    return render_template('data_download.html',
            user_id=user_id,
            test_or_user=x,
            files=files)

def state_estimation_preproc(user_id, path=None):
    """
    Preprocessing for state estimation - generates summary statistics,
    etc...
    """
    #try:
    data, output_filename, init = load_upload_data(path)
    is_gz = False
    if output_filename.endswith('.gz'):
        is_gz = True
    #except:
    #    return error('Error: no file found', 400)
    if path is None:
        path = os.path.join('/tmp/uncurl/', user_id)
    gene_names = load_gene_names(path)
    summary = Summary(data, path, is_gz)
    script, div = summary.visualize()
    summary.preprocessing_params()

def state_estimation_thread(user_id, gene_names=None, init=None, path=None, preprocess=None):
    """
    Uses a new process to do state estimation. Assumes that the input data is already saved in a directory named /tmp/uncurl/<user_id>/.

    Args:
        user_id (str)
        gene_names (str or array, optional): path to a list of gene names, or an array of gene names. Default: None
        init (array, optional): init_means- shape is (genes, k). Default: None.
        path (str, optional): Path where data and results are saved.
        preprocess (dict): dict containing additional parameters: min_reads, max_reads, normalize, is_sparse, is_gz, disttype, genes_frac, cell_frac, vismethod, baseline_vismethod
    """
    if path is None:
        path = os.path.join('/tmp/uncurl/', user_id)
    # get correct data names
    data = os.path.join(path, 'data.mtx')
    if preprocess['is_gz']:
        data += '.gz'
    if not os.path.exists(data):
        data = os.path.join(path, 'data.txt')
        if preprocess['is_gz']:
            data += '.gz'
    # TODO: it's really confusing where the param names come from in
    # the preprocess dict - they come from the input ids in
    # state_estimation_user.html.
    dist_type = preprocess['disttype']
    if dist_type=='Poisson':
        pass
    elif dist_type=='Negative binomial':
        dist_type = 'nb'
    elif dist_type == 'Log-Normal':
        dist_type = 'lognorm'
    uncurl_args = app.config['UNCURL_ARGS']
    uncurl_args['dist'] = dist_type
    uncurl_args['init_means'] = init
    k = int(preprocess['k'])
    vismethod = preprocess['vismethod']
    gene_frac = float(preprocess['genes_frac'])
    min_reads = int(preprocess['min_reads'])
    max_reads = int(preprocess['max_reads'])
    cell_frac = float(preprocess['cell_frac'])
    if 'normalize' in preprocess:
        normalize = True
    baseline_vismethod = preprocess['baseline_vismethod']
    # TODO: save params as json instead of having to pass them...
    # actually save params.json and pass params lol
    generate_uncurl_analysis(data, path, clusters=k, gene_names=gene_names,
            gene_sub=True,
            dim_red_option=vismethod,
            min_reads=min_reads,
            max_reads=max_reads,
            frac=gene_frac,
            cell_frac=cell_frac,
            normalize=normalize,
            baseline_dim_red=baseline_vismethod,
            **uncurl_args)

@app.route('/lineage')
def lineage():
    return render_template('lineage.html')

@app.route('/lineage/input', methods=['POST'])
def lineage_input():
    """
    Note: how do we implement this? is it a view from the state estimation folder? do we have a previous
    """
    if 'useridinput' in request.form and request.form['useridinput']:
        user_id = request.form['useridinput']
        # Try to load m/w data, or return an error if you can't
        if not os.path.exists(os.path.join('/tmp/uncurl/', user_id, 'm.txt')):
            return error('Data for user id not found', 400)
        P = Process(target=lineage_thread, args=(None, None, None, None, user_id))
        P.start()
        return redirect(url_for('lineage_input_user_id', user_id=user_id))
    elif 'fileinput' in request.files:
        data, output_filename, init = load_upload_data()
        k = request.form['k']
        user_id = str(uuid.uuid4())
        P = Process(target=lineage_thread, args=(data, k, None, None, user_id))
        P.start()
        return redirect(url_for('lineage_input_user_id', user_id=user_id))
    elif 'mfileinput' in request.files and 'wfileinput' in request.files:
        m_file = request.files['mfileinput']
        w_file = request.files['wfileinput']
        M = np.loadtxt(m_file)
        W = np.loadtxt(w_file)
        user_id = str(uuid.uuid4())
        P = Process(target=lineage_thread, args=(None, None, M, W, user_id))
        P.start()
        return redirect(url_for('lineage_input_user_id', user_id=user_id))
    else:
        return error('Missing data input', 400)

def lineage_thread(data, k, M, W, user_id):
    """
    Thread to do lineage calculation...
    """
    path = os.path.join('/tmp/uncurl/', user_id)
    if data is not None:
        # run state estimation then lineage estimation
        os.mkdir(path)
        M, W, ll = uncurl.poisson_estimate_state(data, k, max_iters=10, inner_max_iters=400, disp=False)
        np.savetxt(os.path.join(path, 'm.txt'), M)
        np.savetxt(os.path.join(path, 'w.txt'), W)
        curve_params, smoothed_data, edges, clusters = uncurl.lineage(M, W)
    elif M is not None and W is not None:
        # M and W provided as arguments
        os.mkdir(path)
        np.savetxt(os.path.join(path, 'm.txt'), M)
        np.savetxt(os.path.join(path, 'w.txt'), W)
        curve_params, smoothed_data, edges, clusters = uncurl.lineage(M, W)
    else:
        # try to read M, W from user_id
        M = np.loadtxt(os.path.join(path, 'm.txt'))
        W = np.loadtxt(os.path.join(path, 'w.txt'))
        curve_params, smoothed_data, edges, clusters = uncurl.lineage(M, W)
    # TODO: save curve params
    np.savetxt(os.path.join(path, 'smoothed_data.txt'), smoothed_data)
    with open(os.path.join(path, 'edges.txt'), 'w') as f:
        f.write(repr(edges))
    with open(os.path.join(path, 'clusters.txt'), 'w') as f:
        f.write(repr(clusters))
    vis.vis_lineage(M, W, smoothed_data, edges, clusters, user_id)

@app.route('/lineage/results/<user_id>')
def lineage_input_user_id(user_id):
    """
    Lineage input from a user's perspective
    """
    if os.path.exists(os.path.join('/tmp/uncurl/', user_id, 'smoothed_data.txt')):
        try:
            visualization = open(os.path.join('/tmp/uncurl/', user_id, 'vis_lineage.html')).read()
        except:
            visualization = ''
        visualization = Markup(visualization)
        return render_template('lineage_user.html',
                user_id=user_id, has_result=True,
                visualization=visualization)
    else:
        return render_template('lineage_user.html',
                user_id=user_id, has_result=False)

@app.route('/qual2quant')
def qual2quant():
    return render_template('qual2quant.html')

@app.route('/qual2quant/input')
def qual2quant_input():
    if 'fileinput' not in request.files or 'qualinput' not in request.files:
        return error('Missing data input', 400)
    cell_file = request.files['fileinput']
    qual_file = request.files['qualinput']
    cell_data = np.loadtxt(cell_file)
    qual_data = np.loadtxt(qual_file)
    user_id = str(uuid.uuid4())
    P = Process(target=qual2quant_thread, args=(cell_data, qual_data, user_id))
    P.start()
    return redirect(url_for('qual2quant_result', user_id=user_id))

def qual2quant_thread(data, qual, user_id):
    centers = uncurl.qualNorm(data, qual)
    path = os.path.join('/tmp/uncurl/', user_id)
    with open(os.join(path, 'qual2quant_centers.txt'), 'w') as f:
        np.savetxt(f, centers)

@app.route('/qual2quant/results/<user_id>')
def qual2quant_result(user_id):
    if os.path.exists(os.path.join('/tmp/uncurl/', user_id, 'qual2quant_centers.txt')):
       return render_template('qual2quant_user.html',
                user_id=user_id, has_result=True,
                visualization=None)
    else:
        return render_template('qual2quant_user.html',
                user_id=user_id, has_result=False)


def error(msg, code):
    return render_template('error.html', msg=msg), code
