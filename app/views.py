from multiprocessing import Process
import os
import pickle
import uuid

from flask import Markup, render_template, request, redirect, send_from_directory, url_for
from werkzeug import secure_filename

import numpy as np
import uncurl

from app import app

import vis

def load_input_data():
    if 'fileinput' in request.files:
        f = request.files['fileinput']
        output_filename = secure_filename(f.filename)
    else:
        return
    k = int(request.form['k'])
    data = np.loadtxt(f)
    return data, k, output_filename

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/cluster')
def cluster():
    return render_template('clustering.html')

@app.route('/cluster/input', methods=['POST'])
def cluster_input():
    try:
        data, k, output_filename = load_input_data()
    except:
        return error('Error: no file found', 400)
    user_id = str(uuid.uuid4())
    P = Process(target=cluster_thread, args=(data, k, user_id))
    P.start()
    return redirect(url_for('clustering_result', user_id=user_id))

def cluster_thread(data, k, user_id, init=None):
    """
    Thread for performing the clustering operation - currently unused.
    """
    assignments, centers = uncurl.poisson_cluster(data, k, init)
    path = os.path.join('/tmp/', user_id)
    try:
        os.mkdir(path)
    except:
        pass
    with  open(os.path.join(path, 'assignments.txt'), 'w') as output_file:
        np.savetxt(output_file, assignments, fmt='%1.0f')
    with open(os.path.join(path, 'centers.txt'), 'w') as output_file:
        np.savetxt(output_file, centers)
    vis.vis_clustering(data, assignments, centers, user_id)

@app.route('/clustering/results/<user_id>')
def clustering_result(user_id):
    if os.path.exists(os.path.join('/tmp/', user_id, 'centers.txt')):
        try:
            visualization = open(os.path.join('/tmp/', user_id, 'vis_clustering.html')).read()
        except:
            visualization = ''
        visualization = Markup(visualization)
        return render_template('clustering_user.html',
                user_id=user_id, has_result=True,
                visualization=visualization)
    else:
        return render_template('clustering_user.html',
                user_id=user_id, has_result=False)

@app.route('/state_estimation')
def state_estimation():
    return render_template('state_estimation.html')

@app.route('/state_estimation/input', methods=['POST'])
def state_estimation_input():
    try:
        data, k, output_filename = load_input_data()
    except:
        return error('Error: no file found', 400)
    user_id = str(uuid.uuid4())
    P = Process(target=state_estimation_thread, args=(data, k, user_id))
    P.start()
    return redirect(url_for('state_estimation_result', user_id=user_id))

@app.route('/state_estimation/results/<user_id>')
def state_estimation_result(user_id):
    if os.path.exists(os.path.join('/tmp/', user_id, 'm.txt')):
        try:
            visualization = open(os.path.join('/tmp/', user_id, 'vis_state_estimation.html')).read()
        except:
            visualization = ''
        visualization = Markup(visualization)
        return render_template('state_estimation_user.html',
                user_id=user_id, has_result=True,
                visualization=visualization)
    else:
        return render_template('state_estimation_user.html',
                user_id=user_id, has_result=False)

@app.route('/<x>/results/<user_id>/<filename>')
def state_estimation_file(x, user_id, filename):
    path = os.path.join('/tmp/', user_id)
    return send_from_directory(path, filename)

def state_estimation_thread(data, k, user_id):
    """
    Uses a new process to do state estimation
    """
    path = os.path.join('/tmp/', user_id)
    os.mkdir(path)
    # if debugging, set max_iters to 1, inner_max_iters to 1... should be in config
    M, W = uncurl.poisson_estimate_state(data, k, max_iters=10, inner_max_iters=400, disp=False)
    np.savetxt(os.path.join(path, 'm.txt'), M)
    np.savetxt(os.path.join(path, 'w.txt'), W)
    vis.vis_state_estimation(data, M, W, user_id)

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
        if not os.path.exists(os.path.join('/tmp/', user_id, 'm.txt')):
            return error('Data for user id not found', 400)
        P = Process(target=lineage_thread, args=(None, None, None, None, user_id))
        P.start()
        return redirect(url_for('lineage_input_user_id', user_id=user_id))
    elif 'fileinput' in request.files:
        data, k, output_filename = load_input_data()
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
    path = os.path.join('/tmp/', user_id)
    if data is not None:
        # run state estimation then lineage estimation
        os.mkdir(path)
        M, W = uncurl.poisson_estimate_state(data, k, max_iters=10, inner_max_iters=400, disp=False)
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
    if os.path.exists(os.path.join('/tmp/', user_id, 'smoothed_data.txt')):
        try:
            visualization = open(os.path.join('/tmp/', user_id, 'vis_lineage.html')).read()
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
    centers = uncurl.qual2quant(data, qual)
    path = os.path.join('/tmp/', user_id)
    with open(os.join(path, 'qual2quant_centers.txt'), 'w') as f:
        np.savetxt(f, centers)

@app.route('/qual2quant/results/<user_id>')
def qual2quant_result(user_id):
    if os.path.exists(os.path.join('/tmp/', user_id, 'qual2quant_centers.txt')):
       return render_template('qual2quant_user.html',
                user_id=user_id, has_result=True,
                visualization=None)
    else:
        return render_template('qual2quant_user.html',
                user_id=user_id, has_result=False)


def error(msg, code):
    return render_template('error.html', msg=msg), code
