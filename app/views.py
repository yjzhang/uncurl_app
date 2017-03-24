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
    assignments, centers = uncurl.poisson_cluster(data, k)
    user_id = str(uuid.uuid4())
    with  open(os.path.join('/tmp/', output_filename), 'w') as output_file:
        np.savetxt(output_file, assignments, fmt='%1.0f', newline=' ')
        output_file.write('\n')
        np.savetxt(output_file, centers)
    return send_from_directory('/tmp/', output_filename)

def cluster_thread(data, k, user_id, output_filename):
    """
    Thread for performing the clustering operation - currently unused.
    """
    assignments, centers = uncurl.poisson_cluster(data, k)
    with  open(os.path.join('/tmp/', user_id, output_filename), 'w') as output_file:
        np.savetxt(output_file, assignments, fmt='%1.0f', newline=' ')
        output_file.write('\n')
        np.savetxt(output_file, centers)

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
    # TODO: check if results have been completed...
    if os.path.isdir(os.path.join('/tmp/', user_id)):
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

@app.route('/state_estimation/results/<user_id>/<filename>')
def state_estimation_file(user_id, filename):
    path = os.path.join('/tmp/', user_id)
    return send_from_directory(path, filename)

def state_estimation_thread(data, k, user_id):
    """
    Uses a new process to do state estimation
    """
    # if debugging, set max_iters to 1, inner_max_iters to 1... should be in config
    M, W = uncurl.poisson_estimate_state(data, k, max_iters=10, inner_max_iters=400, disp=False)
    path = os.path.join('/tmp/', user_id)
    os.mkdir(path)
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
    if 'useridinput' in request.form:
        user_id = request.form['useridinput']
        # Try to load m/w data, or return an error if you can't
        if not os.path.exists(os.path.join('/tmp/', user_id, 'm.txt')):
            return error('Data for user id not found', 400)
        P = Process(target=lineage_thread, args=(None, None, None, None, user_id))
        P.run()
        return redirect(url_for('lineage_input_user_id', user_id=user_id))
    elif 'fileinput' in request.files:
        data, k = load_input_data()
        user_id = str(uuid.uuid4())
    elif 'mfileinput' in request.files and 'wfileinput' in request.files:
        m_file = request.files['mfileinput']
        w_file = request.files['wfileinput']
        M = np.loadtxt(m_file)
        W = np.loadtxt(w_file)
        user_id = str(uuid.uuid4())
        P = Process(target=lineage_thread, args=(None, None, M, W, user_id))
        P.run()
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
    with open(os.path.join(path, 'edges.pkl'), 'w') as f:
        pickle.dump(edges, f)
    vis.vis_lineage(M, W, smoothed_data, edges, clusters, user_id)

@app.route('/lineage/results/<user_id>')
def lineage_input_user_id(user_id):
    """
    Lineage input from a user's perspective
    """
    # TODO
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


def error(msg, code):
    return render_template('error.html', msg=msg), code
