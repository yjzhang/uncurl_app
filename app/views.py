from multiprocessing import Process
import os
import uuid

from flask import Markup, render_template, request, redirect, send_from_directory, url_for
from werkzeug import secure_filename

import numpy as np
import uncurl

from app import app

import vis

def load_input_data():
    string_data = ''
    if 'textarea' in request.form:
        string_data = request.form['textarea'].split('\n')
        output_filename = str(uuid.uuid4()) + '.txt'
    elif 'fileinput' in request.files:
        f = request.files['fileinput']
        string_data = f.readlines()
        output_filename = secure_filename(f.filename)
    else:
        return
    k = int(request.form['k'])
    line1 = string_data[0].split()
    data = np.zeros((len(string_data), len(line1)))
    for i, line in enumerate(string_data):
        #print line
        data[i,:] = np.fromstring(line, sep='\t')
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
    with  open(os.path.join('/tmp/', output_filename), 'w') as output_file:
        np.savetxt(output_file, assignments, fmt='%1.0f', newline=' ')
        output_file.write('\n')
        np.savetxt(output_file, centers)
    return send_from_directory('/tmp/', output_filename)

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
    string_data = ''
    if 'useridinput' in request.form:
        user_id = request.form['useridinput']
        # TODO: try to load m/w data, or return an error if you can't
    elif 'fileinput' in request.files:
        f = request.files['fileinput']
        string_data = f.readlines()
        output_filename = secure_filename(f.filename)
    elif 'mfileinput' in request.files and 'wfileinput' in request.files:
        pass
    else:
        return error('Missing data input', 400)
 

@app.route('/lineage/input/<user_id>')
def lineage_input_user_id(user_id):
    """
    Lineage input from a user's perspective
    """
    m = np.loadtxt(os.path.join('/tmp/', user_id, 'm.txt'))
    w = np.loadtxt(os.path.join('/tmp/', user_id, 'w.txt'))
    # TODO

def lineage_thread(m, w, user_id):
    pass

def error(msg, code):
    return render_template('error.html', msg=msg), code
