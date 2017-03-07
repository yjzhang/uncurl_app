from multiprocessing import Process
import os
import random

from flask import render_template, request, redirect, send_from_directory
from werkzeug import secure_filename

import numpy as np
import uncurl

from app import app

def load_input_data():
    string_data = ''
    if 'textarea' in request.form:
        string_data = request.form['textarea'].split('\n')
        output_filename = ''.join(chr(random.randint(0,25)+ord('a')) for i in range(10)) + '.txt'
    elif 'fileinput' in request.files:
        f = request.files['fileinput']
        string_data = f.readlines()
        output_filename = secure_filename(f.filename)
    else:
        return 'error: no file provided'
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
    data, k, output_filename = load_input_data()
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
    data, k, output_filename = load_input_data()
    user_id = ''.join(chr(random.randint(0,25)+ord('a')) for i in range(10))
    P = Process(target=state_estimation_thread, args=(data, k, user_id))
    P.start()
    return redirect('/state_estimation/results/'+user_id)

@app.route('/state_estimation/results/<user_id>')
def state_estimation_result(user_id):
    # TODO: check if results have been completed...
    if os.path.isdir(os.path.join('/tmp/', user_id)):
        return render_template('state_estimation_user.html',
                user_id=user_id, has_result=True)
    else:
        return render_template('state_estimation_user.html',
                user_id=user_id, has_result=False)

@app.route('/state_estimation/results/<user_id>/<filename>')
def state_estimation_file(user_id, filename):
    path = os.path.join('/tmp/', user_id)
    return send_from_directory(path, filename)

def state_estimation_thread(data, k, user_id):
    """
    Uses a greenlet to estimate stuff...
    """
    M, W = uncurl.poisson_estimate_state(data, k, disp=False)
    path = os.path.join('/tmp/', user_id)
    os.mkdir(path)
    np.savetxt(os.path.join(path, 'm.txt'), M)
    np.savetxt(os.path.join(path, 'w.txt'), W)
