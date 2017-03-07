import os
import random

from flask import render_template, request, send_from_directory
from werkzeug import secure_filename
import numpy as np
import uncurl
from app import app

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/cluster')
def cluster():
    return render_template('clustering.html')

@app.route('/cluster/input', methods=['POST'])
def cluster_input():
    string_data = ''
    if 'textarea' in request.form:
        string_data = request.form['textarea'].split('\n')
        output_filename = ''.join(chr(random.randint(0,25)+ord('a') for i in range(10))) + '.txt'
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
    assignments, centers = uncurl.poisson_cluster(data, k)
    with  open('/tmp/' + output_filename, 'w') as output_file:
        np.savetxt(output_file, assignments, fmt='%1.0f', newline=' ')
        output_file.write('\n')
        np.savetxt(output_file, centers)
    return send_from_directory('/tmp/', output_filename)

@app.route('/state_estimation')
def state_estimation():
    return render_template('state_estimation.html')

