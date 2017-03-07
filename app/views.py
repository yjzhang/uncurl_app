from flask import render_template, request
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
def cluster_input(data):
    string_data = request.form['textarea'].split('\n')
    k = int(request.form['k'])
    line1 = string_data[0].split()
    data = np.zeros((len(string_data), len(line1)))
    for i, line in enumerate(string_data):
        data[i,:] = np.fromstring(line, sep='\t')
    assignments, centers = uncurl.poisson_cluster(data, k)

@app.route('/state_estimation')
def state_estimation():
    return render_template('state_estimation.html')

