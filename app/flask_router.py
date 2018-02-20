import os

import dash
from flask import Flask, render_template, redirect
from flask_bootstrap import Bootstrap

import dash_cluster_view

from app import app

# TODO: map dash urls using flask?

app.dash_apps = {}
app.user_dirs = []
app.test_dirs = []

def get_test_dirs(base='test_data'):
    subdirectories = os.listdir(base)
    return subdirectories

def initialize():
    try:
        test_dirs = get_test_dirs('test_data')
    except:
        test_dirs = []
    try:
        user_dirs = get_test_dirs('/tmp/uncurl')
    except:
        user_dirs = []
    for d in test_dirs:
        app.dash_apps[d] = deploy_dash_app(#os.path.join('test', d),
                '/test_dash/'+d)
    for d in user_dirs:
        app.dash_apps[d] = deploy_dash_app(#os.path.join('user', d),
                '/user_dash/'+d)
    app.user_dirs = user_dirs
    app.test_dirs = test_dirs

@app.route('/data')
def data_index():
    return render_template('list_view.html',
            user_dirs=app.user_dirs,
            test_dirs=app.test_dirs)

def deploy_dash_app(url):
    dash_app = dash.Dash(name='Cluster view', sharing=True, server=app,
            url_base_pathname=url)
    dash_cluster_view.initialize_layout(dash_app)
    #dash_cluster_view.initialize(dash_app, data_dir)
    return dash_app

@app.route('/user/<user_id>')
def route_user(user_id):
    test_dir = os.path.join('/tmp/uncurl', user_id)
    if not app.dash_apps[user_id].initialized:
        dash_cluster_view.initialize(app.dash_apps[user_id], test_dir,
                '/user/'+user_id)
    return redirect('/user_dash/'+str(user_id))

@app.route('/test/<test_id>')
def route_test(test_id):
    test_dir = os.path.join('test_data', test_id)
    if not app.dash_apps[test_id].initialized:
        dash_cluster_view.initialize(app.dash_apps[test_id], test_dir,
                '/test/'+test_id)
    return redirect('/test_dash/'+str(test_id))

initialize()

if __name__ == '__main__':
    app.run(debug=True)
