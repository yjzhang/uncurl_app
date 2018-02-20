import os

import dash
from flask import Flask, render_template, redirect
from flask_bootstrap import Bootstrap

import dash_cluster_view

# TODO: map dash urls using flask?

flask = Flask(__name__)
Bootstrap(flask)
flask.dash_apps = {}
flask.user_dirs = []
flask.test_dirs = []

def get_test_dirs(base='test'):
    subdirectories = os.listdir(base)
    return subdirectories

def initialize():
    try:
        test_dirs = get_test_dirs('test')
    except:
        test_dirs = []
    try:
        user_dirs = get_test_dirs('/tmp/uncurl')
    except:
        user_dirs = []
    for d in test_dirs:
        flask.dash_apps[d] = deploy_dash_app(#os.path.join('test', d),
                '/test_dash/'+d)
    for d in user_dirs:
        flask.dash_apps[d] = deploy_dash_app(#os.path.join('user', d),
                '/user_dash/'+d)
    flask.user_dirs = user_dirs
    flask.test_dirs = test_dirs

@flask.route('/')
def index():
    return render_template('list_view.html',
            user_dirs=flask.user_dirs,
            test_dirs=flask.test_dirs)

def deploy_dash_app(url):
    dash_app = dash.Dash(name='Cluster view', sharing=True, server=flask,
            url_base_pathname=url)
    dash_cluster_view.initialize_layout(dash_app)
    #dash_cluster_view.initialize(dash_app, data_dir)
    return dash_app

@flask.route('/user/<user_id>')
def route_user(user_id):
    test_dir = os.path.join('/tmp/uncurl', user_id)
    if not flask.dash_apps[user_id].initialized:
        dash_cluster_view.initialize(flask.dash_apps[user_id], test_dir)
    return redirect('/user_dash/'+str(user_id))

@flask.route('/test/<test_id>')
def route_test(test_id):
    test_dir = os.path.join('test', test_id)
    if not flask.dash_apps[test_id].initialized:
        dash_cluster_view.initialize(flask.dash_apps[test_id], test_dir)
    return redirect('/test_dash/'+str(test_id))

initialize()

if __name__ == '__main__':
    flask.run(debug=True)
