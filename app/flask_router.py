import os

from flask import Flask, render_template, redirect
from flask_bootstrap import Bootstrap


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
    #for d in test_dirs:
    #    if d not in app.dash_apps:
    #        app.dash_apps[d] = deploy_dash_app(#os.path.join('test', d),
    #                '/test_dash/'+d)
    #for d in user_dirs:
    #    if d not in app.dash_apps:
    #        app.dash_apps[d] = deploy_dash_app(#os.path.join('user', d),
    #                '/user_dash/'+d)
    app.user_dirs = user_dirs
    app.test_dirs = test_dirs

@app.route('/data')
def data_index():
    initialize()
    return render_template('list_view.html',
            user_dirs=app.user_dirs,
            test_dirs=app.test_dirs)

initialize()

if __name__ == '__main__':
    app.run(debug=True)
