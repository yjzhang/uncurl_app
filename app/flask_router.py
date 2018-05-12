import os

from flask import render_template

from app import app

app.user_dirs = []
app.test_dirs = []

def get_test_dirs(base='test_data'):
    subdirectories = os.listdir(base)
    good_dirs = []
    # allow dir to be accessed if uncurl has ran? if some visualization
    # is available?
    for s in subdirectories:
        if os.path.exists(os.path.join(base, 'mds_data.txt')):
            good_dirs.append(s)
    return good_dirs

def initialize():
    try:
        test_dirs = get_test_dirs('test_data')
    except:
        test_dirs = []
    try:
        user_dirs = get_test_dirs('/tmp/uncurl')
    except:
        user_dirs = []
    app.user_dirs = user_dirs
    app.test_dirs = test_dirs

@app.route('/data')
def data_index():
    initialize()
    return render_template('list_view.html',
            user_dirs=app.user_dirs,
            test_dirs=app.test_dirs)


if __name__ == '__main__':
    app.run(debug=True)
