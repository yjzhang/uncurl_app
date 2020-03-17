import os

from flask import render_template, current_app, Blueprint

flask_router = Blueprint('flask_router', __name__, template_folder='templates')


flask_router.user_dirs = []
flask_router.test_dirs = []

def get_test_dirs(base='test_data'):
    subdirectories = os.listdir(base)
    good_dirs = []
    # allow dir to be accessed if uncurl has ran? if some visualization
    # is available?
    for s in subdirectories:
        if os.path.exists(os.path.join(base, s, 'mds_data.txt')):
            good_dirs.append(s)
    return good_dirs

def initialize():
    try:
        test_dirs = get_test_dirs(current_app.config['TEST_DATA_DIR'])
    except:
        test_dirs = []
    try:
        user_dirs = get_test_dirs(current_app.config['USER_DATA_DIR'])
    except:
        user_dirs = []
    flask_router.user_dirs = user_dirs
    flask_router.test_dirs = test_dirs

@flask_router.route('/data')
def data_index():
    initialize()
    if not current_app.config['SHOW_ALL_RESULTS']:
        flask_router.user_dirs = []
    return render_template('list_view.html',
            user_dirs=flask_router.user_dirs,
            test_dirs=flask_router.test_dirs)

