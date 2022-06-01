import os
from flask import Flask, render_template
from flask_bootstrap import Bootstrap

from . import interaction_views, views, flask_router, db_query, report

from .cache import cache

def create_app(config_filename=None):
    app = Flask(__name__)
    Bootstrap(app)

    # maximum file length is 250MB
    if 'MAX_CONTENT_LENGTH' in os.environ:
        app.config['MAX_CONTENT_LENGTH'] = int(os.environ['MAX_CONTENT_LENGTH'])
    else:
        app.config['MAX_CONTENT_LENGTH'] = 250 * 1024 * 1024
    # default args to pass to uncurl.run_state_estimation
    app.config['UNCURL_ARGS'] = {
            'threads': 2,
            'max_iters': 20,
            'inner_max_iters': 50
    }
    app.config['NMF_ARGS'] = {
    }
    # set the test data dir correctly
    # find current directory, go up
    if 'TEST_DATA_DIR' in os.environ:
        app.config['TEST_DATA_DIR'] = os.environ['TEST_DATA_DIR']
    else:
        app.config['TEST_DATA_DIR'] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'test_data')
    if 'USER_DATA_DIR' in os.environ:
        app.config['USER_DATA_DIR'] = os.environ['USER_DATA_DIR']
    else:
        app.config['USER_DATA_DIR'] = '/tmp/uncurl/'
    if 'SECONDARY_USER_DATA_DIR' in os.environ:
        app.config['SECONDARY_USER_DATA_DIR'] = os.environ['SECONDARY_USER_DATA_DIR']
    else:
        app.config['SECONDARY_USER_DATA_DIR'] = '/tmp/uncurl_2/'
    app.config['BULK_DATA_DIR'] = 'bulk_data/'

    app.config['CACHE_TYPE'] = 'redis'

    if 'SHOW_ALL_RESULTS' in os.environ:
        show_all = os.environ['SHOW_ALL_RESULTS']
        if !show_all or show_all.lower() == 'false' or show_all == '0':
            app.config['SHOW_ALL_RESULTS'] = False
        else:
            app.config['SHOW_ALL_RESULTS'] = True
    else:
        app.config['SHOW_ALL_RESULTS'] = False

    # register blueprints
    app.register_blueprint(interaction_views.interaction_views)
    app.register_blueprint(views.views)
    app.register_blueprint(flask_router.flask_router)
    app.register_blueprint(db_query.db_query)
    @app.route('/')
    def index():
        return render_template('index.html')
    # handle errors
    from werkzeug.exceptions import HTTPException
    import traceback
    @app.errorhandler(Exception)
    def handle_exception(e):
        # pass through HTTP errors
        print(e)
        text = traceback.format_exc()
        if isinstance(e, HTTPException):
            return e
        # now you're handling non-HTTP exceptions only
        return render_template("error.html", msg=str(e) + '\n\n' + text), 500
    return app

def create_app_split_seq(data_dir='./', config_filename=None):
    # TODO: create an app starting at the data preprocessing view, where
    # data_dir is the DGE folder output from split-seq.
    from .split_seq_input import process_split_seq
    from flask import redirect, url_for, Blueprint
    if not os.path.exists(os.path.join(data_dir, 'preprocess.json')):
        process_split_seq(data_dir)
    app = Flask(__name__)
    Bootstrap(app)
    if 'MAX_CONTENT_LENGTH' in os.environ:
        app.config['MAX_CONTENT_LENGTH'] = int(os.environ['MAX_CONTENT_LENGTH'])
    else:
        app.config['MAX_CONTENT_LENGTH'] = 1000 * 1024 * 1024
    # default args to pass to uncurl.run_state_estimation
    app.config['UNCURL_ARGS'] = {
            'threads': 2,
            'max_iters': 20,
            'inner_max_iters': 50
    }
    app.config['NMF_ARGS'] = {
    }
    if 'TEST_DATA_DIR' in os.environ:
        app.config['TEST_DATA_DIR'] = os.environ['TEST_DATA_DIR']
    else:
        app.config['TEST_DATA_DIR'] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'test_data')
    # go up 1 from path
    app.config['USER_DATA_DIR'] = os.path.dirname(data_dir)
    app.config['CACHE_TYPE'] = 'simple'
    cache.config = {'CACHE_TYPE': 'simple'}
    app.config['SHOW_ALL_RESULTS'] = True
    app.register_blueprint(interaction_views.interaction_views)
    app.register_blueprint(views.views)
    app.register_blueprint(flask_router.flask_router)
    app.register_blueprint(db_query.db_query)
    # redirect index to data page
    user_id = os.path.basename(data_dir)
    print('user_id: ', user_id)
    @app.route('/index')
    def index():
        return render_template('index.html')

    bp2 = Blueprint('new_index', __name__)
    @bp2.route('/')
    def index_redirect():
        return redirect(url_for('views.state_estimation_result', user_id=user_id))
    app.register_blueprint(bp2)
    return app
