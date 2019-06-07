import os
from flask import Flask, render_template
from flask_bootstrap import Bootstrap

from . import interaction_views, views, flask_router, db_query

from .cache import cache

def create_app(config_filename=None):
    app = Flask(__name__)
    Bootstrap(app)

    # maximum file length is 250MB
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
    app.config['TEST_DATA_DIR'] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'test_data')
    app.config['USER_DATA_DIR'] = '/tmp/uncurl/'
    app.config['BULK_DATA_DIR'] = 'bulk_data/'

    app.config['CACHE_TYPE'] = 'redis'

    app.config['SHOW_ALL_RESULTS'] = False

    # register blueprints
    app.register_blueprint(interaction_views.interaction_views)
    app.register_blueprint(views.views)
    app.register_blueprint(flask_router.flask_router)
    app.register_blueprint(db_query.db_query)
    @app.route('/')
    def index():
        return render_template('index.html')
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
    app.config['MAX_CONTENT_LENGTH'] = 1000 * 1024 * 1024
    # default args to pass to uncurl.run_state_estimation
    app.config['UNCURL_ARGS'] = {
            'threads': 2,
            'max_iters': 20,
            'inner_max_iters': 50
    }
    app.config['NMF_ARGS'] = {
    }
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
