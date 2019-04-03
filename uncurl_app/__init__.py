import os
from flask import Flask
from flask_bootstrap import Bootstrap

from . import interaction_views, views, flask_router

def create_app(config_filename=None):
    app = Flask(__name__)
    Bootstrap(app)

    # maximum file length is 1000MB
    app.config['MAX_CONTENT_LENGTH'] = 1000 * 1024 * 1024
    # default args to pass to uncurl.run_state_estimation
    app.config['UNCURL_ARGS'] = {
            'threads': 2,
            'max_iters': 20,
            'inner_max_iters': 50
    }
    app.config['NMF_ARGS'] = {
    }
    # TODO: set the test data dir correctly
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
    return app

from .cache import cache
