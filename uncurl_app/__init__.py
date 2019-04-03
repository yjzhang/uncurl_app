from flask import Flask
from flask_bootstrap import Bootstrap


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
import os
app.config['TEST_DATA_DIR'] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                'test_data')
app.config['USER_DATA_DIR'] = '/tmp/uncurl/'
app.config['BULK_DATA_DIR'] = 'bulk_data/'

app.config['CACHE_TYPE'] = 'redis'

app.config['SHOW_ALL_RESULTS'] = False

from . import views, flask_router, interaction_views
from .cache import cache