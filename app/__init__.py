from flask import Flask
from flask_bootstrap import Bootstrap


app = Flask(__name__)
Bootstrap(app)

# maximum file length is 20MB
app.config['MAX_CONTENT_LENGTH'] = 20 * 1024 * 1024
# default args to pass to uncurl.run_state_estimation
app.config['UNCURL_ARGS'] = {
        'threads': 2,
        'max_iters': 20,
        'inner_max_iters': 50
}
app.config['TEST_DATA_DIR'] = 'test_data/'
app.config['USER_DATA_DIR'] = '/tmp/uncurl/'
app.config['BULK_DATA_DIR'] = 'bulk_data/'

from app import views, flask_router
