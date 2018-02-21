from flask import Flask
from flask_bootstrap import Bootstrap


app = Flask(__name__)
Bootstrap(app)

# maximum file length is 20MB
app.config['MAX_CONTENT_LENGTH'] = 20 * 1024 * 1024
# number of threads to use in uncurl
app.config['UNCURL_THREADS'] = 2
app.config['UNCURL_ARGS'] = {
        'threads': 2,
        'max_iters': 20,
        'inner_max_iters': 50
}

from app import views, flask_router
