import os
# TODO: user should change this (default should be test_data)
os.environ['TEST_DATA_DIR'] = '/cse/web/research/uncurl/uncurl_test'

from uncurl_app import create_app
from uncurl_app import cache

app = create_app()
app.config['BOOTSTRAP_SERVE_LOCAL'] = True
app.config['DEPLOY'] = False
app.config['SHOW_ALL_RESULTS'] = True

cache.config = {'CACHE_TYPE': 'simple'}
cache.init_app(app)

if __name__=='__main__':
    app.run(debug=True)
