from app import app
from app import cache

app.config['BOOTSTRAP_SERVE_LOCAL'] = True
app.config['DEPLOY'] = False
app.config['SHOW_ALL_RESULTS'] = True

cache.config = {'CACHE_TYPE': 'simple'}
cache.init_app(app)

if __name__=='__main__':
    app.run(debug=True)
