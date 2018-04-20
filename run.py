from app import app
from flask_caching import Cache

app.config['BOOTSTRAP_SERVE_LOCAL'] = True
app.config['DEPLOY'] = False

app.cache = Cache()
app.cache.config = {'CACHE_TYPE': 'simple'}
cache.init_app(app)

if __name__=='__main__':
    app.run(debug=True)
