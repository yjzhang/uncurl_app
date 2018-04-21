from app import app
from app import cache

app.config['DEPLOY'] = True

cache.config = {'CACHE_TYPE': 'redis',
                'CACHE_REDIS_HOST': '127.0.0.1',
                'CACHE_REDIS_PORT': 6379,
                'CACHE_KEY_PREFIX': 'uncurl'}
cache.init_app(app)

if __name__ == '__main__':
    app.run()
