
from flask_caching import Cache

cache = Cache()

def clear_cache_user_id(user_id):
    """
    Deletes all keys that match the given user_id in redis.
    """
    import redis
    r = redis.StrictRedis(host='localhost', port=6379, db=0)
    # TODO: this doesn't work due to hashing
    for key in r.scan_iter(match='*{0}*'.format(user_id)):
        r.delete(key)
