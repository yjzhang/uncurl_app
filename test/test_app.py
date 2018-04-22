import os
import unittest

from app import app
from app import cache


class UncurlAppTest(unittest.TestCase):

    def setUp(self):
        app.config['DEPLOY'] = False
        cache.config = {'CACHE_TYPE': 'simple'}
        self.app = app.test_client
        cache.init_app(app)
        pass

if __name__ == '__main__':
    unittest.main()
