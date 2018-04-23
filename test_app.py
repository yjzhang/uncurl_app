import json
import os
import unittest

from app import app
from app import cache


class UncurlAppTest(unittest.TestCase):

    def setUp(self):
        app.config['DEPLOY'] = False
        cache.config = {'CACHE_TYPE': 'simple'}
        self.app = app.test_client()
        cache.init_app(app)

    def test_pages(self):
        index = self.app.get('/')
        self.assertEqual(index.status, '200 OK')
        self.assertTrue('UNCURL' in index.data)
        se_results = self.app.get('/state_estimation')
        self.assertTrue('UNCURL' in se_results.data)
        self.assertTrue('State Estimation' in se_results.data)
        self.assertEqual(se_results.status, '200 OK')
        data_results = self.app.get('/data')
        self.assertEqual(data_results.status, '200 OK')
        print(index)
        print(se_results)

    def test_data_get(self):
        page = self.app.get('/user/test_10x_400_new/view')
        self.assertEqual(page.status, '200 OK')
        barplot = self.app.post('/user/test_10x_400_new/view/update_barplot',
                data={'top_or_bulk': 'top',
                      'input_value': 0,
                      'num_genes': 10})
        self.assertEqual(barplot.status, '200 OK')
        barplot_data = json.loads(barplot.data)
        self.assertTrue('data' in barplot_data)
        self.assertTrue('layout' in barplot_data)
        self.assertTrue(len(barplot_data['data'][0]['x']) == 10)
        self.assertTrue(len(barplot_data['data'][0]['y']) == 10)
        barplot = self.app.post('/user/test_10x_400_new/view/update_barplot',
                data={'top_or_bulk': 'top',
                      'input_value': 4,
                      'num_genes': 20})
        self.assertEqual(barplot.status, '200 OK')
        barplot_data = json.loads(barplot.data)
        self.assertTrue('data' in barplot_data)
        self.assertTrue('layout' in barplot_data)
        self.assertTrue(len(barplot_data['data'][0]['x']) == 20)
        self.assertTrue(len(barplot_data['data'][0]['y']) == 20)

    def test_get_scatterplot(self):
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Means',
                      'cell_color': 'cluster'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data)
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Cells',
                      'cell_color': 'cluster'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data)
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Baseline',
                      'cell_color': 'cluster'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data)
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)


if __name__ == '__main__':
    unittest.main()
