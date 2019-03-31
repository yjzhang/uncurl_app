import json
import os
import unittest

from uncurl_app import generate_analysis
from uncurl_app import app
from uncurl_app import cache


class UncurlAppTest(unittest.TestCase):

    def setUp(self):
        app.config['DEPLOY'] = False
        cache.config = {'CACHE_TYPE': 'simple'}
        self.app = app.test_client()
        cache.init_app(app)

    def test_pages(self):
        index = self.app.get('/')
        self.assertEqual(index.status, '200 OK')
        index_data = index.data.decode('utf-8')
        self.assertTrue('UNCURL' in index_data)
        se_results = self.app.get('/state_estimation')
        se_results_data = se_results.data.decode('utf-8')
        self.assertTrue('UNCURL' in se_results_data)
        self.assertTrue('State Estimation' in se_results_data)
        self.assertEqual(se_results.status, '200 OK')
        data_results = self.app.get('/data')
        self.assertEqual(data_results.status, '200 OK')
        print(index)
        print(se_results)

    def test_data_get(self):
        """
        Test getting barplot data
        """
        page = self.app.get('/user/test_10x_400_new/view')
        self.assertEqual(page.status, '200 OK')
        barplot = self.app.post('/user/test_10x_400_new/view/update_barplot',
                data={'top_or_bulk': 'top',
                      'input_value': 0,
                      'num_genes': 10})
        self.assertEqual(barplot.status, '200 OK')
        barplot_data = json.loads(barplot.data.decode('utf-8'))
        self.assertTrue('data' in barplot_data)
        self.assertTrue('layout' in barplot_data)
        self.assertTrue(len(barplot_data['data'][0]['x']) == 10)
        self.assertTrue(len(barplot_data['data'][0]['y']) == 10)
        barplot = self.app.post('/user/test_10x_400_new/view/update_barplot',
                data={'top_or_bulk': 'top',
                      'input_value': 4,
                      'num_genes': 20})
        self.assertEqual(barplot.status, '200 OK')
        barplot_data = json.loads(barplot.data.decode('utf-8'))
        self.assertTrue('data' in barplot_data)
        self.assertTrue('layout' in barplot_data)
        self.assertTrue(len(barplot_data['data'][0]['x']) == 20)
        self.assertTrue(len(barplot_data['data'][0]['y']) == 20)
        barplot = self.app.post('/user/test_10x_400_new/view/update_barplot',
                data={'top_or_bulk': 'top_1_vs_rest',
                      'input_value': 4,
                      'num_genes': 15})
        self.assertEqual(barplot.status, '200 OK')
        barplot_data = json.loads(barplot.data.decode('utf-8'))
        self.assertTrue('data' in barplot_data)
        self.assertTrue('layout' in barplot_data)
        self.assertTrue(len(barplot_data['data'][0]['x']) == 15)
        self.assertTrue(len(barplot_data['data'][0]['y']) == 15)



    def test_get_scatterplot(self):
        """
        Test getting scatterplot data
        """
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Means',
                      'cell_color': 'cluster'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data.decode('utf-8'))
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)
        self.assertTrue(scatterplot_data['data'][0]['name'] == 'cluster 0')
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Cells',
                      'cell_color': 'cluster'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data.decode('utf-8'))
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)
        self.assertTrue(scatterplot_data['data'][1]['name'] == 'cluster 1')
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Baseline',
                      'cell_color': 'cluster'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data.decode('utf-8'))
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Cells',
                      'cell_color': 'gene',
                      'gene_name': 'CD8B'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data.decode('utf-8'))
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)
        # gene-gene
        scatterplot = self.app.post('/user/test_10x_400_new/view/update_scatterplot',
                data={'scatter_type': 'Gene-gene',
                      'cell_color': 'cluster',
                      'gene_name_1': 'CD8B',
                      'gene_name_2': 'REST'})
        self.assertEqual(scatterplot.status, '200 OK')
        scatterplot_data = json.loads(scatterplot.data.decode('utf-8'))
        self.assertTrue('data' in scatterplot_data)
        self.assertTrue('layout' in scatterplot_data)
        self.assertTrue(len(scatterplot_data['data']) == 8)


    def test_helper_functions(self):
        current_task, time_remaining = generate_analysis.get_progress('test_data/10x_400_new')
        print(current_task)
        print(time_remaining)

    # TODO: test submitting a file?


if __name__ == '__main__':
    unittest.main()
