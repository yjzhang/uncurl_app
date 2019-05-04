import json
import os
import time
import unittest

from flask_testing import LiveServerTestCase
from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.by import By

from uncurl_app import create_app
from uncurl_app import cache


class UncurlFrontendTest(LiveServerTestCase):

    def create_app(self):
        app = create_app()
        app.config['DEPLOY'] = False
        app.config['LIVESERVER_PORT'] = 0
        cache.config = {'CACHE_TYPE': 'simple'}
        cache.init_app(app)
        return app

    def setUp(self):
        # webdriver download: https://github.com/mozilla/geckodriver/releases
        self.driver = webdriver.Firefox()
        self.driver.get(self.get_server_url())

    def test_submit(self):
        """
        Test submitting a new dataset to uncurl, as well as interacting
        with the dataset.
        """
        self.driver.find_element_by_id('state-estimation-link').click()
        time.sleep(1)
        self.driver.find_element_by_id('fileinput').send_keys(
                os.path.join(os.getcwd(), 'uncurl_app/static/data_400_cells.mtx.gz'))
        self.driver.find_element_by_id('genenames').send_keys(
                os.path.join(os.getcwd(), 'uncurl_app/static/gene_names_400.tsv'))
        self.driver.find_element_by_id('username').send_keys('testing')
        self.driver.find_element_by_id('submit').click()
        time.sleep(2)
        self.assertTrue('testing' in self.driver.current_url)
        self.driver.get(self.driver.current_url)
        num_cells = self.driver.find_element_by_id('num-cells').text
        self.assertTrue('400' in num_cells)
        num_genes = self.driver.find_element_by_id('num-genes').text
        self.assertTrue('19848' in num_genes)
        self.driver.find_element_by_id('vismethod')
        self.driver.find_element_by_id('baseline-vismethod')
        gene_frac = self.driver.find_element_by_id('genes-frac')
        gene_frac.clear()
        gene_frac.send_keys('0.5')
        cell_frac = self.driver.find_element_by_id('cell-frac')
        cell_frac.clear()
        cell_frac.send_keys('0.5')
        k = self.driver.find_element_by_id('k')
        k.clear()
        k.send_keys('8')
        time.sleep(1)
        self.driver.find_element_by_id('submit').click()
        time.sleep(60)
        # initial processing be completed by now
        print(self.driver.current_url)
        # TODO: refresh page
        self.driver.refresh()
        self.assertTrue('view' in self.driver.current_url)
        # testing scatterplot options
        self.driver.find_element_by_id('visualization')
        self.driver.find_element_by_id('means-scatter-plot')
        self.driver.find_element_by_id('Cells').click()
        time.sleep(1)
        self.driver.find_element_by_id('Baseline').click()
        time.sleep(1)
        # test all barplot options
        select = Select(self.driver.find_element_by_id('top-or-bulk'))
        select.select_by_value('pval')
        time.sleep(0.5)
        select.select_by_value('top_1_vs_rest')
        time.sleep(0.5)
        select.select_by_value('pval_1_vs_rest')
        time.sleep(0.5)
        select.select_by_value('top_pairwise')
        cluster1_select = Select(self.driver.find_element_by_id('barplot_cluster_select_1'))
        cluster1_select.select_by_value('4')
        cluster2_select = Select(self.driver.find_element_by_id('barplot_cluster_select_2'))
        cluster2_select.select_by_value('1')
        time.sleep(1.0)
        # test cell color options
        select = Select(self.driver.find_element_by_id('cell-color'))
        select.select_by_value('entropy')
        time.sleep(1)
        select.select_by_value('gene')
        # test coloring by gene
        self.driver.find_element_by_id('gene_name').send_keys('CD3D')
        self.driver.find_element_by_id('gene_name_submit').click()
        time.sleep(1)
        self.driver.find_element_by_id('gene_name').clear()
        self.driver.find_element_by_id('gene_name').send_keys('CD33,CD34')
        self.driver.find_element_by_id('gene_name_submit').click()
        time.sleep(1)
        # upload custom color map
        select.select_by_value('new')
        self.driver.find_element_by_id('color_track_file').send_keys(
                os.path.join(os.getcwd(), 'uncurl_app/static/labels_400_cells.txt'))
        select = Select(self.driver.find_element_by_id('color_track_type'))
        select.select_by_value('discrete')
        time.sleep(1)
        self.driver.find_element_by_id('submit_color_track').click()
        time.sleep(2)
        print(self.driver.current_url)
        print('color track upload done')
        time.sleep(2)
        # check that custom color track has been uploaded
        select = Select(self.driver.find_element_by_id('cell-color'))
        select.select_by_value('labels_400_cells')
        # test pairwise custom cluster
        select = Select(self.driver.find_element_by_id('top-or-bulk'))
        select.select_by_value('top_pairwise')
        cluster1_select = Select(self.driver.find_element_by_id('barplot_cluster_select_1'))
        cluster1_select.select_by_value('2')
        cluster2_select = Select(self.driver.find_element_by_id('barplot_cluster_select_2'))
        cluster2_select.select_by_value('0')
        time.sleep(1.0)
        # test split cluster
        self.driver.execute_script('window.all_selected_clusters = [0];')
        self.driver.find_element_by_id('recluster_toggle').click()
        self.driver.find_element_by_id('split').click()
        alert = self.driver.switch_to_alert()
        alert.accept()
        time.sleep(45)
        select = Select(self.driver.find_element_by_id('top-or-bulk'))
        select.select_by_value('top_1_vs_rest')
        time.sleep(1)

    def test_submit_2(self):
        # submit using Log-Normal
        self.driver.find_element_by_id('state-estimation-link').click()
        time.sleep(1)
        self.driver.find_element_by_id('fileinput').send_keys(
                os.path.join(os.getcwd(), 'uncurl_app/static/GSE60361_sub.mtx.gz'))
        self.driver.find_element_by_id('username').send_keys('testing-2')
        self.driver.find_element_by_id('submit').click()
        time.sleep(2)
        self.assertTrue('testing-2' in self.driver.current_url)
        self.driver.get(self.driver.current_url)
        num_cells = self.driver.find_element_by_id('num-cells').text
        self.assertTrue('753' in num_cells)
        num_genes = self.driver.find_element_by_id('num-genes').text
        self.assertTrue('3990' in num_genes)
        self.driver.find_element_by_id('vismethod')
        self.driver.find_element_by_id('baseline-vismethod')
        gene_frac = self.driver.find_element_by_id('genes-frac')
        gene_frac.clear()
        gene_frac.send_keys('0.5')
        cell_frac = self.driver.find_element_by_id('cell-frac')
        cell_frac.clear()
        cell_frac.send_keys('0.5')
        k = self.driver.find_element_by_id('k')
        k.clear()
        k.send_keys('7')
        select = Select(self.driver.find_element_by_id('disttype'))
        select.select_by_visible_text('Log-Normal')
        time.sleep(1)
        self.driver.find_element_by_id('submit').click()
        time.sleep(60)
        # should be completed by now
        print(self.driver.current_url)
        self.driver.refresh()
        self.assertTrue('view' in self.driver.current_url)
        self.driver.find_element_by_id('visualization')
        self.driver.find_element_by_id('means-scatter-plot')
        self.driver.find_element_by_id('Cells').click()
        time.sleep(1)
        self.driver.find_element_by_id('Baseline').click()
        time.sleep(1)
        select = Select(self.driver.find_element_by_id('top-or-bulk'))
        select.select_by_value('pval')
        time.sleep(1)
        self.driver.execute_script('window.all_selected_clusters = [0];')
        self.driver.find_element_by_id('recluster_toggle').click()
        self.driver.find_element_by_id('split').click()
        alert = self.driver.switch_to_alert()
        alert.accept()
        time.sleep(45)

    def test_reanalyze(self):
        """
        Test using an existing dataset
        """
        # TODO: test some more features - change visualizations,...
        # test all scatterplot and barplot options
        self.driver.find_element_by_link_text('Example results').click()
        self.driver.find_element_by_link_text('10x_400_new').click()
        self.driver.find_element_by_id('Baseline').click()
        time.sleep(1)
        # test CellMarker
        select = Select(self.driver.find_element_by_id('database-select'))
        select.select_by_value('cellmarker')
        self.driver.find_element_by_id('cellmarker-submit').click()
        time.sleep(1)
        # test cellmesh
        select = Select(self.driver.find_element_by_id('database-select'))
        select.select_by_value('cellmesh')
        self.driver.find_element_by_id('cellmesh-submit').click()
        time.sleep(1)
        # test mw
        select = Select(self.driver.find_element_by_id('cell-color'))
        select.select_by_value('gene')
        self.driver.find_element_by_id('gene_name').send_keys('CD3D')
        select = Select(self.driver.find_element_by_id('use_mw_gene'))
        select.select_by_value('1')
        self.driver.find_element_by_id('gene_name_submit').click()
        time.sleep(1)
        # TODO: test some custom color map options
        # TODO - this test is not currently working
        self.driver.execute_script('window.current_selected_cells = [' + ','.join("'" + str(x) + "'" for x in range(50)) + '];')
        self.driver.find_element_by_id('reanalyze').click()
        self.driver.find_element_by_id('subset_cells').click()
        alert = self.driver.switch_to_alert()
        alert.accept()
        time.sleep(10)

    def tearDown(self):
        self.driver.quit()


if __name__ == '__main__':
    unittest.main()
