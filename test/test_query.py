import unittest

import numpy as np
from uncurl_analysis import dense_matrix_h5
import mouse_cell_query

class QueryTest(unittest.TestCase):

    def setUp(self):
        self.data = dense_matrix_h5.H5Dict('mouse_cell_query/data/cell_type_means.h5')
        self.gene_names = np.loadtxt('mouse_cell_query/data/gene_names.txt', dtype=str)

    def test_query(self):
        results = mouse_cell_query.search_db(self.data['fibroblast'], self.gene_names)
        print(results)
        self.assertTrue(results[0][0] == 'fibroblast')

if __name__ == '__main__':
    unittest.main()
