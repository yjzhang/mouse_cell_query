# TODO: query by neural network classification?
# TODO: how to determine the gene set?

from keras.models import Sequential, load_model
from keras.layers import Dense

import numpy as np
from scipy import sparse

def sparse_batch_iterator(data, labels, batch_size=50, shuffle=True):
    """
    batch iterator for sparse matrices. Converts data to dense matrix.
    see: https://stackoverflow.com/questions/37609892/keras-sparse-matrix-issue
    """
    n_batches = int(np.ceil(data.shape[0]/batch_size))
    counter = 0
    shuffle_index = np.arange(data.shape[0])
    if shuffle:
        np.random.shuffle(shuffle_index)
        # shuffle data and labels
        data = data[shuffle_index, :]
        labels = labels[shuffle_index]
    while True:
        index_batch = shuffle_index[batch_size*counter:batch_size*(counter+1)]
        X_batch = data[index_batch,:].todense()
        y_batch = labels[index_batch]
        counter += 1
        yield (np.array(X_batch), y_batch)
        if (counter >= n_batches):
            if shuffle:
                np.random.shuffle(shuffle_index)
        counter=0

def sparse_batch_data_iterator(data, batch_size=50, shuffle=True):
    """
    batch iterator for sparse matrices. Converts data to dense matrix.
    see: https://stackoverflow.com/questions/37609892/keras-sparse-matrix-issue
    """
    n_batches = int(np.ceil(data.shape[0]/batch_size))
    counter = 0
    shuffle_index = np.arange(data.shape[0])
    # shuffle data and labels
    if shuffle:
        np.random.shuffle(shuffle_index)
        data = data[shuffle_index, :]
    while True:
        index_batch = shuffle_index[batch_size*counter:batch_size*(counter+1)]
        X_batch = data[index_batch,:].todense()
        counter += 1
        yield np.array(X_batch)
        if (counter >= n_batches):
            if shuffle:
                np.random.shuffle(shuffle_index)
            counter = 0


class Classifier(object):

    def __init__(self, gene_list, num_classes, model=None):
        # TODO: layer parameters
        self.genes = gene_list
        if model is not None:
            self.model = model
            self.model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
        else:
            self.model = Sequential()
            self.model.add(Dense(500, activation='tanh', input_dim=len(gene_list)))
            self.model.add(Dense(200, activation='tanh'))
            self.model.add(Dense(num_classes, activation='softmax'))
            self.model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])

    def fit(self, data, labels, n_epochs=20, batch_size=50, validation_data=None, validation_labels=None):
        """
        Args:
            data (array): dense or sparse array of shape (cells, genes).
            labels (array or list): cell labels - should be ints.
        """
        # TODO: validation set?
        if sparse.issparse(validation_data):
            validation_iterator = sparse_batch_iterator(validation_data, validation_labels, batch_size=batch_size)
        if sparse.issparse(data):
            iterator = sparse_batch_iterator(data, labels, batch_size=batch_size)
            self.model.fit_generator(iterator, steps_per_epoch=int(np.ceil(data.shape[0]/batch_size)), epochs=n_epochs)
        else:
            self.model.fit(data, labels, epochs=n_epochs)

    def predict(self, data, genes):
        """
        Args:
            data (array): dense array of shape (cells, genes)
            genes: array or list of strings

        Returns:
            model prediction for all cells
        """
        # map data gene set onto classifier gene set
        genes = list(map(lambda x: x.upper(), genes))
        self_genes = list(map(lambda x: x.upper(), self.genes))
        genes_map = np.zeros(len(self_genes), dtype=int)
        model_gene_index = {gene: i for i, gene in enumerate(self_genes)}
        genes_set = set(genes)
        zero_genes = [i for gene, i in model_gene_index.items() if gene not in genes_set]
        for i, gene in enumerate(genes):
            if gene in model_gene_index:
                genes_map[model_gene_index[gene]] = i
        print(genes_map)
        data_new = data[:, genes_map]
        data_new[:, zero_genes] = 0
        if sparse.issparse(data_new):
            batch_size = 50
            iterator = sparse_batch_data_iterator(data, batch_size=batch_size, shuffle=False)
            results = self.model.predict_generator(iterator, steps=int(np.ceil(data.shape[0]/batch_size)))
        else:
            results = self.model.predict(data)
        return results

    def save(self, model_filename, genes_filename=None):
        """
        save to filename - model_filename should be an h5 file.
        """
        self.model.save(model_filename)
        if genes_filename is None:
            genes_filename = model_filename.split('.')[0] + '_genes.txt'
        np.savetxt(genes_filename, self.genes, fmt='%s')

    @classmethod
    def load_from_file(cls, model_filename, genes_filename=None):
        print(model_filename)
        model = load_model(model_filename)
        if genes_filename is None:
            genes_filename = model_filename.split('.')[0] + '_genes.txt'
        genes = np.loadtxt(genes_filename, dtype=str, delimiter='##')
        classifier = Classifier(genes, 0, model=model)
        return classifier
