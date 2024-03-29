# TODO: query by neural network classification?
# TODO: how to determine the gene set?

from keras.models import Sequential, load_model
from keras.layers import Dense

import numpy as np
from scipy import sparse


def data_preprocess_csc(data, normalize_counts=True, log_transform=True, **params):
    """
    run data preprocessing
    """
    from uncurl import preprocessing
    if normalize_counts:
        data = preprocessing.cell_normalize(data)
    if log_transform:
        data = preprocessing.log1p(data)
    return data

def sparse_batch_iterator(data, labels, batch_size=50, shuffle=True):
    """
    batch iterator for sparse matrices. Converts data to dense matrix.
    This returns an iterator over data, labels.
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
    This returns an iterator over the data.
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

    def __init__(self, gene_list, num_classes, class_names=None, model=None, layers=None):
        # TODO: layer parameters
        self.genes = gene_list
        self.class_names = class_names
        if model is not None:
            self.model = model
            self.model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
        else:
            self.model = Sequential()
            if layers is None:
                layers = [500, 200]
            for i, layer in enumerate(layers):
                if i == 0:
                    self.model.add(Dense(layer, activation='tanh', input_dim=len(gene_list)))
                else:
                    self.model.add(Dense(layer, activation='tanh'))
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
        if sparse.issparse(data):
            iterator = sparse_batch_iterator(data, labels, batch_size=batch_size)
            if sparse.issparse(validation_data):
                validation_iterator = sparse_batch_iterator(validation_data, validation_labels, batch_size=batch_size)
                self.model.fit_generator(iterator, steps_per_epoch=int(np.ceil(data.shape[0]/batch_size)), epochs=n_epochs,
                        validation_data=validation_iterator)
            else:
                self.model.fit_generator(iterator, steps_per_epoch=int(np.ceil(data.shape[0]/batch_size)), epochs=n_epochs)
        else:
            self.model.fit(data, labels, epochs=n_epochs)

    def predict(self, data, genes):
        """
        Args:
            data (array): dense array of shape (cells, genes)
            genes: array or list of strings

        Returns:
            model prediction for all cells - array of shape (cells, classes) of real values
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

    def results_to_labels(self, results):
        """
        Converts results array to list of cell type labels
        """
        labels = results.argmax(1)
        cell_labels = np.array([self.class_names[i] for i in labels])
        return cell_labels

    def save(self, model_filename, genes_filename=None):
        """
        save to filename - model_filename should be an h5 file.
        """
        self.model.save(model_filename)
        if genes_filename is None:
            genes_filename = model_filename.split('.')[0] + '_genes.txt'
        np.savetxt(genes_filename, self.genes, fmt='%s')
        class_filename = model_filename.split('.')[0] + '_classes.txt'
        np.savetxt(class_filename, self.class_names, fmt='%s')

    @classmethod
    def load_from_file(cls, model_filename, genes_filename=None):
        print(model_filename)
        model = load_model(model_filename)
        if genes_filename is None:
            genes_filename = model_filename.split('.')[0] + '_genes.txt'
        genes = np.loadtxt(genes_filename, dtype=str, delimiter='##')
        class_filename = model_filename.split('.')[0] + '_classes.txt'
        try:
            class_names = np.loadtxt(class_filename, dtype=str, delimiter='##')
        except:
            class_names = None
        classifier = Classifier(genes, 0, model=model, class_names=class_names)
        return classifier

import os
PATH = os.path.dirname(__file__)
DEFAULT_MODEL_PATH = os.path.join(PATH, 'models', 'tm_combined_model_200_200.h5')
LOADED_MODEL = None

def predict_using_default_classifier(data, genes, model_filename=None):
    """
    Returns cell type names using a default classifier model.

    Args:
        data (array): dense array of shape (cells, genes)
        genes: array or list of strings

    Returns:
        cell_names, cell_probs, cell_name_indices
    """
    model_filename = DEFAULT_MODEL_PATH if model_filename is None else model_filename
    global LOADED_MODEL
    LOADED_MODEL = Classifier.load_from_file(model_filename)
    results = LOADED_MODEL.predict(data, genes)
    cell_names = LOADED_MODEL.results_to_labels(results)
    return cell_names, results, LOADED_MODEL.class_names








class ClassifierAutoencoder(object):

    def __init__(self, gene_list, num_classes, model=None, layers=None, decoder_layers=None):
        self.genes = gene_list
        if model is not None:
            self.model = model
            self.model.compile(optimizer='adam',
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
        else:
            # TODO: use functional api
            self.model = Sequential()
            if layers is None:
                layers = [500, 200]
            if decoder_layers is None:
                decoder_layers = [200, 500]
            for i, layer in enumerate(layers):
                if i == 0:
                    self.model.add(Dense(layer, activation='tanh', input_dim=len(gene_list)))
                else:
                    self.model.add(Dense(layer, activation='tanh'))
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
        print('genes_map:', genes_map)
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
