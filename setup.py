from setuptools import setup, find_packages

install_requires = [
        'numpy',
        'scipy',
        'keras',
        'tensorflow',
        'uncurl_seq',
]

setup(
    name='mouse_cell_query',
    version='0.0.1',
    author='Yue Zhang',
    author_email='yjzhang@cs.washington.edu',
    url='https://github.com/yjzhang/mouse_cell_query',
    license='MIT',
    packages=find_packages("."),
    install_requires=install_requires,

    zip_safe=False,
    package_data={'mouse_cell_query': ['data/cell_type_means.h5',
        'data/gene_names.txt', 'data/cell_type_ontology_ids.pkl']},
    include_package_data=True,

    test_suite='nose.collector',
    tests_require=['nose', 'flaky'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
    ],
)
