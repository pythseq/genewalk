import sys
from setuptools import setup


def main():
    install_list = ['numpy', 'pandas', 'networkx>=2.1', 'gensim', 'goatools']

    setup(name='genewalk',
          version='0.0.1',
          description='Determine gene function based on network embeddings.',
          long_description=(''),
          author='Robert Ietswaart',
          author_email='robert_ietswaart@hms.harvard.edu',
          url='https://github.com/churchmanlab/GeneWalk',
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
          keywords=['gene function', 'network', 'embedding'],
          packages=['genewalk'],
          install_requires=install_list,
          extras_require={'indra': ['indra']},
          tests_require=['nose'],
          include_package_data=True,
        )


if __name__ == '__main__':
    main()
