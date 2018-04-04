"""
Adds tips to a backbone phylogeny using taxonomy simulated with birth-death models
"""
from setuptools import find_packages, setup

dependencies = ['click', 'dendropy', 'numpy']

setup(
    name='tact',
    version='0.1.0',
    url='https://github.com/jonchang/tact',
    license='BSD',
    author='Jonathan Chang',
    author_email='jonathan.chang@ucla.edu',
    description='Adds tips to a backbone phylogeny using taxonomy simulated with birth-death models',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'tact_build_taxonomic_tree = tact.cli_taxonomy:main',
            'tact_add_taxa = tact.cli_add_taxa:main',
            'tact_check_results = tact.cli_check_trees:main',
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
        'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Operating System :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
