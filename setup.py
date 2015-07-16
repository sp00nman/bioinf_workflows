from setuptools import setup

def readme():
    with open('README.md' as f:
        return f.read()

setup(
    name='bioinf_workflows',
    version='0.2',
    description='Collection of bioinformatics workflows to analyze next-generation sequening data',
    long_description=readme(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords='bioinformatics ngs screen',
    url='https://github.com/sp00nman/bioinf_workflows',
    author='Fiorella Schischlik',
    author_email='fiorella.schischlik@gmail.com',
    license="GPL2",
    install_requires=[
        'pandas'
    ],
    scripts=[
        'bin/gt_screen_workflow.py',
        'bin/count_insertions.R',
        'bin/fisher-test.R',
        'bin/plot_screen.R'
    ],
    include_package_data=True
    ]
)
