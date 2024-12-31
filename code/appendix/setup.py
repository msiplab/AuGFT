from setuptools import setup, find_packages

setup(
    name='py_augft',
    version='0.1',
    packages=find_packages(),
    python_requires='>=3.10',    
    install_requires=[
        'numpy',
        'networkx',
        ],
    author='Shogo MURAMATSU',
    author_email='shogo@eng.niigata-u.ac.jp',
    description='Supplemental materials for article "Realization of DiGraph Filters Via Augmented GFT"',
    url='https://github.com/msiplab/AuGFT'
)
