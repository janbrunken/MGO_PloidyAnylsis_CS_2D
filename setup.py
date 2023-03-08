from setuptools import setup, find_packages

# Setting up
setup(
    name="PloidyAnalysis_2D",
    version="0.1",
    author="Jan Brunken",
    author_email="j.brunken@dkfz-heidelberg.de",
    description='Analyze ploidy from 2D image data.',
    packages=['PloidyAnalysis_2D'],
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'scikit-image',
        'natsort',
        'seaborn',
        'nbformat',
        'plotly',
        'tensorflow-cpu',
        'stardist',
        'cellpose',
        'tifffile',
        'ipython',
        'ipywidgets',
        'ipykernel'
    ]
)
