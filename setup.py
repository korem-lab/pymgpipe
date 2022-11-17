from setuptools import setup, find_packages

setup(
    name='pymgpipe',
    version='0.0.1',
    install_requires=[
        'cobra',
        'optlang',
        'tqdm',
    ],
    url='https://github.com/korem-lab/pymgpipe',
    package_dir={'pymgpipe':'src/pymgpipe'},
    packages=['pymgpipe'],
    package_data = {
        'pymgpipe': ['resources/diets/*.txt'],
    }
)
