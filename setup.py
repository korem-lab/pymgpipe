from setuptools import setup, find_packages

setup(
    name='pymgpipe',
    version='0.0.1',
    install_requires=[
        'micom',
        'cobra',
        'optlang',
        'tqdm',
    ],
    url='https://github.com/korem-lab/pymgpipe',
    packages=find_packages(),
    package_data = {
        'pymgpipe': ['*.mps'],
    }
)
