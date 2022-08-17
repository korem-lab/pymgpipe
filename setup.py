from setuptools import setup

setup(
    name='pymgpipe',
    version='0.0.1',
    install_requires=[
        'micom',
        'cobra',
        'cplex',
        'gurobipy',
        'optlang',
        'tqdm',
    ],
    packages=['pymgpipe','pymgpipe.mseFBA'],
    url='https://github.com/korem-lab/pymgpipe'
)
