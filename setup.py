from setuptools import setup

setup(
    name='pymgpipe',
    version_format='{tag}.dev{commitcount}+{gitsha}',
    setup_requires=['setuptools-git-version'],
    install_requires=[
        'micom',
        'cobra',
        'cplex',
        'gurobipy',
        'optlang'
    ],
)