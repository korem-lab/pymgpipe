from setuptools import setup

setup(
    name='pymgpipe',
    setuptools_git_versioning={
        "enabled": True,
    },
    setup_requires=["setuptools-git-versioning"],
    install_requires=[
        'micom',
        'cobra',
        'cplex',
        'gurobipy',
        'optlang'
    ],
)