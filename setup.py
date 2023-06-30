from setuptools import setup, find_packages

with open('.VERSION') as f:
    version = f.readline().strip()[1:]

setup(
    name="pymgpipe",
    version="0.0.1",
    install_requires=[
        "cobra",
        "optlang",
        "tqdm",
    ],
    url="https://github.com/korem-lab/pymgpipe",
    package_dir={"pymgpipe": "pymgpipe"},
    packages=["pymgpipe"],
    package_data={
        "pymgpipe": [
            "resources/diets/*.txt",
            "resources/models/*",
            "resources/problems/*",
            "resources/miniTaxa/*",
        ],
    },
)
