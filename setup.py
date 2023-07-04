from setuptools import setup, find_packages

setup(
    author="Yoli Meydan, Federico Baldini, Tal Korem",
    author_email="ym2877@cumc.columbia.edu, fb2557@cumc.columbia.edu, tk2829@cumc.columbia.edu",
    name="pymgpipe",
    description="Community level microbiome metabolic modeling in Python",
    version="v0.16.1",
    classifiers=[
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=[
        "cobra",
        "optlang",
        "scikit_bio",
        "optlang",
        "tqdm",
        "seaborn",
        "gurobipy"
    ],
    url="https://github.com/korem-lab/pymgpipe",
    package_dir={"pymgpipe": "pymgpipe"},
    packages=find_packages(include=["pymgpipe"]),
    package_data={
        "pymgpipe": [
            "resources/diets/*.txt",
            "resources/models/*",
            "resources/problems/*",
            "resources/miniTaxa/*",
            "resources/.VERSION"
        ],
    },
    project_urls={
        "Issues": "https://github.com/korem-lab/pymgpipe/issues",
        "Source": "https://github.com/korem-lab/pymgpipe",
        "Readme": "https://github.com/korem-lab/pymgpipe/blob/main/README.md",
    },
    python_requires=">=3",
    license="Apache-2.0",
    license_files=["LICENSE"],
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
)

