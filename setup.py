from setuptools import setup, find_packages

setup(
    name="pybat",
    version="0.1",
    packages=find_packages(exclude=["docs"]),
    install_requires=[
        "pymatgen",
        "numpy",
        "scipy",
        "click",
        "fireworks",
        "custodian"
    ],
    entry_points='''
        [console_scripts]
        pybat=pybat.cli.cli:main
    '''
)