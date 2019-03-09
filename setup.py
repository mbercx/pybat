from setuptools import setup, find_packages

setup(
    name="pybat",
    version="pre-alpha",
    packages=find_packages(exclude=["docs"]),
    install_requires=[
        "pymatgen",
        "click",
        "pymongo",
        "fireworks",
        "custodian",
        "tabulate",
        "icet"
    ],
    entry_points='''
        [console_scripts]
        pybat=pybat.cli.cli:main
    '''
)