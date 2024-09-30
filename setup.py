from setuptools import setup, find_packages

setup(
    name="toegen",
    version="0.1.0",
    author="Orezbm",
    author_email="igemtau2024@gmail.com",
    description="A package for generating toehold switches with homology ranking for gene expression control",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Orezbm/toegen",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython",
        "pandas",
        "joblib",
        "fuzzysearch",
        "viennarna",
    ],
    python_requires='>=3.6',
    include_package_data=True,
    package_data={
        '': ['README.md'],
    },
)
