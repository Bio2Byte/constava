from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="constava",
    version="0.1.0",
    author="Wim Vranken",
    author_email="wim.vranken@vub.be",
    description="This software is used to calculate conformational states probability & conformational state "
                "variability from a protein structure ensemble.",
    license="OSI Approved :: GNU General Public License v3 (GPLv3)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    maintainer="Jose Gavalda-Garcia, David Bickel, Adrian Diaz, Wim Vranken",
    maintainer_email="jose.gavalda.garcia@vub.be, david.bickel@vub.be, adrian.diaz@vub.be, wim.vranken@vub.be",
    url="https://bitbucket.org/bio2byte/constava/",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Development Status :: 5 - Production/Stable"
    ],
    python_requires=">=3.6, <3.11",
    install_requires=[
        "MDAnalysis",
        "numpy",
        "pandas",
        "scikit-learn",
    ],
    entry_points={
        "console_scripts": [
            "constava = constava.__main__:main",
        ],
    },
)
