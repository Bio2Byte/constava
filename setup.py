from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="constava",
    version="0.0.1b",
    author="Bio2Byte",
    author_email="bio2byte@ibsquare.be",
    description="Command line tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/bio2byte/constava/",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
    ],
    python_requires=">=3.6",
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
