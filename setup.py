from setuptools import setup, find_namespace_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()

setup(
    name="constava",
    version="1.2.0b1",
    author="Wim Vranken",
    author_email="wim.vranken@vub.be",
    description="This software is used to calculate conformational states probability & conformational state "
    "variability from a protein structure ensemble.",
    license="GPL-3.0-only",
    license_files=["LICENSE", "authors.md"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    maintainer="Adrián Díaz",
    maintainer_email="bio2byte@vub.be, adrian.diaz@vub.be",
    url="https://github.com/bio2byte/constava/",
    packages=find_namespace_packages(
        include=["constava*"],
    ),
    include_package_data=True,
    classifiers=[
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Development Status :: 5 - Production/Stable",
    ],
    python_requires=">=3.8,<3.15",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "constava = constava.__main__:main",
        ],
    },
    package_data={
        "constava": ["constava/data/*"],
    },
)
