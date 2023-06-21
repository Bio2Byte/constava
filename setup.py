from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="constava",
    version="0.0.1",
    author="Bio2Byte",
    author_email="bio2byte@ibsquare.be",
    description="Command line tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/bio2byte/constava/",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
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
            "constava_template = template.__main__:main",
            "constava_template_2 = template_2.__main__:main",
            "constava = constava.__main__:main",

        ],
    },
)
