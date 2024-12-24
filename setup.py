import setuptools
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setuptools.setup(
    name="residPCA-test",
    author="Shaye Carver",
    version="1.0.4",
    author_email="scarver@g.harvard.edu",
    description="ResidPCA Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/carversh/residPCA",
    project_urls={
        "Bug Tracker": "https://github.com/carversh/residPCA/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    entry_points={
        'console_scripts': [
            'residPCA = residPCA:main',
        ],
    }  # No trailing comma here
)

