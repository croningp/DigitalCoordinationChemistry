import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="carboxy",
    version="0.1.0",
    author="Daniel Kowalski",
    author_email="",
    description="Code for the Carboxylates Project",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    project_urls={
        "Bug Tracker": "",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=[
        "matplotlib>=3.4",
        "numpy>=1.21",
        "pandas>=1.2",
        "scikit-learn>=0.24",
        "scikit-optimize>=0.9",
        "summit>=0.8"
    ],
    python_requires=">=3.7",
)