import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="example-pkg-GNicola", # Replace with your own username
    version="0.0.0.1",
    author="Nicola Giacobbo",
    author_email="giacobbo.nicola@gmail.com",
    description="A package to help build initial conditions for popsynthesis codes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/GiacobboNicola/IC4popsyn",
    project_urls={
        "Bug Tracker": "https://github.com/GiacobboNicola/IC4popsyn/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.8",
)
