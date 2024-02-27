from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='ic4popsyn',
    version='0.3.1',
    author='Nicola Giacobbo, Giuliano Iorio, Guglielmo Costa',
    author_email="giuliano.iorio.astro@gmail.com",
    description="A package to help build initial conditions for population-synthesis codes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/GiacobboNicola/IC4popsyn',
    project_urls={
        "Bug Tracker": "https://github.com/GiacobboNicola/IC4popsyn/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
    ],
    packages=['ic4popsyn',],
    install_requires=['numpy','pandas','progressbar','pytest','scipy'],
    python_requires=">=3.7"
)
