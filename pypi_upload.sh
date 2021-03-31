#!/bin/bash

rm -r dist/
python setup.py sdist
python -m twine upload dist/*