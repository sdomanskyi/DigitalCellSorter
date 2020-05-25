echo Removing any previous distributions
del /q "dist\*.*"

python setup.py sdist bdist_wheel

python -m twine check dist/*