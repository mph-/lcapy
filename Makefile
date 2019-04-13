install:
	python3 setup.py install

package:
	python3 setup.py sdist bdist_wheel

upload-test: package
	python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

test: lcapy/*.py
	nosetests3 --pdb

cover: lcapy/*.py
	nosetests3 --pdb --with-coverage --cover-package=lcapy --cover-html

.PHONY: doc
doc:
	cd doc; make html
