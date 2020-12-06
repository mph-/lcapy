.PHONY: install
install:
	python3 setup.py install

.PHONY: package
package:
	python3 setup.py sdist bdist_wheel

.PHONY: upload-test
upload-test: package
	python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

.PHONY: upload
upload: package
	python3 -m twine upload dist/*

.PHONY: test
test: lcapy/*.py
	nosetests3 --pdb

.PHONY: cover
cover: lcapy/*.py
	nosetests3 --pdb --with-coverage --cover-package=lcapy --cover-html

.PHONY: doc-install
doc-install: doc
	scp -r doc/_build/html/* lcapy.elec.canterbury.ac.nz:/var/www/lcapy/

.PHONY: doc
release: doc push
	cd /tmp; rm -rf lcapy; git clone git@github.com:mph-/lcapy.git; cd lcapy; make test; make upload

.PHONY: release-test
release-test: doc push
	cd /tmp; rm -rf lcapy; git clone git@github.com:mph-/lcapy.git; cd lcapy; make test

.PHONY: style-check
style-check:
	-flake8 lcapy

.PHONY: flake8
flake8:
	-flake8 lcapy

.PHONY: push
push:
	git push
	git push --tags

.PHONY: doc
doc:
	cd doc; make html

.PHONY: clean
clean:
	-rm -rf build lcapy.egg-info dist
