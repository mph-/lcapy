install:
	python3 setup.py install

package:
	python3 setup.py sdist bdist_wheel

upload-test: package
	python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

upload: package
	python3 -m twine upload dist/*

test: lcapy/*.py
	nosetests3 --pdb

cover: lcapy/*.py
	nosetests3 --pdb --with-coverage --cover-package=lcapy --cover-html

doc-install: doc
	scp -r doc/_build/html/* lcapy.elec.canterbury.ac.nz:/var/www/lcapy/

release: doc push
	cd /tmp; rm -rf lcapy; git clone git@github.com:mph-/lcapy.git; cd lcapy; make test; make upload

release-test: doc push
	cd /tmp; rm -rf lcapy; git clone git@github.com:mph-/lcapy.git; cd lcapy; make test

style-check:
	-flake8 lcapy

push:
	git push
	git push --tags

.PHONY: doc
doc:
	cd doc; make html

.PHONY: clean
clean:
	-rm -rf build lcapy.egg-info dist
