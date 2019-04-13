install:
	python3 setup.py install

package:
	python3 setup.py sdist bdist_wheel

test: lcapy/*.py
	nosetests3 --pdb

cover: lcapy/*.py
	nosetests3 --pdb --with-coverage --cover-package=lcapy --cover-html
