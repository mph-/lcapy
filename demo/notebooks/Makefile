IPYNB = $(wildcard *.ipynb)
IPYNB_UPDATE = $(IPYNB:.ipynb=.ipynb.update)
IPYNB_HTML = $(IPYNB:.ipynb=.html)
IPYNB_HTML_FILES = $(IPYNB:.ipynb=_files)
IPYNB_PDF = $(IPYNB:.ipynb=.pdf)
IPYNB_CHECK = $(IPYNB:.ipynb=.ipynb.check)

IPYTHON=ipython

%.html: %.ipynb
	$(IPYTHON) nbconvert --to html $<
	mv $@ $@.bak
	# Tweak file to allow latex math to be displayed correctly
	sed 's|c328740.ssl.cf1.rackcdn.com|cdn.mathjax.org|g' $@.bak > $@
	rm $@.bak

%.pdf: %.ipynb
	$(IPYTHON) nbconvert --to latex --post PDF $<

html: $(IPYNB_HTML)

pdf: $(IPYNB_PDF)

update: $(IPYNB_UPDATE)

%.ipynb.update: %.ipynb
	cp $^ $@
	runipy --overwrite $^

check: $(IPYNB_CHECK)
	echo $^

%.ipynb.check: %.ipynb
	runipy $^ $@

clean:
	-rm $(IPYNB_PDF)
	-rm $(IPYNB_HTML)
	-rm -r $(IPYNB_HTML_FILES)

