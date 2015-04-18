IPYNB = $(wildcard *.ipynb)
IPYNB_UPDATE = $(IPYNB:.ipynb=.ipynb.update)
IPYNB_TEX = $(IPYNB:.ipynb=.tex)
IPYNB_HTML = $(IPYNB:.ipynb=.html)
IPYNB_HTML_FILES = $(IPYNB:.ipynb=_files)
IPYNB_PDF = $(IPYNB:.ipynb=.pdf)
IPYNB_RST = $(IPYNB:.ipynb=.rst)
IPYNB_CHECK = $(IPYNB:.ipynb=.ipynb.check)

IPYTHON=ipython

all: pdf

pdf: $(IPYNB_PDF)

#%.pdf: %.ipynb
#	$(IPYTHON) nbconvert --to latex --post PDF $<

%.pdf: %.tex
	pdflatex $*

rst: $(IPYNB_RST)

%.rst: %.ipynb
	$(IPYTHON) nbconvert --to rst $<

html: $(IPYNB_HTML)

%.html: %.ipynb
	$(IPYTHON) nbconvert --to html $<
	mv $@ $@.bak
	# Tweak file to allow latex math to be displayed correctly
	sed 's|c328740.ssl.cf1.rackcdn.com|cdn.mathjax.org|g' $@.bak > $@
	rm $@.bak

tex: $(IPYNB_TEX)

%.tex: %.ipynb
	$(IPYTHON) nbconvert --to latex  $<
	mv $@ $@.bak
	sed 's|documentclass|documentclass[fleqn]|g' $@.bak > $@
	rm $@.bak

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
	-rm $(IPYNB_TEX)
	-rm $(IPYNB_RST)
	-rm -r $(IPYNB_HTML_FILES)

