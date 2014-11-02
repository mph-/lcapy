all: pic1.pdf pic2.pdf pic3.pdf

%.pytex: %.py
	cd $(dir $@); ipython $(notdir $<) > $(notdir $@) 

%.pdf: %.tex %.pytex
	pdflatex $*

clean:
	rm *.pdf *.pytex

