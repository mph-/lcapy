all: pic1.pdf pic2.pdf pic3.pdf pic4.pdf VRL1.pdf VacRL1.pdf pic5.pdf rnet1.pdf

%.pytex: %.py
	cd $(dir $@); ipython $(notdir $<) > $(notdir $@) 

%.pdf: %.tex %.pytex
	pdflatex $*

clean:
	rm *.pdf *.pytex

