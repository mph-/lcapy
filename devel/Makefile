all: pic1.pdf pic2.pdf pic3.pdf pic4.pdf VRL1.pdf VacRL1.pdf pic5.pdf rnet1.pdf pic6.pdf rnet2.pdf rnet3.pdf rnet4.pdf

%.pytex: %.py
	cd $(dir $@); ipython $(notdir $<) > $(notdir $@) 

%.pytex: %.sch
	cd $(dir $@); ./schtex.py $(notdir $<) > $(notdir $@) 

%.pdf: %.tex %.pytex
	pdflatex $*

clean:
	rm *.pdf *.pytex

