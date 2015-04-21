SCH = $(wildcard *.sch)
PDF = $(SCH:.sch=.pdf)
PNG = $(SCH:.sch=.png)

all: png

png: $(PNG)

pdf: $(PDF)

%.pytex: %.sch
	cd $(dir $@); schtex.py $(notdir $<) $(notdir $@) 

%.pdf: %.sch
	cd $(dir $@); schtex.py $(notdir $<) $(notdir $@) 

%.png: %.sch
	cd $(dir $@); schtex.py $(notdir $<) $(notdir $@) 

clean:
	-rm *.pdf
	-rm *.pytex
	-rm *.png

.PRECIOUS: %.pytex
