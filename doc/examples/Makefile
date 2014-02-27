PYFILES = $(wildcard *.py)
PNGFILES = $(PYFILES:.py=.png)

all: $(PNGFILES)

%.png: %.py
	# Pylab option prevent ipyton hanging on show
	cd $(dir $@); ipython --pylab=gtk $(notdir $<)
