PYFILES = $(wildcard *.py)
PNGFILES = $(PYFILES:.py=.png)

all: $(PNGFILES)

%.png: %.py
	# pylab option prevents ipython hanging on show
	cd $(dir $@); ipython --pylab=gtk $(notdir $<)
