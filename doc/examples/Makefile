PYFILES = $(wildcard *.py)
PNGFILES = $(PYFILES:.py=.png)

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --pylab=gtk $(notdir $<)
