PYFILES = $(wildcard *.py)
SCHFILES = $(wildcard *.sch)
PNGFILES = $(PYFILES:.py=.png) $(SCHFILES:.sch=.png)

all: $(PNGFILES)

%.png: %.py
	# pylab option prevents ipython hanging on show
	cd $(dir $@); ipython --pylab=gtk $(notdir $<)

%.png: %.sch
	cd $(dir $@); schtex.py $(notdir $<) $(notdir $@) 
