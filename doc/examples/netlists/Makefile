PYFILES = $(wildcard *.py)
SCHFILES = $(wildcard *.sch)
PNGFILES = $(PYFILES:.py=.png) $(SCHFILES:.sch=.png)

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --matplotlib=Agg $(notdir $<)

%.png: %.sch
	cd $(dir $@); schtex $(notdir $<) $(notdir $@)
