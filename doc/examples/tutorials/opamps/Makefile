PYFILES = opamp-piezo-amplifier2-asd2.py opamp-piezo-amplifier2-frequency-response1.py vnoise1.py
SCHFILES = $(wildcard *.sch)
PNGFILES = $(PYFILES:.py=.png) $(SCHFILES:.sch=.png)

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --matplotlib=gtk3 $(notdir $<)

%.png: %.sch
	cd $(dir $@); schtex $(notdir $<) $(notdir $@) 
