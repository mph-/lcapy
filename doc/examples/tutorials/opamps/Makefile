PYFILES = opamp-piezo-amplifier2-asd2.py opamp-piezo-amplifier2-frequency-response1.py vnoise1.py inoise1.py opamp-open-loop1.py opamp-voltage-follower-C-load-open-loop1.py opamp-voltage-follower-RC-load-open-loop1.py opamp-voltage-follower-C-load-closed-loop-pole-plot1.py opamp-voltage-follower-RC-load-closed-loop-pole-plot1.py
SCHFILES = $(wildcard *.sch)
PNGFILES = $(PYFILES:.py=.png) $(SCHFILES:.sch=.png)

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --matplotlib=gtk3 $(notdir $<)

%.png: %.sch
	cd $(dir $@); schtex $(notdir $<) $(notdir $@)
