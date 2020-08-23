PNGFILES = sim1.png sim2.png VRC1.png VRL1.png

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --matplotlib=gtk3 $(notdir $<)

%.png: %.sch
	cd $(dir $@); schtex $(notdir $<) $(notdir $@) 

clean:
	-rm $(PNGFILES)
