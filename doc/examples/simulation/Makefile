PNGFILES = sim1.png sim2.png

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --matplotlib=gtk3 $(notdir $<)

clean:
	-rm $(PNGFILES)
