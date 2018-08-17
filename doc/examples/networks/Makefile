PYFILES = $(wildcard *.py)
PNGFILES = $(PYFILES:.py=.png)

all: $(PNGFILES)

%.png: %.py
	cd $(dir $@); ipython --matplotlib=gtk3 $(notdir $<)

clean:
	-rm $(PNGFILES)
