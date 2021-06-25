PYFILES = $(wildcard *.py)
PNGFILES = $(PYFILES:.py=.png)

all: $(PNGFILES)

%.png: %.py function.mplstyle
	ipython --matplotlib=gtk3 $<

clean:
	-rm $(PNGFILES)
