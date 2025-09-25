TMHOME=/Applications/TeXmacs.app/Contents/MacOS/TeXmacs

TARGET = dealii.pdf

all: $(TARGET)

%.pdf: %.tm
	$(TMHOME) -c $*.tm $*.pdf -q

dealii.pdf: dealii.tm
	$(TMHOME) -c dealii.tm dealii.pdf -q
	cp dealii.pdf ~/Dropbox/Teaching/slides/dealii.pdf

clean:
	rm -rf $(TARGET)
