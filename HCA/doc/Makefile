
all:    hcaman.pdf hca.1 clean

hca.1:  hcaman.tex
	latex2man hcaman.tex hca.1

hcaman.pdf:  hcaman.tex
	pdflatex hcaman.tex

clean:
	rm -f hcaman.tex~ hcaman.log hcaman.aux 
