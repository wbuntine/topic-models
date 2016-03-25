
all:    hcaman.pdf hca.1 tcaman.pdf tca.1 hca.pdf clean 

tca.1:  tcaman.tex
	latex2man tcaman.tex tca.1

tcaman.pdf:  tcaman.tex
	pdflatex tcaman.tex

hca.1:  hcaman.tex
	latex2man hcaman.tex hca.1

hcaman.pdf:  hcaman.tex
	pdflatex hcaman.tex

hca.pdf:  hca.tex hca.bbl
	pdflatex hca.tex

clean:
	rm -f hcaman.tex~ hcaman.log hcaman.aux tcaman.tex~ tcaman.log tcaman.aux
	rm -f hca.tex~ hca.log hca.aux hca.blg hca.out
