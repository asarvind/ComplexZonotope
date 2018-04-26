symnames = stylefiles figures

all:
	make clean
	make sym
	make pdf
	bibtex slides
	make pdf
	make pdf


sym: 
	ln -s ../../PhdThesis/stylefiles/ stylefiles	
	ln -s ../../figures/ figures

pdf: 
	pdflatex slides.tex

add: 
	make clean; git add ./


commit:
	git commit -m "updated repo"

upload: 
	make add; make commit;  git push

clean:
	rm -f *.pdf *.aux *.log *~ *.bbl *.blg *.nav *.out *.snm *.toc
	rm -f $(symnames)
	rm -f *#*
