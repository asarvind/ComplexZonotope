symnames = stylefiles figures fig lib.bib thesis

all:
	make clean
	make sym
	make pdf
	bibtex slides
	make pdf
	make pdf


sym: 
	ln -s ../../PhdThesis/stylefiles/ stylefiles;	
	ln -s ../../figures/ figures;
	ln -s ../malthyJan2018/fig fig;
	ln -s ../../PhdThesis/myLib.bib lib.bib
	ln -s ../../PhdThesis/ thesis

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
	rm -f .nfs*
