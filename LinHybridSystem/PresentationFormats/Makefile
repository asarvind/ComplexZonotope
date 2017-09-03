files1 = Makefile macros.tex settings.tex slides.tex 

bibfiles = ref.bib

files = $(files1)

figures = fig

sources = $(files) $(bibfiles) $(figures)

output: $(sources)
	make clean; pdflatex slides.tex; bibtex slides ; pdflatex slides.tex; pdflatex slides.tex;

pdf: $(files)
	pdflatex slides.tex;

add: 
	make clean; git add $(sources)


commit:
	git commit -m "updated repo"

upload: 
	make add; make commit;  git push

clean:
	rm *.pdf *.aux *.log *~ *.bbl *.blg *.nav *.out *.snm *.toc
