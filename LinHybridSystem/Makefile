output: main.tex
	make clean; pdflatex main.tex; bibtex main; pdflatex main.tex; pdflatex main.tex;



add: 
	make clean; git add Makefile main.tex macros.tex ref.bib SetInclusions.tex PositiveInvariance.tex

commit:
	git commit -m "updated repo"

upload: 
	make add; make commit;  git push

clean:
	rm main.pdf *.aux *.log *~ *.bbl *.blg

