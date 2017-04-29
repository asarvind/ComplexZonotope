output: main.tex
	make clean; pdflatex main.tex; bibtex main; pdflatex main.tex; pdflatex main.tex;

pdf: main.tex
	pdflatex main.tex;



add: 
	make clean; git add Makefile main.tex macros.tex ref.bib prelims.tex zonotope.tex cztope.tex PositiveInvariance.tex benchmarks.tex 

add1:
	git add intro.tex contribution.tex conclusion.tex appendix.tex

commit:
	git commit -m "updated repo"

upload: 
	make add; make add1; make commit;  git push

clean:
	rm main.pdf *.aux *.log *~ *.bbl *.blg

