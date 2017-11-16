INPUT_FILE=main.tex
SPELL_CHECK=aspell -c

#submission:
#	latex ${INPUT_FILE}
#	latex ${INPUT_FILE}
#	dvips -o main.ps main.dvi
#	ps2pdf main.ps main.pdf

all:
#	latexmk -f main.tex -pdfdvi
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}
	bibtex main 
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}
#	open main.pdf &

pdf:
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}

TEX_FILES := $(wildcard *.tex chapters/*.tex chapters/*/*.tex)

.PHONY: check

check: 
	for texfile in $(TEX_FILES); do \
		$(SPELL_CHECK) $$texfile; \
		done

clean:  \rm -f *.pdf
	\rm -f  *.dvi
	\rm -f  *.brf
	\rm -f  *.gz
	\rm -f  *.out
	\rm -f  *.log
	\rm -f  *.bak
	\rm -f	*.aux
	\rm -f	*.blg
	\rm -f	*.bbl
	\rm -f	*.*~
	\rm -f  *.mtc*
	\rm -f  *.maf
	\rm -f  *.toc
	\rm -f  *.ptc
	\rm -f  *.tex.backup
	\rm -f  chapters/*.tex.backup	
	\rm -f  chapters/*/*.tex.backup		
	\rm -f  chapters/*.tex~	
	\rm -f  chapters/*/*.tex~		
	#\rm -f Figures/*.pdf

upload:
	make clean; git add .; git commit -m "updated"; git push