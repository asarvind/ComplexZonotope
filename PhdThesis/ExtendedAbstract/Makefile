INPUT_FILE_LONG=summary.tex
INPUT_FILE_SHORT=ShortSummary.tex
SPELL_CHECK=aspell -c

#submission:
#	latex ${INPUT_FILE}
#	latex ${INPUT_FILE}
#	dvips -o main.ps main.dvi
#	ps2pdf main.ps main.pdf

all:
#	latexmk -f main.tex -pdfdvi
	make clean
	ln -s ../../PhdThesis
	ln -s ../stylefiles
	ln -s ../../figures
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE_LONG}
	bibtex summary
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE_LONG}
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE_LONG}
#	open main.pdf &

short_pdf:
#xs	latexmk -f main.tex -pdfdvi
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE_SHORT}
	bibtex ShortSummary
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE_SHORT}
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE_SHORT}
#	open main.pdf &

TEX_FILES := $(wildcard *.tex chapters/*.tex chapters/*/*.tex)

.PHONY: check

check: 
	for texfile in $(TEX_FILES); do \
		$(SPELL_CHECK) $$texfile; \
		done

clean:  
	\rm -f *.pdf
	\rm -f  *.dvi
	\rm -f  *.brf
	\rm -f  *.gz
	\rm -f  *.out
	\rm -f  *.log
	\rm -f  *.bak
	\rm -f  *.aux
	\rm -f  *.blg
	\rm -f  *.bbl
	\rm -f  *.*~
	\rm -f  *.mtc*
	\rm -f  *.maf
	\rm -f  *.toc
	\rm -f  *.ptc
	\rm -f  *.tex.backup
	\rm -f  chapters/*.tex.backup	
	\rm -f  chapters/*/*.tex.backup		
	\rm -f  chapters/*.tex~	
	\rm -f  chapters/*/*.tex~
	\rm -f *~
	\rm -f  Phd*
	\rm -f -f stylefiles
	\rm -f figures
	\rm -f *syntex*
	# git rm -f Figures/*.pdf

upload:
	make clean; git commit -m "cleaned"; git push; git add .; git commit -m "updated"; git push
