INPUT_FILE=main.tex
SPELL_CHECK=aspell -c

# Symbolic link

SYMBOLIC_LINK_NAME1=reflis.bib
SYMBOLIC_LINK_NAME2=fig

#submission:
#	latex ${INPUT_FILE}
#	latex ${INPUT_FILE}
#	dvips -o main.ps main.dvi
#	ps2pdf main.ps main.pdf

all:
#	latexmk -f main.tex -pdfdvi
	make clean
	make sym
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}
	bibtex main 
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}
#	open main.pdf &

pdf:
	pdflatex -shell-escape -synctex=1 ${INPUT_FILE}

TEX_FILES := $(wildcard *.tex chapters/*.tex chapters/*/*.tex)

sym:
	\rm -f ${SYMBOLIC_LINK_NAME1} ${SYMBOLIC_LINK_NAME2}
	ln -s ~/Dropbox/Invariants/complex-zonotope/ref.bib ${SYMBOLIC_LINK_NAME1}
	ln -s ../figures ${SYMBOLIC_LINK_NAME2}


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
	\rm -f  ${SYMBOLIC_LINK_NAME1}
	\rm -f  ${SYMBOLIC_LINK_NAME2}
	# git rm -f Figures/*.pdf

upload:
	make clean; git commit -m "cleaned"; git push; git add .; git commit -m "updated"; git push
