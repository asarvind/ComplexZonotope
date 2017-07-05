source_files1 = IEEEtran.cls ieeeconf.cls intro.tex main.tex mprt.bib Makefile aliascnt.sty conclusion.tex domain.bib

source_files2 = fig/ lis-exp.tex lis.tex litr.tex llncs-files/ llncs.cls llncs2e.zip llncsdoc.sty ref.bib rellis.tex 

source_files3 = relswitched.tex remreset.sty sw-exp.tex system.tex temp-comp-zon.tex

source_files4 = packages.tex settings.tex macros.tex

source_files = $(source_files1) $(source_files2) $(source_files3) $(source_files4)



output: $(source_files)
	make clean; pdflatex main.tex; bibtex main; pdflatex main.tex; pdflatex main.tex;

pdf: $(source_files)
	pdflatex main.tex;



add: 
	make clean; git add $(source_files)


commit:
	git commit -m "updated repo"

upload: 
	make add; make commit;  git push

clean:
	rm *.pdf *.aux *.log *~ *.bbl *.blg
