# nechajte iba jeden z main.pdf a main-en.pdf
all: main-en.pdf

main.pdf: main.tex *.tex *.bib images/*
	pdflatex main
	biber main
	pdflatex main
	pdflatex main


main-en.pdf: main-en.tex *.tex *.bib images/*
	pdflatex main-en
	biber main-en
	pdflatex main-en
	pdflatex main-en
