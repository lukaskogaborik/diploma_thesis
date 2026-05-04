# nechajte iba jeden z main.pdf a main-en.pdf
all: main-en.pdf

main.pdf: main.tex *.tex *.bib images/*
	pdflatex escape main
	biber main
	pdflatex main
	pdflatex main


main-en.pdf: main-en.tex *.tex *.bib images/*
	pdflatex --shell-escape main-en
	biber main-en
	pdflatex --shell-escape main-en
	pdflatex --shell-escape main-en
