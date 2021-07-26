maintext.pdf: maintext.md bibliography.bib
	pandoc maintext.md -o maintext.pdf --bibliography bibliography.bib