maintext.pdf: maintext.md bibliography.bib
	pandoc -F pandoc-crossref maintext.md -o maintext.pdf --bibliography bibliography.bib