maintext.pdf: maintext.md bibliography.bib
	pandoc -V fontsize=11pt --include-in-header=header.tex -F pandoc-crossref maintext.md -o maintext.pdf --bibliography bibliography.bib