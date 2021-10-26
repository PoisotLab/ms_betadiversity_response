ALL: reviews/responses.pdf reviews/diff.pdf

reviews/responses.pdf: reviews/01_PCI_round1.md
	pandoc $< -o $@

reviews/diff.tex: reviews/01_draft.tex reviews/01_revised.tex
	latexdiff -t BOLD -s COLOR $^ > $@

reviews/diff.pdf: reviews/diff.tex
	pdflatex $< 
	mv diff* reviews/
