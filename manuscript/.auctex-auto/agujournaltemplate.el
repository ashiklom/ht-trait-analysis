(TeX-add-style-hook
 "agujournaltemplate"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("agujournal2019" "draft")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("trackchanges" "inline")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "agujournal2019"
    "agujournal201910"
    "url"
    "lineno"
    "trackchanges"
    "soul")
   (LaTeX-add-bibliographies
    " enter your bibtex bibliography filename here "))
 :latex)

