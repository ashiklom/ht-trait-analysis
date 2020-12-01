(TeX-add-style-hook
 "manuscript"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "authblk")
   (LaTeX-add-labels
    "sec:introduction"
    "sec:methods"
    "sec:results"))
 :latex)

