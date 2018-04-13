(TeX-add-style-hook
 "Avance"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "xetex" "mathserif" "serif")))
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "physics"))
 :latex)

