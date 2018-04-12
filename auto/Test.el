(TeX-add-style-hook
 "Test"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("report" "11pt" "spanish")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("babel" "spanish") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "report"
    "rep11"
    "babel"
    "inputenc"
    "mathtools"
    "physics"
    "qcircuit"
    "graphicx"
    "float"
    "dsfont"
    "tkz-graph")
   (TeX-add-symbols
    '("qwxo" ["argument"] 1)))
 :latex)

