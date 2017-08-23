(TeX-add-style-hook
 "afm-cs"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("ieeeconf" "letterpaper" "10pt" "conference")))
   (TeX-run-style-hooks
    "latex2e"
    "ieeeconf"
    "ieeeconf10"
    "array"
    "amsmath"
    "amssymb"
    "bbold"
    "mathtools"
    "graphicx"
    "epstopdf"
    "caption"
    "subcaption"
    "color"
    "multirow"
    "url"
    "booktabs"
    "placeins"
    "cases"
    "empheq"
    "picture"
    "calc"
    "xspace")
   (TeX-add-symbols
    '("pic" 2)
    "na"
    "eg"
    "ie"
    "Dm")
   (LaTeX-add-labels
    "sec:pack")))

