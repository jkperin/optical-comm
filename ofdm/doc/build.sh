latex --shell-escape main.tex

latex -shell-escape -halt-on-error-interaction=batchmode -jobname "figs/main-figure0" "\def \tikzexternalrealjob {main}\input {main}"; dvips -o "figs/main-figure0".ps "figs/main-figure0".dvi; ps2eps "figs/main-figure0".ps
latex -shell-escape -halt-on-error-interaction=batchmode -jobname "figs/main-figure1" "\def \tikzexternalrealjob {main}\input {main}"; dvips -o "figs/main-figure1".ps "figs/main-figure1".dvi; ps2eps "figs/main-figure1".ps
latex -shell-escape -halt-on-error-interaction=batchmode -jobname "figs/main-figure2" "\def \tikzexternalrealjob {main}\input {main}"; dvips -o "figs/main-figure2".ps "figs/main-figure2".dvi; ps2eps "figs/main-figure2".ps
latex -shell-escape -halt-on-error-interaction=batchmode -jobname "figs/main-figure3" "\def \tikzexternalrealjob {main}\input {main}"; dvips -o "figs/main-figure3".ps "figs/main-figure3".dvi; ps2eps "figs/main-figure3".ps
latex -shell-escape -halt-on-error-interaction=batchmode -jobname "figs/main-figure4" "\def \tikzexternalrealjob {main}\input {main}"; dvips -o "figs/main-figure4".ps "figs/main-figure4".dvi; ps2eps "figs/main-figure4".ps
latex -shell-escape -halt-on-error-interaction=batchmode -jobname "figs/main-figure5" "\def \tikzexternalrealjob {main}\input {main}"; dvips -o "figs/main-figure5".ps "figs/main-figure5".dvi; ps2eps "figs/main-figure5".ps
