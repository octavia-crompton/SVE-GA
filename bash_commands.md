ssh octavia@hpc.brc.berkeley.edu

transfer to savio:

    scp -r model octavia@dtn.brc.berkeley.edu:/global/home/users/octavia/model


transfer from savio:

    scp -r  octavia@dtn.brc.berkeley.edu:/global/home/users/octavia/model model


sublime alias:
alias subl="/Applications/Sublime\ Text.app/Contents/SharedSupport/bin/subl"


pattern="/Users/octavia/Dropbox/SVEpattern"

hide .pyc files:

    find -x /path/to/folder -name \*.pyc -exec chflags hidden {} +
    find -x ~/Dropbox/ -name \*.pyc -exec chflags hidden {} +


vim dryR.for -c ":hardcopy > abc.ps" -c ":q"
ps2pdf abc.ps ~/Downloads/dryR.pdf
