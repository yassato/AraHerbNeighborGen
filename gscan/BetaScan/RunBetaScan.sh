fpath="../output/"
fprefix="BetaScanChr5"
input=$fpath$fprefix"in.txt"
output=$fpath$fprefix"out.txt"

# BetaScan.py is available via https://github.com/ksiewert/BetaScan
python2 ./BetaScan/BetaScan.py -i $input -o $output
