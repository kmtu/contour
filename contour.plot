#!/bin/bash
# This is a script to generate gnuplot codes for contour plot and draw it
if [ "$1" == "" ]; then
   echo "Please give me the data filename for plotting!"
   exit 1
fi
OUTFILE=$1".ps"

cat <<EOF > contour.gpt
set terminal postscript color
set output "$OUTFILE"
set pm3d map
set size square
set palette rgbformulae 22,13,10
#set xrange [0:19]
#set yrange [0:19]
#set cbrange [0:500]
splot "$1"
EOF

gnuplot contour.gpt
exit 0
