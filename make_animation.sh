#!/bin/tcsh
#foreach x (TNdiff_??.eps DENratio_??.eps DENratio_alt_??.eps)
#  convert -density 96 $x $x:r.gif
#end

gifsicle -U --disposal=previous --colors=256 --delay=10 --loop=0 -O2 TNdiff_??.gif --delay=400 TNdiff_72.gif > tutorial_01.gif
gifsicle -U --disposal=previous --colors=256 --delay=10 --loop=0 -O2 DENratio_??.gif --delay=400 DENratio_72.gif > tutorial_02.gif
gifsicle -U --disposal=previous --colors=256 --delay=10 --loop=0 -O2 DENratio_alt_??.gif --delay=400 DENratio_alt_72.gif > tutorial_03.gif

#--transparent="#FFFFFF"
