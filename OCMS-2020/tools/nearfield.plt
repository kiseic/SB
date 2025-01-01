#!/usr/bin/gnuplot -persist

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.

#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 6    last modified September 2014
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2014
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')

set terminal png font "Sans,12"
set output "dat.nearfield.png"

set format x "% g"
set format y "% g"

set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000

set offsets 0, 0, 0, 0
set pointsize 1

set size ratio 0 1,1
set origin 0,0

set xtics autofreq  norangelimit
set ytics autofreq  norangelimit

set title "" 

set xlabel "The arc length along the cavity boundary" 
set xrange [ * : * ] noreverse nowriteback

set ylabel "Intensity [arb. units]" 
set yrange [ * : * ] noreverse nowriteback

plot "dat.nearfield" w lp lc "red"

