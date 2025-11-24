#!/usr/bin/env gnuplot
# ==========================================
# 3D Plot Script - Single Configuration
# ==========================================

filename = '../src/3D_coordinates.txt'

# Lecture du nombre d'atomes depuis la première ligne
firstline = system("head -n 1 ".filename)
ncolonnes = words(firstline)
n_atomes = int(ncolonnes / 3)

# Réglage général
param = 10   # définit l'étendue des axes
set terminal pngcairo size 800,800 enhanced font 'Helvetica,10'
set output '3D_plot.png'

# Palette rose-violet lissée
set palette defined (\
    0  0.1647 0.0314 0.2706, \
    1  0.4157 0.0196 0.4471, \
    2  0.6706 0.0510 0.5490, \
    3  0.8902 0.2902 0.6235, \
    4  0.9725 0.6471 0.7608 )

set palette color

# Paramètres du graphe
unset key
set size ratio -1
set xlabel "x" textcolor rgb "black"
set ylabel "y" textcolor rgb "black"
set zlabel "z" textcolor rgb "black"
set grid lc rgb "gray"
set xrange [-param:param]
set yrange [-param:param]
set zrange [-param:param]
set title "Configuration unique" textcolor rgb "black"
set border lc rgb "black"
set tics textcolor rgb "black"
set view 90, 30, 1.0, 1.0  # Vue 3D

# --- 3D plot ---
splot for [j=(n_atomes-1):0:-1] filename every ::0::0 using (column(3*j+1)):(column(3*j+2)):(column(3*j+3)):(j) \
 with points pt 7 ps 1 lc palette

unset output
