# plot_batched.gp

filename = 'Donnees_minimisation_vdw.txt' 
stats filename nooutput
nblignes = STATS_records

# Lecture du nombre d'atomes depuis la première ligne
firstline = system("head -n 1 ".filename)
ncolonnes = words(firstline)
n_atomes = int(ncolonnes / 2)

# Réglage général
param = 10 # <- Va définir la distance entre l'extrema d'un axe et le centre du repère, peut avoir besoin d'être adapté
set terminal pngcairo size 800,800 enhanced font 'Helvetica,10'

# Palette rose-violet lissée
set palette defined (\
    0  0.1647 0.0314 0.2706, \
    1  0.4157 0.0196 0.4471, \
    2  0.6706 0.0510 0.5490, \
    3  0.8902 0.2902 0.6235, \
    4  0.9725 0.6471 0.7608 )

set palette color

# Nombre d'images par dossier
batch_size = 200

# Plot boucle


do for [i=0:nblignes-1] {
    batch = int(i / batch_size)
    dossier = sprintf("plots_batch_%04d", batch)
    set output sprintf("%s/plot_%04d.png", dossier, i)

    unset key
    set size ratio -1
    set xlabel "x" textcolor rgb "black"
    set ylabel "y" textcolor rgb "black"
    set grid lc rgb "gray"
    set xrange [-param:param]
    set yrange [-param:param]
    set title sprintf("Configuration %d", i) textcolor rgb "black"
    set border lc rgb "black"
    
    set tics textcolor rgb "black"


    plot for [j=(n_atomes-1):0:-1] filename every ::i::i using (column(2*j+1)):(column(2*j+2)):(j) \
     with points pt 7 ps 1 lc palette
}

unset output

