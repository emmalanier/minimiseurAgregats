#!/bin/bash

# Nettoyage des anciens résultats
rm -rf plots_batch_* batch_*.mp4 animation_finale.mp4 liste_videos.txt

#
file="Donnees_minimisation_cov.txt"

nb_lignes=$(wc -l < "$file")

nb_batch=$(( ($nb_lignes % 200) + 1 ))


# Génération les dossiers à l'avance
for i in $(seq 0 $((nb_batch - 1))); do
    mkdir -p plots_batch_$(printf "%04d" $i)
done

# Partie Gnuplot
gnuplot gnuplot_cov.gp

# Génération d'un GIF par batch (si le batch n'est pas vide), GIF stocké dans le dossier parent
for dir in plots_batch_*/; do
    if compgen -G "${dir}/*.png" > /dev/null; then
        echo "Génération du MP4 pour $dir"
        batchname=$(basename "$dir")
        ffmpeg -framerate 40 -pattern_type glob -i "${dir}/plot_*.png" \
               -c:v libx264 -pix_fmt yuv420p -y "${batchname}.mp4"
    fi
done

echo "Vidéos générées"

# Génération de la liste des vidéos
for f in plots_batch_*.mp4; do
    echo "file '$f'" >> liste_videos.txt
done

# Assemblage en une vidéo finale
ffmpeg -f concat -safe 0 -i liste_videos.txt -c copy animation_finale.mp4

echo "Vidéo finale générée"

