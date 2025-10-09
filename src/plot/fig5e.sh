hicPlotTADs --tracks tracks.ini -o hic_track.pdf --region Chr01:1-5835065

########tracks.ini
'''
[x-axis]
where = top

[hic matrix]
file = hic_corrected.h5
title = Hi-C data
# depth is the maximum distance plotted in bp. In Hi-C tracks
# the height of the track is calculated based on the depth such
# that the matrix does not look deformed
depth = 700000
transform = log1p
file_type = hic_matrix

[tads]
file = ontad/OnTAD_chr15.bed
file_type = domains
border_color = black
color = none
# the tads are overlay over the hic-matrix
# the share-y options sets the y-axis to be shared
# between the Hi-C matrix and the TADs.
overlay_previous = share-y

[spacer]

[genes]
file = Dc.bed
height = 7
title = Genes
style = flybase
fontsize = 10

#[test bedgraph]
#file = geneDensity.bedgraph
#color = #1E90FF
#height = 5
#title = Gene density
#max_value = 7

[test bedgraph]
file = CopiaDensity.bedgraph
color = orange
height = 5
title = Copia
max_value = 10

[test bedgraph]
file = GypsyDensity.bedgraph
color = #D6AFB9
height = 5
title = Gypsy
max_value = 10
'''