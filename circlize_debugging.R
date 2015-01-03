# Useful for adjusting circlize's circos plotting parameters
# circos.par("cell.padding") # [1] 0.02 0.00 0.02 0.00; bottom, left, top, right; 1st & 3rd are % radius, 2nd & 4th are degrees
# circos.par("track.height") # [1] 0.2 , i.e. 20% of radius
# circos.par("track.margin") # [1] 0.01 0.01 # bottom and top margins, left and right are controlled by gap.degree
# circos.par("gap.degree") # [1] 1
# circos.info()

# Useful for debugging inside trackPlotRegions
#         print(get.cell.meta.data("sector.index"))
#         print(get.cell.meta.data("xlim"))
#         print(get.cell.meta.data("ylim"))
#         print(region)
#         print(value)