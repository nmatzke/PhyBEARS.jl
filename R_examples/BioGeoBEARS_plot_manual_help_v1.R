res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, show.tip.label=TRUE, tipcol="white", xlims=c(0,8))

# From: 
# https://stackoverflow.com/questions/43378199/italics-and-regular-text-in-phylogeny-tip-labels-in-r
for(i in seq_along(tr$tip.label)){
    tip <- unlist(strsplit(tr$tip.label[i], "_"))

    tiplabels(
        bquote(italic(.(paste(tip[-length(tip)], collapse = ' '))) ~ .(tip[length(tip)])),
        tip = i, adj = c(0, 0.5), frame = "n", bg = "white", offset=0.45
    )

}
