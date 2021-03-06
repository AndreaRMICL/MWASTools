## INTERNAL FUNCTIONS ##

### SK_scores ####
SK_scores = function(pvalue, estimate) {
    logpvalue = -log10(pvalue)

    if (estimate < 0) {
        logpvalue = -logpvalue
    }
    return(logpvalue)
}

### SK_color ####
SK_color = function(score, scale_color = scale_color, alpha_th = alpha_th) {
    th = -log10(alpha_th)
    if (score < (-th)) {
        color = scale_color[2]
    } else if (score > th) {
        color = scale_color[3]
    } else {
        color = scale_color[1]
    }
    return(color)
}

### sig_color ####
sig_color = function(score, alpha_th = alpha_th) {
    th = -log10(alpha_th)
    if (score < (-th)) {
        color = 1
    } else if (score > th) {
        color = 2
    } else {
        color = 0
    }
    return(color)
}


## EXTERNAL FUNCTIONS ##

### MWAS_skylineNMR ####
MWAS_skylineNMR = function(metabo_SE, MWAS_matrix, ref_sample, alpha_th = 0.05,
                           output = "all", xlab = "ppm", ylab1 = "sign*log(pFDR)",
                           ylab2 = "intensity", pch = 20, marker_size = 1,
                           scale_color = c("black", "cornflowerblue", "red"),
                           size_lab = 12, size_axis = 12, xlim = NULL, ylim1 = NULL,
                           ylim2 = NULL, guide_type = "legend", xbreaks = waiver(),
                           xnames = waiver(), ybreaks1 = waiver(), ybreaks2 = waiver(),
                           ynames1 = waiver(), ynames2 = waiver()) {

    ## Check that the input data are correct
    if (class(metabo_SE)[1] != "SummarizedExperiment") {
        stop("metabo_SE must be a SummarizedExperiment object")
    }
    metabo_matrix = t(assays(metabo_SE)$metabolic_data)
    ppm = rownames(metabo_SE)
    if (suppressWarnings(is.na(as.numeric(ppm[1])))) {
        stop("metabo_SE rownames seem not to correspond to a ppm scale")
    } else {
        ppm = as.numeric(ppm)
    }
    if (ref_sample %in% colnames(metabo_SE) == FALSE) {
        stop ("ref_sample is not included in metabo_SE colnames")
    } else {
        metabo_vector = metabo_matrix[ref_sample, ]
    }

    if (!is.matrix (MWAS_matrix) | !is.numeric(MWAS_matrix)) {
        stop ("MWAS_matrix must be a numeric matrix. Check MWAS_stats()")
    }
    if (ncol(MWAS_matrix) < 3 ) {
        stop ("MWAS_matrix does not seem to have the right format. Check MWAS_stats")
    }

    if (sum(is.na(MWAS_matrix))) {
        stop ("NA values in MWAS_matrix are not allowed")
    }
    estimates = MWAS_matrix[,1]
    pvalues = MWAS_matrix[,3]

    if (length(pvalues) != length(ppm)) {
        stop("metabo_SE and MWAS_matrix are not consistent")
    }

    ## Sort xlim in decreasing order (ppm is plotted in inverse scale)
    xlim = suppressWarnings(sort(xlim, decreasing = TRUE))

    if (length(scale_color) != 3) {
        stop("scale_color must have 3 color values")
    }

    scores = unlist(mapply(SK_scores, pvalues, estimates, SIMPLIFY = FALSE))
    assoc = as.numeric(sapply(scores, sig_color, alpha_th = alpha_th))

    nb = sort(as.numeric(unique(assoc)), decreasing = FALSE)
    legend_labels = as.character(nb)
    legend_labels[legend_labels == "0"] = "unchanged"
    legend_labels[legend_labels == "1"] = "downregulated"
    legend_labels[legend_labels == "2"] = "upregulated"
    col_breaks = nb
    col_breaks[col_breaks == 0] = scale_color[1]
    col_breaks[col_breaks == 1] = scale_color[2]
    col_breaks[col_breaks == 2] = scale_color[3]

    # Adjust scale for spectrum

    if (class(ybreaks2) == "waiver" & class(ynames2) == "waiver") {#default values
        factor = nchar(round(max(metabo_vector), 0)) - 1
          if (factor > 1) {
            metabo_vector = metabo_vector / 10^factor
            ylab2 = paste(ylab2, " (x", 10, "^", factor, ")", sep = "")
            #ylab2 = bquote(~ .(ylab2) ~ '(' ~ '10'^(.factor))
          }
        }

    data_SK = data.frame(ppm = ppm, scores = scores, assoc = assoc, metabo_vector = metabo_vector)

    # Do plots

    figure_SK = ggplot(data_SK, aes(ppm, scores, color = assoc)) +
                geom_point(shape = pch, size = marker_size) +
                scale_colour_gradientn(colours = col_breaks, breaks = nb, space = "Lab",
                                       guide = guide_type, labels = legend_labels) +
                scale_x_reverse(limits = xlim, breaks = xbreaks, labels = xnames) +
                scale_y_continuous(limits = ylim1, breaks = ybreaks1, labels = ynames1) +
                theme_bw() + labs(x = xlab, y = ylab1) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      axis.text = element_text(size = size_axis),
                      axis.title = element_text(size = size_lab, vjust = 0)) +
                #geom_hline(yintercept = 0, col = scale_color[1]) +
                      {
                        if (1 %in% nb) geom_hline(yintercept = +log10(alpha_th),
                                                  col = scale_color[2],
                                                  linetype="dashed")
                      } + {
                        if (2 %in% nb) geom_hline(yintercept = -log10(alpha_th),
                                                  col = scale_color[3],
                                                  linetype="dashed")
                      }


    figure_spectrum = ggplot(data_SK, aes(ppm, metabo_vector, color = assoc)) +
                      geom_line() +
                      scale_colour_gradientn(colours = col_breaks, breaks = nb, space = "Lab",
                                             guide = guide_type, labels = legend_labels) +
                      scale_x_reverse(limits = xlim, breaks = xbreaks, labels = xnames) +
                      scale_y_continuous(limits = ylim2, breaks = ybreaks2, labels = ynames2) +
                      theme_bw() + labs(x = xlab, y = ylab2) +
                      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                            axis.text = element_text(size = size_axis),
                            axis.title = element_text(size = size_lab, vjust = 0))

    if (output == "skyline") {
        plot(figure_SK)
        return(figure_SK)

    } else if (output == "spectrum") {
        plot(figure_spectrum)
        return(figure_spectrum)

    }  else {
      figure_SK_gtable = ggplot_gtable(ggplot_build(figure_SK))
      figure_spectrum_gtable = ggplot_gtable(ggplot_build(figure_spectrum))
      maxWidth = unit.pmax(figure_SK_gtable$widths[2:3], figure_spectrum_gtable$widths[2:3])

      figure_SK_gtable$widths[2:3] = maxWidth
      figure_spectrum_gtable$widths[2:3] = maxWidth

      biplot = arrangeGrob(figure_SK_gtable, figure_spectrum_gtable)
      plot(biplot)
      return(biplot)
    }
}


### MWAS_barplot ####
MWAS_barplot = function(MWAS_matrix, alpha_th = 0.05, width = NULL, scale_color
                        = c("darkgray", "cornflowerblue", "firebrick1"),
                        legend_labs = c("unchanged", "downregulated", "upregulated"),
                        ylab = "sign*log(pFDR)", size_yaxis = 12, size_ylab = 12,
                        size_names = 10, angle_names = 45, sort = TRUE) {

    ## Check that the input data are correct
    if (!is.matrix (MWAS_matrix) | !is.numeric(MWAS_matrix)) {
      stop ("MWAS_matrix must be a numeric matrix. Check MWAS_stats()")
    }
    if (ncol(MWAS_matrix) < 3 ) {
      stop ("MWAS_matrix does not seem to have the right format. Check MWAS_stats")
    }
    if (sum(is.na(MWAS_matrix))) {
      stop ("NA values in MWAS_matrix are not allowed")
    }

    estimates = as.numeric(MWAS_matrix[,1])
    pvalues = as.numeric(MWAS_matrix[,3])

    metabo_ids = rownames(MWAS_matrix)

    if (length(pvalues) != length(metabo_ids)) {
        stop("metabo_ids length must be consistent with MWAS_matrix dimension")
    }

    ## Get scores and color scale
    scores = unlist(mapply(SK_scores, pvalues, estimates, SIMPLIFY = FALSE))
    names(scores) = metabo_ids
    scores_col = sapply(scores, SK_color, scale_color = scale_color, alpha_th = alpha_th)

    col_ref = scores_col
    col_ref[col_ref == scale_color[1]] = legend_labs[1]
    col_ref[col_ref == scale_color[2]] = legend_labs[2]
    col_ref[col_ref == scale_color[3]] = legend_labs[3]

    scores_col_values = unique(scores_col)
    names(scores_col_values) = unique(col_ref)
    assoc = factor(col_ref)

    scoresM = data.frame(assoc = assoc, metabo_ids = factor(metabo_ids), scores = scores)

    ## Set positions
    if(sort == TRUE){
        positions = names(sort(scores))
    } else{
        positions = names(scores)
    }

    ## Do barplot
    figure_barplot = ggplot(data = scoresM, aes(x = metabo_ids, y = scores, fill = assoc)) + {
        if (is.null(width)) geom_bar(stat = "identity", position = position_dodge())
        else geom_bar(stat = "identity", position = position_dodge(), width = width)
      } +
        scale_x_discrete(limits = positions) +
        scale_fill_manual(values = scores_col_values) +
        theme_bw() +
        labs(x = " ", y = ylab) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text = element_text(size = size_yaxis),
              axis.text.x = element_text(angle = angle_names, hjust = 1, size = size_names),
              axis.title = element_text(size = size_ylab, vjust = 0 ))
    plot(figure_barplot)
    return(figure_barplot)
}
