## EXTERNAL FUNCTION ##

plot_spectraNMR = function(metabo_matrix, ppm, type = "l", lty = 1,
                           xlab = "ppm", ylab = "intensity", xlim = NULL, ...) {

    ## Check that input data are correct
    if (is.matrix(metabo_matrix) == FALSE & is.vector(metabo_matrix) == FALSE) {
        stop("metabo_matrix needs to be a numeric matrix or vector")
    }
    if (is.vector(metabo_matrix)) {
        if (length(metabo_matrix) != length(ppm)) {
            stop("ppm length is not consistent with metabo_matrix length")
        }
        metabo_matrix = matrix(metabo_matrix, ncol = length(ppm))
    }
    if (!is.numeric(ppm)) {
        stop ("ppm must be a numeric vector")
    }

    if (is.null(xlim)) {
        xlim = rev(range(ppm))
    } else {
        ind1 = grep(xlim[2], ppm, fixed = TRUE)[1]
        ind2 = grep(xlim[1], ppm, fixed = TRUE)[1]
        ppm = ppm[ind1:ind2]
        metabo_matrix = metabo_matrix[, ind1:ind2]
    }

    metabo_matrix = t(metabo_matrix)
    matplot(ppm, metabo_matrix, type = type, xlim = xlim, xlab = xlab,
        ylab = ylab, lty = lty, ...)
}

