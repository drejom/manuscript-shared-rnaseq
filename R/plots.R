# plot_libsize() stolen from:
# https://github.com/elsayed-lab/hpgltools/blob/6a9674c02a771368f841b48041db88c89057d0c7/R/plot_bar.r
#' Make a ggplot graph of library sizes.
#'
#' It is often useful to have a quick view of which samples have more/fewer
#' reads.  This does that and maintains one's favorite color scheme and tries to
#' make it pretty!
#'
#' @param data Expt, dataframe, or expressionset of samples.
#' @param condition vector of sample condition names.
#' @param colors Color scheme if the data is not an expt.
#' @param text Add the numeric values inside the top of the bars of the plot?
#' @param order  Explicitly set the order of samples in the plot?
#' @param title Title for the plot.
#' @param yscale Whether or not to log10 the y-axis.
#' @param expt_names  Design column or manually selected names for printing sample names.
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param ... More parameters for your good time!
#' @return a ggplot2 bar plot of every sample's size
#' @seealso \pkg{ggplot2}
#'  \code{\link[ggplot2]{geom_bar}} \code{\link[ggplot2]{geom_text}}
#'  \code{\link{prettyNum}} \code{\link[ggplot2]{scale_y_log10}}
#' @examples
#' \dontrun{
#'  libsize_plot <- plot_libsize(expt=expt)
#'  libsize_plot  ## ooo pretty bargraph
#' }
#' @export
plot_libsize <- function(data, condition=NULL, colors=NULL,
                         text=TRUE, order=NULL, title=NULL,  yscale=NULL,
                         expt_names=NULL, label_chars=10,
                         ...) {
  arglist <- list(...)
  if (is.null(text)) {
    text <- TRUE
  }
  
  ## In response to Keith's recent comment when there are more than 8 factors
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }
  design <- NULL
  data_class <- class(data)[1]
  if (data_class == "expt") {
    design <- pData(data)
    condition <- design[["condition"]]
    colors <- data[["colors"]]
    mtrx <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    design <- pData(data)
    condition <- design[["condition"]]
    mtrx <- exprs(data)
  } else if (data_class == "matrix" | data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    mtrx <- as.data.frame(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }
  
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(ncol(mtrx),
                               chosen_palette))(ncol(mtrx))
  }
  
  ## Get conditions
  if (is.null(condition)) {
    stop("Missing condition label vector.")
  }
  
  values <- as.numeric(mtrx)
  integerp <- all.equal(values, as.integer(values))
  
  colors <- as.character(colors)
  sum <- NULL
  
  if (!is.null(expt_names) & class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      colnames(mtrx) <- make.names(design[[expt_names]], unique=TRUE)
    } else {
      colnames(mtrx) <- expt_names
    }
  }
  if (!is.null(label_chars) & is.numeric(label_chars)) {
    colnames(mtrx) <- abbreviate(colnames(mtrx), minlength=label_chars)
  }
  
  libsize_df <- data.frame("id" = colnames(mtrx),
                           "sum" = colSums(mtrx),
                           "condition" = condition,
                           "colors" = as.character(colors))
  summary_df <- data.table::setDT(libsize_df)[, list("min"=min(sum),
                                                     "1st"=quantile(x=sum, probs=0.25),
                                                     "median"=median(x=sum),
                                                     "mean"=mean(sum),
                                                     "3rd"=quantile(x=sum, probs=0.75),
                                                     "max"=max(sum)),
                                              by="condition"]
  libsize_plot <- plot_sample_bars(libsize_df, condition=condition, colors=colors,
                                   text=text, order=order, title=title, integerp=integerp,
                                   yscale=yscale, ...)
  ##libsize_plot <- plot_sample_bars(libsize_df, condition=condition, colors=colors,
  ##                                 text=text, order=order, title=title, integerp=integerp,
  ##                                 yscale=yscale)
  retlist <- list(
    "plot" = libsize_plot,
    "table" = libsize_df,
    "summary" = summary_df)
  return(retlist)
}




#' Generates module-specific expression profile plots 
#'
#' Generates a grid of expression profile plots for each module, up to a
#' specified limit. If there are more modules than the maximum allowed, then
#' a randomly sub-sampled collection of modules will be plotted.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts_long long format count table with columns: "id", "condition",
#' "expression", and "cluster"
#' @param module_colors module assisgnments for each gene
#' @param module_order order to use when plotting modules
#' @param max_plots maximum number of plots to display 
#' @param modules_per_plot number of modules to include in each plot (default: 9)
#' @param ncols Number of columns to display for grid plots
#'
#' @return None
module_expression_profile_plot <- function (counts_long, module_colors,
                                            module_order=NA, max_plots=50,
                                            modules_per_plot=9,
                                            ncols=3) {
    # determine order of plots
    if (any(is.na(module_order))) {
        module_order <- unique(module_colors)
    }

    num_modules <- length(module_order)
    max_modules_to_plot <- max_plots * modules_per_plot

    # Select modules to display
    if (num_modules > (max_plots * modules_per_plot)) {
        module_order <- module_order[approx(1:num_modules, n=max_modules_to_plot)$x]
    }

    # Split modules into groups of n or less for plotting
    module_groups <- split(module_order, ceiling(seq_along(module_order) /
                                                 modules_per_plot))

    for (i in head(names(module_groups), max_plots)) {
        plots <- list()
        counter <- 1
        for (color in module_groups[[i]]) {
            # choose color
            #if (length(module_order) <= 435) {
            #    line_color <- color
            #} else {
            #    line_color <- '#333333'
            #}
            #line_color <- color

            # For now, use grey for default line color; other colors may be
            # difficult to see.
            line_color <- '#212121'

            plots[[counter]] <- module_profile_plot(counts_long, color,
                                                    line_color=line_color)
            counter <- counter + 1
        }
        print(do.call("grid.arrange", c(plots, ncol=ncols)))
    }

}

#' Generates an gene expression plot for a single co-expression module.
#'
#' @param counts_long Data frame of gene counts in long format
#' @param module Name of module to be plotted
#' @param expr_var Name of column in input data containing expression values
#' @param line_width Line width to use for plotting
#' @param line_color Line color to use in plot
#' @param font_color Font color to use for plt
#' @param highlight_group Optional set of genes to be highlighted in a
#'        different color
#' @param scale_color_values Colors to use for highlighted genes, passed to
#'        scale_color_manual function.
#' @param scale_size_values Sizes to use for highlighted genes, passed to
#'        scale_size_manual function.
#' @param scale_linetype_values linetypes to use for highlighted genes, passed to
#'        scale_linetype_manual function.
#' @param include_title Whether or not to include a plot title
#' @param ylimits Options two-value vector containining min and max y-limits
#' @param xlabel Label for x-axis.
#' @param ylabel Label for y-axis.
#'
#' @return ggplot instance
module_profile_plot <- function(counts_long, module, expr_var='expression',
                                line_width=0.2, 
                                line_color='#333333',
                                font_color='#333333', 
                                highlight_group=NULL,
                                scale_color_values=NULL,
                                scale_size_values=NULL,
                                scale_linetype_values=NULL,
                                include_title=TRUE,
                                ylimits=NULL,
                                xlabel='Condition',
                                ylabel='Expression level (log2-CPM)') {
    # Get the module genes
    cluster_genes <- counts_long[counts_long$cluster == module,]

    # Determine average variance of the module genes (used in title)
    cluster_genes_wide <- dcast(cluster_genes, gene_id ~ condition,
                                value.var=expr_var)
    cluster_var <- median(
        apply(cluster_genes_wide[,2:ncol(cluster_genes_wide)], 1, var)
    )

    # choose scale
    if (is.numeric(counts_long$condition)) {
        x_breaks <- c(0, sort(unique(counts_long$condition)))
        plot_scale <- scale_x_continuous(breaks=x_breaks)
    } else {
        # use same order 
        plot_scale <- scale_x_discrete(
            limits=as.character(unique(counts_long$condition)),
            expand=c(0, 0)
        )
    }

    # base plot
    plt <- ggplot(data=cluster_genes,
            aes_string(x='condition', y=expr_var, group='gene_id'))

    if (include_title) {
        plt <- plt + ggtitle(sprintf("%s (n=%d)", module, 
                                     nrow(cluster_genes_wide)))
    }

    # line plot
    if (!is.null(highlight_group)) {
        # dynamic size and color (can use to highlight specific genes)
        plt <- plt + geom_line(aes_string(colour=highlight_group, 
                                          order=highlight_group,
                                          size=highlight_group,
                                          linetype=highlight_group)) +
                     scale_color_manual(values=scale_color_values) + 
                     scale_size_manual(values=scale_size_values) +
                     scale_linetype_manual(values=scale_linetype_values)
    } else {
        # static color and size
        plt <- plt + geom_line(size=I(line_width), colour=line_color) + 
            scale_size(guide='none')
    }

    # y-axis limits (optional)
    if (!is.null(ylimits)) {
        plt <- plt + scale_y_continuous(limits=ylimits)
    }

    # remaining options
    plt <- plt +
        xlab(xlabel) +
        ylab(ylabel) +
        plot_scale

        #theme(axis.text=element_text(colour=font_color), 
        #      plot.margin=unit(c(5.5, 12.5, 5.5, 7.5), "pt"))

    return(plt)
}

#' Plot module subnetwork
#'
#' Generates a network plot for the genes in a single module
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param module_color module color
#' @param module_gene_ids List of genes contained in the module
#' @param adj_matrix Adjacency matrix
#'
#' @return None
plot_module_network <- function(module_color, module_gene_ids, adj_matrix) {
    # Color network plot according to module color
    rgb_col <- col2rgb(module_color)
    node_color <- rgb(rgb_col[1], rgb_col[2], rgb_col[3], maxColorValue=255)

    # igraph plot options
    igraph.options(vertex.size=10, 
                vertex.color=node_color, 
                edge.arrow.mode='-',
                edge.color='#333333')

    genes_sim <- similarity_matrix[rownames(similarity_matrix) %in% module_gene_ids,
                                colnames(similarity_matrix) %in% module_gene_ids]

    # grab up to the 150 most highly correlated gene pairs in the module
    sorted_correlations <- sort(abs(genes_sim), decreasing=TRUE)
    num_to_include <- min(150, length(sorted_correlations))
    ind <- which(genes_sim >= sorted_correlations[num_to_include], arr.ind=TRUE)

    # create an edge list
    edge_list1 <- matrix(
        c(rownames(genes_sim)[as.numeric(ind[,'row'])],
        colnames(genes_sim)[as.numeric(ind[,'col'])]),
        ncol=2
    )

    # Convert to an edge list
    g1 <- graph.edgelist(edge_list1)

    # force-directed layout
    g1 <- permute.vertices(g1, rank(V(g1)$name))
    coords <- layout.fruchterman.reingold(g1)

    # plot module genes with correlation-based edges
    plot(g1, layout=coords)
    title(sprintf("%s module", module_color))
}

#' Useful function for plotting a colorbar
#'
#' Source:
#' http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
plot_color_bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale <- (length(lut) - 1) / (max - min)

    dev.new(width=1.75, height=5)

    plot(c(0,10), c(min,max), 
         type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)

    for (i in 1:(length(lut)-1)) {
        y <- (i - 1) / scale + min
        rect(0, y, 10, y + 1 / scale, col=lut[i], border=NA)
    }
    dev.off()
}

#' Display sample PCA plot
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param dat Expression matrix with samples for columns
#' @param sample_ids Sample ids
#' @param condition factor of conditions
#' @param batch factor of batches
#' @param num_pcs Number of PCs to generate pairwise plots for (default=2)
#' @param main plot title
#' @param include_legend Whether or not to include a plot legend. 
#' @param scale Whether or not to perform a scaled PCA plot
#'
#' @return None
plot_sample_pca <- function(dat, sample_ids, condition, batch, num_pcs=2, 
                            main="", include_labels=TRUE, include_legend=TRUE, scale=FALSE) {
    # check to make sure request number of PC's to plot is valid
    if (num_pcs > ncol(dat)) {
        stop("Invalid number of PCs requested.")
    }

    # remove any zero-variance genes before performing PCA
    dat_filtered <- dat[apply(dat, 1, var) != 0,]

    prcomp_results <- prcomp(t(dat_filtered), scale=scale)

    # iterate over pairwise combinations of requested PCs
    pc_pairs <- combn(1:num_pcs, 2)

    for (i in 1:ncol(pc_pairs)) {
        pc1 <- pc_pairs[1,i]
        pc2 <- pc_pairs[2,i]

        # Percent variance explained
        var_explained <- round(summary(prcomp_results)$importance[2,] * 100, 2)

        xl <- sprintf("PC%d (%.2f%% variance)", pc1, var_explained[pc1])
        yl <- sprintf("PC%d (%.2f%% variance)", pc2, var_explained[pc2])

        # Dataframe for PCA plot
        df <- data.frame(sample_id=sample_ids,
                         pc1=prcomp_results$x[,pc1], pc2=prcomp_results$x[,pc2],
                         condition=condition, batch=batch)
        # PCA plot
        if (include_legend) {
          # with legend
          plt <- ggplot(df, aes(pc1, pc2, color=condition, shape=batch)) +
            theme_bw()
        } else {
          # without legend
          plt <- ggplot(df, aes(pc1, pc2)) +
            theme_bw() +
            theme(legend.position="none")
        }

        # scatter plot and axes labels
        plt <- plt + 
          geom_point(stat="identity",size=5) +
          xlab(xl) + ylab(yl) +
          ggtitle(sprintf("PCA: %s", main)) +
          theme(axis.ticks=element_blank(), axis.text.x=element_text(angle=-90))

        # text labels
        if (include_labels) {
          plt <- plt + geom_text(aes(label=sample_id), angle=45, size=4,vjust=2)
        }

        plot(plt)
    }
}


#' Display a sample heatmap
#'
#' Plots a sample similarity heatmap with sample condition plotted along the
#' x-axis and one or more sample covariates plotted along the y-axis. 
#' 
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts Count matrix
#' @param xlabels Factor to use for x-axis labels and colorbar (e.g. condition)
#' @param ylabels Factor to use for y-axis labels (e.g. batch)
#' @param covariates Dataframe containing one or more sample covariates
#' @param metric Distance metric to use when comparing samples
#' @param main Title for heatmap plot.
#'
#' @return None
plot_sample_heatmap <- function(counts, xlabels, ylabels, covariates=NULL, 
                                metric='dist', main='Sample Heatmap', ...) {
    # Convert to float to avoid issues with cor()
    counts <- counts * 1.0

    # Choose dissimilarity metric to use
    if (metric == 'dist') {
        heatmap_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100)
        dist_matrix <- as.matrix(dist(t(counts)))
    } else if (metric %in% c('pearson', 'spearman')) {
        heatmap_colors <- rev(RColorBrewer::brewer.pal(9, 'YlOrRd'))
        dist_matrix <- cor(counts, method=metric)
    }

    # discard any unused factor levels for labels
    xlabels <- factor(xlabels)
    ylabels <- factor(ylabels)

    # If no covariates supplied, just use ylabels
    if (is.null(covariates)) {
        covariates <- ylabels 
    }

    # If only a single covariate supplied as factor, convert to 1d data.frame
    if (!is.data.frame(covariates)) {
        covariates <- as.data.frame(covariates)
    }

    # Condition colormap
    pal <- colorRampPalette(rev(brewer.pal(n=9, name='Set1')))(length(unique(xlabels)))
    cond_colors = pal[as.integer(xlabels)]

    num_cond_cols <- max(9, max(as.integer(as.factor(xlabels))))
    cond_pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_cond_cols)
    cond_colors <- cond_pal[as.integer(xlabels)]

    # Create a matrix of covariate color assignments for use in heatmap plot
    adjust_var_colors <- RColorBrewer::brewer.pal(8, 'Set2')[as.factor(covariates[[names(covariates)[1]]])]

    # Labels
    colnames(dist_matrix) <- sprintf("%s (%s)", colnames(counts), xlabels)
    rownames(dist_matrix) <- ylabels

    # If multiple covariates are defined, use heatmap.plus to display colorbars
    # for each covariate
    if (ncol(covariates) > 1) {
        for (i in 2:ncol(covariates)) {
            cols <- RColorBrewer::brewer.pal(8, 'Set2')[as.factor(covariates[[names(covariates)[i]]])]
            adjust_var_colors <- cbind(adjust_var_colors, cols)
        }
        adjust_var_colors <- as.matrix(adjust_var_colors)
        colnames(adjust_var_colors) <- names(covariates)

        # heatmap.plus doesn't like 1d matrices for either RowSideColors or ColSideColors
        cond_colors_mat <- as.matrix(cbind(cond_colors, cond_colors))
        colnames(cond_colors_mat) <- c("Condition", "")

        heatmap.plus::heatmap.plus(dist_matrix,
                                  ColSideColors=cond_colors_mat,
                                  RowSideColors=adjust_var_colors,
                                  margin=c(13, 13), revC=FALSE,
                                  xlab='Condition', ylab='Covariates', main=main,
                                  col=heatmap_colors, ...)

    } else {
        # Otherwise, for single colorbars, we can just use heatmap.2
        gplots::heatmap.2(dist_matrix, 
                          ColSideColors=cond_colors,
                          RowSideColors=adjust_var_colors, 
                          margin=c(13, 13), revC=FALSE,
                          xlab='Condition', ylab='Covariates', main=main,
                          trace='none', col=heatmap_colors, ...)
    }
}


#' Plot median pairwise sample correlations
#'
#' Plots the median correlation of each sample with every other sample.
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param counts ExpressionSet count matrix
#' @param condition Vector of sample conditions
#' @param batch Vector of sample batches
#' @param main Plot title.
#' @param method Correlation method to use (default: pearson)
#' @param mar Vector of margin dimensions to use for plot.
#'
#' @return None
plot_sample_correlations <- function (counts, condition, batch, main="",
                                      method='pearson', mar=c(12,6,4,6)) {
    # Convert if type is ExpressionSet
    if (class(counts) == "ExpressionSet") {
        counts <- exprs(counts)
    }

    # Compute pairwise sample correlations
    median_pairwise_cor <- rowMedians(cor(counts * 1.0, method=method))

    quantiles <- quantile(median_pairwise_cor, probs=c(0.25, 0.75))
    iqr <- diff(quantiles)

    #outlimit
    cutoff <- quantiles[1] - 1.5 * iqr

    ylimit <- c(pmin(min(median_pairwise_cor), cutoff), max(median_pairwise_cor))

    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(max(as.integer(condition)))
    cond_colors <- pal[as.integer(condition)]

    # sample labels
    sample_labels <- sprintf("%s (%s / %s)", colnames(counts), condition, batch)

    # render plot
    par(mar=mar)
    plot(median_pairwise_cor, xaxt="n", ylim=ylimit,
         ylab="Median Pairwise Correlation", xlab="", main=main,
         col=cond_colors, pch=16, cex=2.2)
    axis(side=1, at=seq(along=median_pairwise_cor),
        labels=sample_labels, las=2)
    abline(h=cutoff, lty=2)
    abline(v=1:length(condition), lty=3, col="black")
}

