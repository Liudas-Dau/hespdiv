#' Plot hespdiv dendrogram
#'
#' @description This function displays a dendrogram of polygons produced by hespdiv split-lines. Branch length is proportional to difference. If performance of split-lines is a similarity measure, it is internally converted to difference.
#' @param obj A hespdiv object.
#' @param poly.scheme ggplot2 object produced with poly_scheme function. Provide if you want identical colors for polygons in both plots.
#' @param color color vector used for dendrogram nodes and branches.
#' @param performance.col color vector used for text, displaying difference values between polygons
#' @param labels.col color vector used for text, displaying polygon IDs.
#' @param offset.factor numeric value used to scale the offset distance of displayed polygon IDs and performance values from a dendrogram node. Adjust experimentally, if you don't like the current distance.
#' @param arrange logical. Plot hespdiv dendrogram above polygon scheme?
#' @param grob logical. Convert plot to grob? Must be true, if you want to arrange polygon scheme and the dendrogram in a single plot.
#' @return A `grob` or `TableGrob` object if `grob = TRUE`, otherwise `NULL`.
#' @importFrom igraph graph_from_data_frame V E ends vcount layout_as_tree
#' @importFrom gridExtra grid.arrange
#' @importFrom grid grid.grabExpr grid.draw grid.newpage
#' @importFrom gridGraphics grid.echo
#' @author Liudas Daumantas
#' @note If you want to transform similarity to difference externally, before applying dendro, change maximize to TRUE in the call info of obj.
#' @family HespDiv visualization options
#' @examples
#' ## DO NOT RUN:
#' ## display dendrogram
#' # dendro(HDData::hd)
#' ## use the same colors as in poly_scheme
#' # scheme <- poly_scheme(HDData::hd, seed = 2)
#' # scheme
#' # dendro(HDData::hd, scheme, arrange = FALSE)
#' ## When arranging plots, you may need to experimentally adjust offset.factor.
#' # dendro(HDData::hd, scheme, offset.factor = 2)
#' ## End(Not run)
#' @export
dendro <- function(obj, poly.scheme = NULL, color = 1, performance.col = "blue",
                   labels.col = 1, offset.factor = 1, arrange = TRUE,
                   grob = TRUE){
  if (arrange & !grob) {
    stop("Can only arrange polygon scheme and dendrogram in a single plot,
         when dendrogram is grob.")
  }

  if (grob) {

    # Wrapper function to call .base_dendro
    .call_base <- function() {
      .base_dendro(obj, poly.scheme = poly.scheme, color = color,
                   performance.col = performance.col, labels.col = labels.col,
                   offset.factor = offset.factor)
    }

    grob <- grid::grid.grabExpr(gridGraphics::grid.echo(
      .call_base))

    if (arrange & !is.null(poly.scheme)){
      gridExtra::grid.arrange(grob, poly.scheme, ncol = 1)
    } else {
      grid::grid.newpage()
      grid::grid.draw(grob)
      grob
    }
  } else {
    .base_dendro(obj, poly.scheme = poly.scheme, color = color,
                 performance.col = performance.col, labels.col = labels.col,
                 offset.factor = offset.factor)
  }

}

# generate the dendrogram using base R graphics
#' @noRd
.base_dendro <- function(obj, poly.scheme = NULL, color = 1,
                         performance.col = "blue",
                         labels.col = 1, offset.factor = 1 ){
  # extract essential data
  pols <- obj$poly.stats[, c("plot.id", "root.id", "has.split")]
  pols$performance <- ifelse(
    test = obj$poly.stats$has.split,
    yes = obj$split.stats$performance[match(obj$poly.stats$plot.id,
                                            obj$split.stats$plot.id)],
    no = NA
  )

  if (!is.null(poly.scheme)) {
    pols$color <- poly.scheme$layers[[4]]$aes_params$colour
  } else {
    if (is.null(color)){
      stop("Either poly.scheme or color must be provided")
    }
    pols$color <- color
  }


  # Step 1: Calculate branch lengths for child nodes
  if (!obj$call.info$Call_ARGS$maximize) {
    if (all(pols$performance >= 0, na.rm = TRUE) &
        all(pols$performance <= 1, na.rm = TRUE)) {
      # assuming that maximum theoretical value is 1
      pols$branch_length <- 1 - pols$performance
    } else {
      if (any(pols$performance < 0, na.rm = TRUE)) {
        # removing negative values
        pols$performance <- pols$performance - min(pols$performance)
      }
      # Assuming that maximum similarity is the one recorded
      pols$branch_length <- max(pols$performance) - pols$performance
    }
  } else {
    if (any(pols$performance < 0, na.rm = TRUE)) {
      # removing negative values
      pols$performance <- pols$performance - min(pols$performance)
    }
    pols$branch_length <- pols$performance
  }

  # Step 2: Create edges with branch lengths
  edges <- data.frame(
    from = pols$root.id[pols$root.id != 0],  # Parent cluster
    to = pols$plot.id[pols$root.id != 0],    # Child cluster
    weight = pols$branch_length[match(pols$root.id[pols$root.id != 0],
                                      pols$plot.id)]  # Branch length
  )

  # Step 3: Build the graph with weights and vertex attributes
  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = pols)

  # Add `plot_id` and `color` attributes to vertices
  igraph::V(g)$plot_id <- pols$plot.id
  igraph::V(g)$color <- pols$color

  # Assign edge colors based on the child node's color
  vertex_map <- stats::setNames(seq_along(igraph::V(g)), names(igraph::V(g)))

  igraph::E(g)$color <- igraph::V(g)$color[
    vertex_map[as.character(igraph::ends(g,igraph::E(g))[, 2])]]

  # Step 4: Calculate node positions manually
  node_positions <- numeric(length(vertex_map))  # Initialize positions for all nodes
  for (edge in igraph::E(g)) {
    parent <- vertex_map[as.character(igraph::ends(g, edge)[1])]  # Map parent to index
    child <- vertex_map[as.character(igraph::ends(g, edge)[2])]  # Map child to index
    node_positions[child] <- node_positions[parent] + igraph::E(g)$weight[edge]
  }

  # Step 5: Create a layout reflecting the computed node positions
  layout <- igraph::layout_as_tree(g)  # Use tree layout for x-coordinates
  layout[, 2] <- -node_positions  # Set y-coordinates to the computed positions (invert to go downward)

  # Step 6: Prepare data for lateral-first edges
  edge_segments <- data.frame(
    x0 = layout[vertex_map[as.character(igraph::ends(g, igraph::E(g))[, 1])], 1],  # Parent x
    y0 = layout[vertex_map[as.character(igraph::ends(g, igraph::E(g))[, 1])], 2],  # Parent y
    x1 = layout[vertex_map[as.character(igraph::ends(g, igraph::E(g))[, 2])], 1],  # Child x
    y1 = layout[vertex_map[as.character(igraph::ends(g, igraph::E(g))[, 2])], 2],  # Child y
    weight = igraph::E(g)$weight  # Edge weight for segment length
  )
  # Step 7: Plot the dendrogram with lateral-first edges
  plot(
    NULL,
    xlim = range(layout[, 1])*1.1,  # Set x-axis limits
    ylim = range(layout[, 2])*1.1,  # Set y-axis limits
    xlab = "",          # Label for x-axis
    ylab = "",            # Label for y-axis
    main = "",
    xaxt = "n",                 # Suppress default x-axis
    yaxt = "n",                       # Suppress default y-axis
    bty = "n"                         # Remove the box around the plot
  )

  # Draw lateral-first edges
  for (i in 1:nrow(edge_segments)) {
    # Horizontal edge from parent to the x-position of the child
    lines(
      x = c(edge_segments$x0[i], edge_segments$x1[i]),
      y = c(edge_segments$y0[i], edge_segments$y0[i]),
      col = igraph::E(g)$color[i]
    )
    # Vertical edge from horizontal line down to the child node
    lines(
      x = c(edge_segments$x1[i], edge_segments$x1[i]),
      y = c(edge_segments$y0[i], edge_segments$y1[i]),
      col = igraph::E(g)$color[i]
    )
    # Get parent x- and y-coordinates
    parent_x <- edge_segments$x0[i]
    parent_y <- edge_segments$y0[i]
    # Add segment length only for left branches
    if (edge_segments$x1[i] < edge_segments$x0[i]) {
      graphics::text(
        x = parent_x,
        y = parent_y - 0.2 * offset.factor * stats::dist(range(layout[, 2]) *
                                                           1.1) /
          length(unique(obj$split.stats$rank)),  # Slightly below the parent node
        labels = round(edge_segments$weight[i], 2),
        col = performance.col,
        cex = 0.8
      )
    }
  }

  # Add nodes (circles)
  points(
    layout[, 1],
    layout[, 2],
    pch = 21,
    bg = igraph::V(g)$color,
    cex = 2
  )

  # Add cluster labels below the nodes
  graphics::text(
    layout[, 1],
    layout[, 2] - 0.1 * offset.factor * stats::dist(range(layout[, 2]) * 1.1) /
      length(unique(obj$split.stats$rank)),  # Slightly below the nodes
    labels = igraph::V(g)$plot_id,
    col = labels.col,
    pos = NULL,  # Place them directly below
    cex = 0.8
  )

  rounded_y_ticks <- round(seq(min(layout[, 2]), max(layout[, 2]), length.out = 5), 2)

  graphics::axis(
    side = 2,                                # Left-side Y-axis
    at = rounded_y_ticks,                    # Rounded tick positions
    labels = -rounded_y_ticks,                # Rounded tick labels
    las = 1,                                 # Rotate labels to be horizontal
    line = -0.5                              # Adjust spacing of axis to prevent overlap
  )
}
