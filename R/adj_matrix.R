#' @title adj_matrix_sp
#' @description
#' Create neighbours list with spatial weights for the chosen coding scheme for the given sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param neighbour_type default queen, may also be rook. \code{\link[spdep]{cell2nb}}
#' @param code_style style can take values W, B, C, and S. \code{\link[spdep]{nb2mat}}
#' @return a value will be assigned to the adj_matrix_list slot corresponding to the sp object.
#' @export
#' @rdname adj_matrix_sp
adj_matrix_sp <- function(object,neighbour_type = "queen",code_style="W"){
  individual <- unique(object@sample_features[["individual"]])
  sample_features <- object@sample_features
  sample_features <- sample_features[sample_features$`samples` %in% colnames(object@clean_data),]
  adj_matrix <- lapply(individual,adj,location = sample_features,neighbour_type = neighbour_type,code_style=code_style)
  names(adj_matrix) <- individual
  object@adj_matrix_list <- adj_matrix
  object
}

adj <- function(i,location,neighbour_type = "queen",code_style="W"){

  valid_indices <- location[location$`individual` == i,]
  valid_indices <- valid_indices[,c("x","y","samples")]
  valid_indices$y <- (max(valid_indices$y)+1) - valid_indices$y

  matrix_df <- valid_indices %>%
    pivot_wider(
      id_cols = y,
      names_from = x,
      values_from = samples  # 用sample填充值
    ) %>%
    arrange(y)

  raster_matrix <- as.matrix(matrix_df[,-1])
  rownames(raster_matrix) <- matrix_df$y

  sorted_idx <- order(as.numeric(colnames(raster_matrix)))

  raster_matrix<- raster_matrix[, sorted_idx]

  colnames(raster_matrix) <- 1:ncol(raster_matrix)

  valid_cells <- !is.na(raster_matrix)

  n_valid <- sum(valid_cells)

  number_cells <- 1:length(raster_matrix)
  dim(number_cells) <- dim(raster_matrix)
  valid_cell_numbers <- number_cells[valid_cells]
  id_mapping <- data.frame(
    original_id = valid_cell_numbers,
    new_id = 1:n_valid
  )

  nb_all <- cell2nb(nrow(raster_matrix), ncol(raster_matrix), type = neighbour_type,legacy = T)
  nb_valid <- lapply(valid_cell_numbers, function(id) {
    neighbors <- nb_all[[id]]
    valid_neighbors <- neighbors[neighbors %in% valid_cell_numbers]
    id_mapping$new_id[match(valid_neighbors, id_mapping$original_id)]
  })

  nb_valid <- lapply(nb_valid, function(x) {
    if (is.null(x) || length(x) == 0) {
      0L
    } else {
      x
    }
  })

  class(nb_valid) <- "nb"
  adj_matrix <- nb2mat(nb_valid, style = code_style, zero.policy = T)
  k <- na.omit(as.vector(raster_matrix))
  dimnames(adj_matrix) <- list(k,k)
  adj_matrix <- adj_matrix[valid_indices[["samples"]],valid_indices[["samples"]]]

  return(adj_matrix)
}
