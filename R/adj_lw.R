#' @title adj_lw_sp
#' @description
#' Create neighbours list with spatial weights for the chosen coding scheme for the given sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param neighbour_type default queen, may also be rook. \code{\link[spdep]{cell2nb}}
#' @param lw_style style can take values “W”, “B”, “C”, “U”, “minmax” and “S”. \code{\link[spdep]{nb2listw}}
#' @param combine default FALSE, if TRUE all individuals' neighbours list will be consolidated into one large list.
#' @return if combine=F, a value will be assigned to the list_neighbour slot corresponding to the sp object.
#' if combine=T, a value will be assigned to the combine_listw slot corresponding to the sp object.
#' @export
#' @rdname adj_lw_sp
adj_lw_sp <- function(object,neighbour_type = "queen",lw_style = "B",combine = F){

  individual <- unique(object@sample_features[["individual"]])
  sample_features <- object@sample_features
  sample_features <- sample_features[sample_features$`samples` %in% colnames(object@clean_data),]
  lws <- lapply(individual,adj_lw,location = sample_features,neighbour_type = neighbour_type)

  if(combine == F){

    result_vector <- map_dbl(lws, ~ .x[[2]])
    fi <- list()

    for (i in 1:length(lws)){
      my_list <- lws[[i]][[1]]
      new_list <- map(my_list, ~as.integer(.x))
      class(new_list) <- "nb"
      lw <- nb2listw(new_list, style=lw_style,zero.policy=TRUE)
      fi[[i]] <- lw
    }
    names(fi) <- individual
    object@list_neighbour <- fi

  }else{

    result_vector <- map_dbl(lws, ~ .x[[2]])
    fi <- list()
    for (i in 1:length(lws)) {
      if(i-1 == 0){
        add_value <- 0
      }else{
        add_value <- sum(result_vector[1:(i-1)])
      }
      my_list <- lws[[i]][[1]]
      new_list <- map(my_list, ~as.integer(.x + add_value))
      fi[[i]] <- new_list
    }

    combined_list <- do.call(c, fi)
    class(combined_list) <- "nb"
    lw <- nb2listw(combined_list, style=lw_style)
    object@combine_listw <- lw
  }
  return(object)
}

adj_lw <- function(i,location,neighbour_type = "queen"){

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
    samples = na.omit(as.vector(raster_matrix)),
    new_id1 = 1:n_valid
  )

  order_idx <- match(valid_indices$samples,id_mapping$samples)
  id_mapping <- id_mapping[order_idx, ]
  id_mapping <- mutate(id_mapping,new_id2 = 1:n_valid)


  nb_all <- cell2nb(nrow(raster_matrix), ncol(raster_matrix), type = neighbour_type,legacy = T)
  nb_valid <- lapply(valid_cell_numbers, function(id) {
    neighbors <- nb_all[[id]]
    valid_neighbors <- neighbors[neighbors %in% valid_cell_numbers]
    id_mapping$new_id2[match(valid_neighbors, id_mapping$original_id)]
  })

  nb_valid <- nb_valid[id_mapping$new_id1]
  nb_valid <- lapply(nb_valid, function(x) {
    if (is.null(x) || length(x) == 0) {
      0L
    } else {
      x
    }
  })
  class(nb_valid) <- "nb"
  z <- list(nb_valid,n_valid)

  return(z)
}
