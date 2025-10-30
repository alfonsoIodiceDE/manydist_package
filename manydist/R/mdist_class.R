#' MDist: R6 class for mixed-type distance objects
#'
#' Stores the distance matrix, preset information, and computation parameters.
#' @export
MDist <- R6::R6Class("MDist",
                     public = list(
                       data     = NULL,
                       distance = NULL,
                       preset   = NULL,
                       params   = NULL,
                       
                       initialize = function(distance, preset, data = NULL, params = list()) {
                         self$data     <- data
                         self$distance <- distance
                         self$preset   <- preset
                         self$params   <- params
                       },
                       
                       # Optional: edit the per-preset view at runtime
                       set_param_view = function(preset, keys) {
                         stopifnot(is.character(preset), is.character(keys))
                         private$param_keys_by_preset[[preset]] <- keys
                         invisible(self)
                       },
                       
                       # Print with preset-specific parameter selection
                       print = function(show = NULL, ...) {
                         m <- as.matrix(self$distance)
                         cat("MDist object\n")
                         cat("  Preset :", self$preset, "\n")
                         cat("  Size   :", paste(dim(m), collapse = " x "), "\n")
                         
                         # Determine which params to show
                         keys <- show
                         if (is.null(keys)) keys <- private$param_keys_by_preset[[self$preset]]
                         if (is.null(keys)) keys <- names(self$params)  # fallback: show all
                         
                         if (!is.null(self$params) && length(self$params) && length(keys)) {
                           cat("  Parameters:\n")
                           for (k in keys) {
                             if (k %in% names(self$params)) {
                               cat("    - ", k, ": ", private$format_val(self$params[[k]]), "\n", sep = "")
                             }
                           }
                         }
                         
                         # Small preview (optional)
                         if (nrow(m) > 0 && ncol(m) > 0) {
                           cat("\n  Preview of distance matrix:\n")
                           print(round(m[1:min(3, nrow(m)), 1:min(3, ncol(m))], 3))
                         }
                         invisible(self)
                       },
                       
                       to_dist = function() {
                         if (inherits(self$distance, "dist")) return(self$distance)
                         as.dist(self$distance)
                       }
                     ),
                     
                     private = list(
                       # Which params to display for each preset (edit as you like)
                       param_keys_by_preset = list(
                         gower              = c("commensurable"),
                         custom             = c("distance_cont", "distance_cat", "scaling_cont", "commensurable", "ncomp", "threshold"),
                         euclidean_onehot   = c("scaling_cont"),
                         unbiased_dependent = c("distance_cont", "distance_cat", "scaling_cont", "commensurable"),
                         dkss               = character(0),  # or set your own
                         gudmm              = character(0),
                         mod_gower          = character(0)
                       ),
                       
                       format_val = function(v) {
                         if (is.null(v)) return("NULL")
                         if (is.logical(v)) return(ifelse(v, "TRUE", "FALSE"))
                         if (length(v) > 1) return(paste(v, collapse = ", "))
                         as.character(v)
                       }
                     )
)
