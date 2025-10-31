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
                       
                       print = function(show = NULL, ...) {
                         m <- as.matrix(self$distance)
                         cat("MDist object\n")
                         cat("  Preset :", self$preset, "\n")
                         cat("  Size   :", paste(dim(m), collapse = " x "), "\n")
                         
                         # Determine which parameters to show
                         keys <- show
                         if (is.null(keys)) keys <- private$param_keys_by_preset[[self$preset]]
                         if (is.null(keys)) keys <- names(self$params)
                         
                         if (!is.null(self$params) && length(self$params) && length(keys)) {
                           cat("  Parameters:\n")
                           for (k in keys) {
                             if (k %in% names(self$params)) {
                               label <- private$pretty_labels[[k]] %||% k  # fallback to key name
                               val   <- private$format_val(self$params[[k]])
                               cat("    - ", label, ": ", val, "\n", sep = "")
                             }
                           }
                         }
                         
                         # Optional small preview of the matrix
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
                       
                       # which parameters to display for each preset
                       param_keys_by_preset = list(
                         gower              = c("commensurable"),
                         custom             = c("distance_cont", "distance_cat", "scaling_cont", "commensurable", "ncomp", "threshold"),
                         euclidean_onehot   = c("scaling_cont"),
                         unbiased_dependent = c("distance_cont", "distance_cat", "scaling_cont", "commensurable"),
                         dkss               = character(0),
                         gudmm              = character(0),
                         mod_gower          = character(0)
                       ),
                       
                       # mapping of internal parameter names to pretty labels
                       pretty_labels = list(
                         distance_cont  = "Continuous distance method",
                         distance_cat   = "Categorical distance method",
                         scaling_cont   = "Scaling for continuous vars",
                         commensurable  = "Commensurability adjustment",
                         ncomp          = "Number of continuous variables",
                         threshold      = "Variance threshold",
                         n_cont         = "Number of continuous variables",
                         n_cat          = "Number of categorical variables"
                       ),
                       
                       # helper for formatting values cleanly
                       format_val = function(v) {
                         if (is.null(v)) return("NULL")
                         if (is.logical(v)) return(ifelse(v, "TRUE", "FALSE"))
                         if (length(v) > 1) return(paste(v, collapse = ", "))
                         as.character(v)
                       }
                     )
)
