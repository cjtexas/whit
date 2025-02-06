#' dif: A custom R6 class
#'
#' @section public methods:
#' \describe{
#'   \item{\code{new()}}{Creates a new object of class `dif`.}
#' }
#' @import MASS
#' @export
dif <- R6::R6Class("dif",
  public = list(
    ikm = NULL,
    ism = NULL,

    initialize = function(d) {
      self$ikm <- self$c_km(d)
      self$ism <- self$s_m()
    },

    c2a = function(c) {
      paste(gsub("[^ -~]", " ", c), collapse = " ")
    },

    c_km = function(d) {
      matrix(sample(1:95, d^2, replace = TRUE), nrow = d)
    },

    s_m = function() {
      oc <- 32:126
      sc <- sample(oc)
      setNames(sc, oc)
    },

    s_c = function(t, sm) {
      tn <- as.numeric(charToRaw(t))
      s <- sm[as.character(tn)]
      rawToChar(as.raw(s))
    },

    s_c_inv = function(t, sm) {
      tn <- as.numeric(charToRaw(t))
      inv <- as.numeric(names(sm))[match(tn, sm)]
      rawToChar(as.raw(inv))
    },

    bt_o = function(t) {
      tn <- as.numeric(charToRaw(t))
      o <- bitwXor(tn, 0xFF)
      rawToChar(as.raw(o))
    },

    tt_m = function(t) {
      tn <- as.numeric(charToRaw(t)) - 32
      d <- nrow(self$ikm)
      pl <- d^2 - length(tn)

      if (pl > 0) {
        tn <- c(tn, rep(0, pl))
      }

      matrix(tn, nrow = d)
    },

    ec = function(t) {
      st <- self$s_c(t, self$ism)
      ot <- self$bt_o(st)
      tmdt <- self$tt_m(ot)
      em <- as.matrix(self$ikm) %*% tmdt
      em
    },

    dc = function(em) {
      km_inv <- MASS::ginv(self$ikm)
      dm <- km_inv %*% em
      dm <- round(dm)
      dm <- dm[dm > 0]

      dt <- rawToChar(as.raw(dm + 32))
      ot <- rawToChar(as.raw(bitwXor(as.numeric(charToRaw(dt)), 0xFF)))
      ot <- self$s_c_inv(ot, self$ism)
      ot
    }
  )
)
