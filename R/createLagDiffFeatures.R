#' @title Generate lags and differences for feature variables 2
#'
#' @description
#' Replace all variables with their generated lagged and differenced variables.
#' Uses the \code{xts} framework for developing lags and is only available for \code{TimeTasks}.
#'
#' @template arg_taskdf
#' @template arg_taskdf_target
#' @param lag [\code{integer}]\cr
#'   An integer vector of lag lengths.
#' @param difference [\code{integer}]\cr
#'   An integer of the order of differencing
#' @param cols [\code{character}]\cr
#'   A character vector of columns to create lag features for.
#'    Default is to use all columns. NOTE: For forecast regression tasks, it is not
#'    a good idea to make lags of your target variable. So if cols are not specied by
#'    the user, createLagDiffFeatures will return a regr task.
#' @param frequency [\code{integer}]\cr
#'   An integer representing the periodicity in the time series. If frequency is declared in the task,
#'   the task frequency will be used.
#' @param na.pad [\code{logical}]\cr
#'   A logical to denote whether the data should be padded to the original size with NAs
#' @param difference.lag [\code{integer}]\cr
#'   An integer denoting the period to difference over
#' @param seasonal.difference.lag [\code{integer}]\cr
#'   An integer denoting the period to seasonaly difference over
#' @param return.nonlag [\code{logical}]\cr
#'   A logical to denote whether the original unlagged features should be returned
#' @export
#' @family eda_and_preprocess
#' @examples
#' set.seed(1234)
#' dat = arima.sim(model = list(ar = c(.5,.2), ma = c(.4), order = c(2,0,1)), n = 200)
#' times = (as.POSIXlt("1992-01-14")) + lubridate::days(1:200)
#' dat = xts::xts(dat,order.by = times)
#' colnames(dat) = c("arma_test")
#' regr.task = makeRegrTask(id = "Lagged ML model", data = as.data.frame(dat), target = "arma_test")
#' regr.task.lag = createLagDiffFeatures(regr.task, lag = 1L:10L, difference = 0L)

createLagDiffFeatures = function(obj, lag = 0L, difference = 0L,
                                   cols = NULL, target = character(0L), frequency = 1L,
                                   na.pad = TRUE, return.nonlag = TRUE, stratify = NULL) {
  ## FIXME: differences only accepts one value, should allow for more
  assertIntegerish(lag,lower = 0L, upper = 1000L)
  assertIntegerish(difference,lower = 0L, upper = 1000L, len = 1L)
  assertLogical(na.pad)
  assert(checkClass(obj, "data.frame"), checkClass(obj, "Task"))
  assertCharacter(target, any.missing = FALSE)
  if (!is.null(cols))
    assertCharacter(cols, any.missing = FALSE)

  UseMethod("createLagDiffFeatures")
}

#' @export
#' @importFrom zoo cbind.zoo
#' @import data.table
createLagDiffFeatures.data.frame = function(obj, lag = 0L, difference = 0L,
                                              cols = NULL, target = character(0L), frequency = 1L,
                                              na.pad = TRUE, return.nonlag = TRUE, stratify = NULL) {
  # Columns to work with
  nums = sapply(obj,is.numeric)
  work.cols = colnames(obj)
  if (!is.null(cols)) {
    assertSubset(cols, work.cols)
    x = obj[,cols]
  } else {
    cols = work.cols
  }

  # Make seasonal lags
  seasonal.lag        = lag * frequency
  seasonal.difference = difference * frequency

  # FIXME: We have to use data.table, otherwise the code is disgusting
  x <- data.table::data.table(obj)

  # Make names of lags and seasonal lags
  if (length(cols) > 1){
    lag_column_names <- paste0(t(replicate(max(lag),cols)),"_lag", lag)
    seasonal_lag_column_names <- unique(paste0("seasonal_",t(replicate(max(seasonal.lag),cols)),"_lag", seasonal.lag))
  } else {
    lag_column_names <- paste0(replicate(max(lag),cols),"_lag", lag)
    seasonal_lag_column_names <- unique(paste0("seasonal_",replicate(max(seasonal.lag),cols),"_lag", seasonal.difference))
  }

  # Names of differences
  diff_cols = cols[nums]
  if (length(diff_cols > 1)){
    diff_column_names <- paste0(diff_cols,"_diff", difference)
    seasonal_diff_column_names <- unique(paste0("seasonal_",t(replicate(difference,diff_cols)),"_diff", difference))
  } else {
    diff_column_names <- paste0(replicate(difference,diff_cols),"_diff", difference)
    seasonal_diff_column_names <- unique(paste0("seasonal_",replicate(seasonal.difference,diff_cols),"_diff", seasonal.difference))
  }

  # NOTE: The if statement here with get() allows for stratify lags across groups
  if (lag > 0)
    x[,(lag_column_names) := shift(.SD, n = lag), by = if (!is.null(stratify)){
      get(stratify)
    } else {
      NULL
    }, .SDcols = cols]

  if (difference > 0)
    x[,(diff_column_names) := lapply(.SD, function(x) c(rep(NA,difference),
                                                        diff(x,difference ))), by = if (!is.null(stratify)){
                                                          get(stratify)
                                                        } else {
                                                          NULL
                                                        }, .SDcols = diff_cols]

  if (frequency > 1){
    if (lag > 0){
      x[,(seasonal_lag_column_names) := shift(.SD, seasonal.lag), by = if (!is.null(stratify)){
        get(stratify)
      } else {
        NULL
      }, .SDcols = cols]
    }
    if (difference > 0){
      x[,(seasonal_diff_column_names) := lapply(.SD, function(x) c(rep(NA,seasonal.difference),
                                                                   diff(x, lag = seasonal.difference ))), by = if (!is.null(stratify)){
                                                                     get(stratify)
                                                                   } else {
                                                                     NULL
                                                                   }, .SDcols = diff_cols]
    }
  }

  if (return.nonlag){
    obj = obj[,c(intersect(work.cols, cols)), drop = FALSE]
    obj = cbind(obj, data.frame(x)[,setdiff(colnames(x), colnames(obj))])
  } else {
    obj = cbind(obj[,target],x)
  }

  if (na.pad == FALSE){
    removeNaPad = 1:(max(lag, seasonal.lag) +
                       (max(difference,
                            seasonal.difference)))
    obj = obj[-removeNaPad,]
  }
  return(obj)
}

#' @export
createLagDiffFeatures.Task = function(obj, lag = 0L, difference = 0L, cols = NULL,
                                        target = character(0L),
                                        frequency = 1L, na.pad = TRUE,
                                        return.nonlag = TRUE, stratify = NULL) {
  target = getTaskTargetNames(obj)
  data = getTaskData(obj)


  if (!is.null(obj$task.desc$frequency) && frequency == 1L)
    frequency = obj$task.desc$frequency
  if (!is.null(cols)){
    if (!all(cols %in% colnames(data)))
      stop("Chosen cols not in data")
    work.cols = intersect(colnames(data),cols)
  } else {
    work.cols = colnames(data)
  }
  # We store the original columns as we need them for forecasting
  data.original = data


  data = data[,work.cols, drop = FALSE]
  row.dates = try(as.POSIXct(rownames(data)), silent = TRUE)
  if (is.error(row.dates))
    stop("The data's rownames must be convertible to POSIXct")

  data = suppressWarnings(createLagDiffFeatures( obj = data, lag = lag, difference = difference,
                                                   cols = cols, target = target,
                                                   frequency = frequency, na.pad = na.pad,
                                                   return.nonlag = return.nonlag, stratify = stratify))


  obj = changeData(obj,data = data)
  data.original = data.original[I(nrow(data) - (max(lag) * frequency + difference * frequency)):I(nrow(data)),,drop=FALSE]
  obj$task.desc$pre.proc$data.original = data.original
  obj$task.desc$pre.proc$par.vals = list(lag = lag, difference = difference,
                                         cols = cols, target = target,
                                         frequency = frequency, na.pad = na.pad,
                                         return.nonlag = return.nonlag, stratify = stratify)
  obj
}


#' @title Update Lagged and Differenced Data
#'
#' @description
#' Update a task modified by lags and differences with new data
#'
#' @param task [\code{\link{Task}}]\cr
#'   The task.
#' @param newdata [\code{data.frame}]\cr
#'   New observations to update the model
#' @param weights [\code{numeric}]\cr
#'   Optional, non-negative case weight vector to be used during fitting.
#'   If given, must be of same length as \code{subset} and in corresponding order.
#'   By default \code{NULL} which means no weights are used unless specified in the task (\code{\link{Task}}).
#'   Weights from the task will be overwritten.
#' @param ... [any]\cr
#'   Currently ignored.
#' @return [\code{\link{WrappedModel}}].
#' @export
#' @importFrom xts rbind.xts
updateLagDiff = function(task, newdata, weights,...) {
  assertClass(task, "Task")
  assertDataFrame(newdata)
  if (missing(weights))
    weights = task$weights



  data = rbind(task$task.desc$pre.proc$data.original, newdata)
  row.dates = try(as.POSIXct(rownames(data)), silent = TRUE)
  if (is.error(row.dates))
    stop("The data's rownames must be convertible to POSIXct")

  data = xts::xts(data, order.by = row.dates)
  lagdiff.func = function(...){
    createLagDiffFeatures(obj = data,...)
  }
  data = do.call(lagdiff.func, task$task.desc$pre.proc$par.vals)
  data = data.frame(data, row.names = index(data))

  obj = changeData(task, data, weights)
  obj
}

#' @export
as.data.frame.xts = function(x, ...){
  data.frame(row.names = index(x), zoo::coredata(x))
}

