context("learners_all_classif2")

detachAllPackages()
R.utils::gcDLLs()

test_that("learners work: classif2 ", {
  detachAllPackages()
  R.utils::gcDLLs()
  # settings to make learners faster and deal with small data size
  hyperpars = list(
    classif.boosting = list(mfinal = 2L),
    classif.cforest = list(mtry = 1L),
    classif.bartMachine = list(verbose = FALSE, run_in_sample = FALSE,
                               # without this (and despite use_missing_data being TRUE), the test with missing data fails with a null point exception, which manifests itself as a completely different rJava error in the test
                               replace_missing_data_with_x_j_bar = TRUE,
                               num_iterations_after_burn_in = 10L),
    classif.bdk = list(ydim = 2L),
    classif.earth = list(degree = 3L, nprune = 2L),
    classif.gbm = list(bag.fraction = 1, n.minobsinnode = 1),
    classif.lssvm = list(kernel = "rbfdot", reduced = FALSE),
    classif.nodeHarvest = list(nodes = 100L, nodesize = 5L),
    classif.xyf = list(ydim = 2L),
    classif.h2o.deeplearning = list(hidden = 2L, seed = getOption("mlr.debug.seed"), reproducible = TRUE),
    classif.h2o.randomForest = list(seed = getOption("mlr.debug.seed"))
  )


  # binary classif
  task = subsetTask(binaryclass.task, subset = c(10:20, 180:190),
                    features = getTaskFeatureNames(binaryclass.task)[12:15])

  # binary classif with prob
  lrns = mylist(binaryclass.task, properties = "prob", create = TRUE)
  lapply(lrns, testBasicLearnerProperties, task = binaryclass.task,
         hyperpars = hyperpars, pred.type = "prob")

})

test_that("learners work: classif2 ", {
  detachAllPackages()
  R.utils::gcDLLs()
  # settings to make learners faster and deal with small data size
  hyperpars = list(
    classif.boosting = list(mfinal = 2L),
    classif.cforest = list(mtry = 1L),
    classif.bartMachine = list(verbose = FALSE, run_in_sample = FALSE,
                               # without this (and despite use_missing_data being TRUE), the test with missing data fails with a null point exception, which manifests itself as a completely different rJava error in the test
                               replace_missing_data_with_x_j_bar = TRUE,
                               num_iterations_after_burn_in = 10L),
    classif.bdk = list(ydim = 2L),
    classif.earth = list(degree = 3L, nprune = 2L),
    classif.gbm = list(bag.fraction = 1, n.minobsinnode = 1),
    classif.lssvm = list(kernel = "rbfdot", reduced = FALSE),
    classif.nodeHarvest = list(nodes = 100L, nodesize = 5L),
    classif.xyf = list(ydim = 2L),
    classif.h2o.deeplearning = list(hidden = 2L, seed = getOption("mlr.debug.seed"), reproducible = TRUE),
    classif.h2o.randomForest = list(seed = getOption("mlr.debug.seed"))
  )


  # binary classif
  task = subsetTask(binaryclass.task, subset = c(10:20, 180:190),
                    features = getTaskFeatureNames(binaryclass.task)[12:15])
  detachAllPackages()
  R.utils::gcDLLs()
  # binary classif with weights
  lrns = mylist("classif", properties = "weights", create = TRUE)
  lapply(lrns, testThatLearnerRespectsWeights, hyperpars = hyperpars,
         task = binaryclass.task, train.inds = binaryclass.train.inds, test.inds = binaryclass.test.inds,
         weights = rep(c(10000L, 1L), c(10L, length(binaryclass.train.inds) - 10L)),
         pred.type = "prob", get.pred.fun = getPredictionProbabilities)



  detachAllPackages()
  R.utils::gcDLLs()
  # classif with missing
  lrns = mylist("classif", properties = "missings", create = TRUE)
  lapply(lrns, testThatLearnerHandlesMissings, task = task, hyperpars = hyperpars)

  detachAllPackages()
  # classif with variable importance
  lrns = mylist("classif", properties = "featimp", create = TRUE)
  lapply(lrns, testThatLearnerCanCalculateImportance, task = task, hyperpars = hyperpars)
  R.utils::gcDLLs()
})

