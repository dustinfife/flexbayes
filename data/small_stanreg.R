# Keep fitted-model deserialization deferred in development and installations.
delayedAssign("small_stanreg", {
  utils::getFromNamespace("flexbayes_model_dataset", "flexbayes")("small_stanreg")
})
