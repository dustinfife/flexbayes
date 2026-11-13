# Keep fitted-model deserialization deferred in development and installations.
delayedAssign("marketing_stanreg", {
  utils::getFromNamespace("flexbayes_model_dataset", "flexbayes")("marketing_stanreg")
})
