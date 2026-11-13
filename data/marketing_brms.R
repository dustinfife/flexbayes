# Keep fitted-model deserialization deferred in development and installations.
delayedAssign("marketing_brms", {
  utils::getFromNamespace("flexbayes_model_dataset", "flexbayes")("marketing_brms")
})
