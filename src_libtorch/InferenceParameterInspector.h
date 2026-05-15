#pragma once

#include "CheckpointModelLoader.h"

#include <string>
#include <vector>

class InferenceParameterInspector {
 public:
  struct Parameter {
    std::string name;
    std::string type;
    bool required = true;
    std::string source;
    std::string note;
  };

  struct Report {
    bool exact = false;
    std::string reason;
    std::vector<Parameter> inputs;
    std::vector<Parameter> outputs;
    std::vector<Parameter> runtime_requirements;
    std::vector<std::string> model_hints;
  };

/*
* @fn InferenceParameterInspector
* @brief Constructs an inspector bound to a checkpoint loader.
* @signature explicit InferenceParameterInspector(const CheckpointModelLoader& loader);
* @param loader: checkpoint loader that owns model and archive metadata.
* @throws None.
* @return None.
*/
  explicit InferenceParameterInspector(const CheckpointModelLoader& loader);

/*
* @fn inspect
* @brief Inspects available model metadata to infer forward inputs, outputs, and runtime requirements.
* @signature Report inspect() const;
* @throws std::logic_error when a TorchScript module is expected but unavailable.
* @return Structured inference-parameter report.
*/
  Report inspect() const;

/*
* @fn summary
* @brief Formats the inference-parameter inspection as human-readable text.
* @signature std::string summary() const;
* @throws std::logic_error when a TorchScript module is expected but unavailable.
* @return Inspection summary text.
*/
  std::string summary() const;

 private:
  const CheckpointModelLoader& loader_;
};
