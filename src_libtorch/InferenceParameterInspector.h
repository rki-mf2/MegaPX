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

  explicit InferenceParameterInspector(const CheckpointModelLoader& loader);

  Report inspect() const;
  std::string summary() const;

 private:
  const CheckpointModelLoader& loader_;
};
