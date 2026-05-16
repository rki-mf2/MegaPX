#include "InferenceParameterInspector.h"

#include <ATen/core/function_schema.h>

#include <algorithm>
#include <sstream>
#include <string>

namespace {

/*
* @fn argument_to_type
* @brief Converts a TorchScript schema argument type to printable text.
* @signature std::string argument_to_type(const c10::Argument& argument);
* @param argument: TorchScript schema argument or return value.
* @throws None.
* @return Type text, or "unknown" when no real type is available.
*/
std::string argument_to_type(const c10::Argument& argument) {
  if (!argument.real_type()) {
    return "unknown";
  }
  return argument.real_type()->repr_str();
}

/*
* @fn default_value_note
* @brief Formats a TorchScript schema argument default value.
* @signature std::string default_value_note(const c10::Argument& argument);
* @param argument: TorchScript schema argument.
* @throws None.
* @return Default-value note text, or an empty string when absent.
*/
std::string default_value_note(const c10::Argument& argument) {
  if (!argument.default_value()) {
    return {};
  }

  std::ostringstream out;
  out << "default=" << *argument.default_value();
  return out.str();
}

/*
* @fn contains_any
* @brief Tests whether any string contains any requested substring.
* @signature bool contains_any(const std::vector<std::string>& values, const std::vector<std::string>& needles);
* @param values: strings to search.
* @param needles: substrings to match.
* @throws None.
* @return True when at least one substring is present.
*/
bool contains_any(const std::vector<std::string>& values,
                  const std::vector<std::string>& needles) {
  return std::any_of(values.begin(), values.end(), [&](const auto& value) {
    return std::any_of(needles.begin(), needles.end(), [&](const auto& needle) {
      return value.find(needle) != std::string::npos;
    });
  });
}

/*
* @fn append_hint
* @brief Appends a model hint only if it has not already been added.
* @signature void append_hint(std::vector<std::string>& hints, std::string hint);
* @param hints: destination hint list.
* @param hint: hint text to append.
* @throws None.
* @return None.
*/
void append_hint(std::vector<std::string>& hints, std::string hint) {
  if (std::find(hints.begin(), hints.end(), hint) == hints.end()) {
    hints.push_back(std::move(hint));
  }
}

/*
* @fn format_parameter
* @brief Formats one inspected parameter entry for summary output.
* @signature std::string format_parameter(const InferenceParameterInspector::Parameter& item);
* @param item: parameter metadata to format.
* @throws None.
* @return Formatted parameter line.
*/
std::string format_parameter(const InferenceParameterInspector::Parameter& item) {
  std::ostringstream out;
  out << "  " << item.name << ": " << item.type;
  out << (item.required ? " required" : " optional");
  if (!item.source.empty()) {
    out << " [" << item.source << ']';
  }
  if (!item.note.empty()) {
    out << " - " << item.note;
  }
  return out.str();
}

}  // namespace

/*
* @fn InferenceParameterInspector
* @brief Constructs an inspector bound to a checkpoint loader.
* @signature InferenceParameterInspector::InferenceParameterInspector(const CheckpointModelLoader& loader);
* @param loader: checkpoint loader that owns model and archive metadata.
* @throws None.
* @return None.
*/
InferenceParameterInspector::InferenceParameterInspector(
    const CheckpointModelLoader& loader)
    : loader_(loader) {}

/*
* @fn inspect
* @brief Inspects available model metadata to infer forward inputs, outputs, and runtime requirements.
* @signature InferenceParameterInspector::Report InferenceParameterInspector::inspect() const;
* @throws std::logic_error when a TorchScript module is expected but unavailable.
* @return Structured inference-parameter report.
*/
InferenceParameterInspector::Report InferenceParameterInspector::inspect()
    const {
  Report report;

  if (loader_.has_torchscript_module()) {
    report.exact = true;
    report.reason = "TorchScript forward schema is available.";

    const auto& module = loader_.module();
    const auto method = module.get_method("forward");
    const auto& schema = method.function().getSchema();

    for (const auto& argument : schema.arguments()) {
      if (argument.name() == "self") {
        continue;
      }

      Parameter parameter;
      parameter.name = argument.name().empty() ? "<unnamed>" : argument.name();
      parameter.type = argument_to_type(argument);
      parameter.required = !argument.default_value();
      parameter.source = "forward schema";
      parameter.note = default_value_note(argument);
      if (argument.kwarg_only()) {
        if (!parameter.note.empty()) {
          parameter.note += "; ";
        }
        parameter.note += "keyword-only";
      }
      report.inputs.push_back(std::move(parameter));
    }

    int output_index = 0;
    for (const auto& result : schema.returns()) {
      Parameter output;
      output.name = result.name().empty()
                        ? "output_" + std::to_string(output_index)
                        : result.name();
      output.type = argument_to_type(result);
      output.required = true;
      output.source = "forward schema";
      report.outputs.push_back(std::move(output));
      ++output_index;
    }

    report.runtime_requirements.push_back(
        {"device", "cpu or cuda device", true, "runtime",
         "must match the device used to load input tensors"});
    report.runtime_requirements.push_back(
        {"dtype", "tensor scalar type", true, "runtime",
         "must match the model's expected input dtype"});
    return report;
  }

  if (loader_.is_python_checkpoint()) {
    report.exact = false;
    report.reason =
        "This is a Python/PyTorch-Lightning checkpoint, not a TorchScript "
        "module. The real forward inputs live in the original Python model "
        "class and cannot be recovered exactly from state_dict weights.";

    report.runtime_requirements.push_back(
        {"python_model_class", "source code", true, "checkpoint",
         "needed to instantiate the architecture before loading state_dict"});
    report.runtime_requirements.push_back(
        {"state_dict", "checkpoint tensor weights", true, "checkpoint",
         "present in archive/data.pkl"});
    report.runtime_requirements.push_back(
        {"preprocessing", "domain-specific transform", true, "user/model code",
         "must match training: scaling, tokenization, padding, masks, ordering"});
    report.runtime_requirements.push_back(
        {"exported_torchscript", ".pt file", true, "conversion step",
         "needed before C++ LibTorch can run inference directly"});

    const auto& names = loader_.sample_state_dict_names();
    if (contains_any(names, {"mz_encoder"})) {
      append_hint(report.model_hints,
                  "state_dict contains mz_encoder weights; inference likely "
                  "uses m/z peak-value tensors.");
    }
    if (contains_any(names, {"int_encoder"})) {
      append_hint(report.model_hints,
                  "state_dict contains int_encoder weights; inference likely "
                  "uses intensity tensors.");
    }
    if (contains_any(names, {"rt_encoder"})) {
      append_hint(report.model_hints,
                  "state_dict contains rt_encoder weights; inference likely "
                  "uses retention-time tensors.");
    }
    if (contains_any(names, {"level_encoder"})) {
      append_hint(report.model_hints,
                  "state_dict contains level_encoder weights; inference likely "
                  "uses MS level or category tensors.");
    }
    if (contains_any(names, {"transformer_encoder"})) {
      append_hint(report.model_hints,
                  "state_dict contains transformer_encoder weights; inference "
                  "likely needs padded sequences and an attention/padding mask.");
    }
    if (loader_.pytorch_lightning_version()) {
      append_hint(report.model_hints,
                  "checkpoint was saved by PyTorch-Lightning " +
                      *loader_.pytorch_lightning_version() + ".");
    }
    return report;
  }

  report.exact = false;
  report.reason =
      "No loadable TorchScript module or recognizable Python checkpoint was "
      "found.";
  report.runtime_requirements.push_back(
      {"model_file", "TorchScript .pt/.pth file", true, "user",
       "export the model with torch.jit.save or provide a compatible file"});
  return report;
}

/*
* @fn summary
* @brief Formats the inference-parameter inspection as human-readable text.
* @signature std::string InferenceParameterInspector::summary() const;
* @throws std::logic_error when a TorchScript module is expected but unavailable.
* @return Inspection summary text.
*/
std::string InferenceParameterInspector::summary() const {
  const auto report = inspect();
  std::ostringstream out;
  out << "Inference parameters: "
      << (report.exact ? "exact from model schema" : "best effort") << '\n';
  out << "Reason: " << report.reason << '\n';

  if (!report.inputs.empty()) {
    out << "Inputs:\n";
    for (const auto& input : report.inputs) {
      out << format_parameter(input) << '\n';
    }
  }

  if (!report.outputs.empty()) {
    out << "Outputs:\n";
    for (const auto& output : report.outputs) {
      out << format_parameter(output) << '\n';
    }
  }

  if (!report.runtime_requirements.empty()) {
    out << "Runtime requirements:\n";
    for (const auto& requirement : report.runtime_requirements) {
      out << format_parameter(requirement) << '\n';
    }
  }

  if (!report.model_hints.empty()) {
    out << "Model hints:\n";
    for (const auto& hint : report.model_hints) {
      out << "  " << hint << '\n';
    }
  }

  return out.str();
}
