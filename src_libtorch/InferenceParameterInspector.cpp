#include "InferenceParameterInspector.h"

#include <ATen/core/function_schema.h>

#include <algorithm>
#include <sstream>
#include <string>

namespace {

std::string argument_to_type(const c10::Argument& argument) {
  if (!argument.real_type()) {
    return "unknown";
  }
  return argument.real_type()->repr_str();
}

std::string default_value_note(const c10::Argument& argument) {
  if (!argument.default_value()) {
    return {};
  }

  std::ostringstream out;
  out << "default=" << *argument.default_value();
  return out.str();
}

bool contains_any(const std::vector<std::string>& values,
                  const std::vector<std::string>& needles) {
  return std::any_of(values.begin(), values.end(), [&](const auto& value) {
    return std::any_of(needles.begin(), needles.end(), [&](const auto& needle) {
      return value.find(needle) != std::string::npos;
    });
  });
}

void append_hint(std::vector<std::string>& hints, std::string hint) {
  if (std::find(hints.begin(), hints.end(), hint) == hints.end()) {
    hints.push_back(std::move(hint));
  }
}

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

InferenceParameterInspector::InferenceParameterInspector(
    const CheckpointModelLoader& loader)
    : loader_(loader) {}

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
