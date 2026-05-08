#include "CascadiaModelConfig.h"
#include "CheckpointModelLoader.h"
#include "InferenceParameterInspector.h"

#include <torch/torch.h>

#include <filesystem>
#include <iostream>

int main(int argc, char** argv) {
  if (argc > 1) {
    const std::filesystem::path input_path(argv[1]);
    std::filesystem::path model_path = input_path;

    if (input_path.extension() == ".toml") {
      try {
        const auto config = CascadiaModelConfig::from_toml(input_path);
        std::cout << config.summary() << std::endl;
        model_path = config.model_path;
      } catch (const std::exception& error) {
        std::cerr << "Config error: " << error.what() << std::endl;
        return 1;
      }
    }

    CheckpointModelLoader loader;
    const bool loaded = loader.load(model_path);
    std::cout << loader.summary() << std::endl;
    InferenceParameterInspector inspector(loader);
    std::cout << inspector.summary() << std::endl;
    return loaded ? 0 : 1;
  }

  std::cout << "\nPass a checkpoint/model path or TOML config:\n"
            << "  torch_loader model/cascadia.ckpt\n"
            << "  torch_loader src_libtorch/cascadia_sequence.example.toml\n";
  return 0;
}
