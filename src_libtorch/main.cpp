#include "CheckpointModelLoader.h"
#include "InferenceParameterInspector.h"
#include <torch/torch.h>

#include <iostream>

int main(int argc, char** argv) {
  if (argc > 1) {
    CheckpointModelLoader loader;
    const bool loaded = loader.load(argv[1]);
    std::cout << loader.summary() << std::endl;
    InferenceParameterInspector inspector(loader);
    std::cout << inspector.summary() << std::endl;
    return loaded ? 0 : 1;
  }

  std::cout << "\nPass a checkpoint path to print a model summary and inference parameters:\n"
            << "  torch_loader src_libtorch/model/cascadia.ckpt\n";
  return 0;
}
