#include "CheckpointModelLoader.h"
#include <torch/torch.h>

#include <iostream>

int main(int argc, char** argv) {
  if (argc > 1) {
    CheckpointModelLoader loader;
    const bool loaded = loader.load(argv[1]);
    std::cout << loader.summary() << std::endl;
    return loaded ? 0 : 1;
  }

  const torch::Tensor tensor = torch::rand({2, 3});
  std::cout << tensor << std::endl;
  std::cout << "\nPass a checkpoint path to print a model summary:\n"
            << "  torch_loader src_libtorch/model/cascadia.ckpt\n";
  return 0;
}
