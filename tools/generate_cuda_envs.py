import os

# Mapping of CUDA versions to their maximum supported GCC versions
# This is based on https://gist.github.com/ax3l/9489132
CUDA_GCC_MAPPING = {
    "13.0": "15",
    "13.1": "15",
    
    "12.8": "14",
    "12.9": "14",

    "12.4": "13.2",
    "12.5": "13.2",
    "12.6": "13.2",

    "12.1": "12.2",
    "12.2": "12.2",
    "12.3": "12.2",

    "12.0": "12.1",

    "11.4.1": "11",
    "11.5": "11",
    "11.6": "11",
    "11.7": "11",
    "11.8": "11",

    "11.1": "10",
    "11.2": "10",
    "11.3": "10",
    "11.4.0": "10",

    "11.0": "9",

    "10.1": "8",
    "10.2": "8",

    "9.2": "7",
    "10.0": "7",

    "9.0": "6",
    "9.1": "6",
}

YAML_TEMPLATE = """name: cuda{cuda_version}
channels:
  - conda-forge
  - nvidia
dependencies:
  - gcc_linux-64={g_version}
  - gxx_linux-64={g_version}
  - cuda-toolkit={cuda_version}
  - cuda-version={cuda_version}
variables:
  OPENMM_CUDA_COMPILER: $CONDA_PREFIX/bin/nvcc
  CUDA_HOST_COMPILER: $CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
"""

def generate_cuda_yamls(target_directory: str) -> None:
    if not os.path.exists(target_directory):
        raise OSError(f"Target directory '{target_directory}' does not exist.")

    for cuda_v, gcc_v in CUDA_GCC_MAPPING.items():
        file_name = f"cuda{cuda_v}.yaml"
        file_path = os.path.join(target_directory, file_name)
        
        content = YAML_TEMPLATE.format(
            cuda_version=cuda_v,
            g_version=gcc_v
        )
        
        try:
            with open(file_path, "w") as f:
                f.write(content)
            print(f"Generated: {file_path}")
        except OSError as e:
            print(f"Error writing {file_name}: {e}")


if __name__ == "__main__":
    generate_cuda_yamls(target_directory="envs")