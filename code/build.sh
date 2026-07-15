#!/usr/bin/env bash
set -e

# 1. Help Menu Function
show_help() {
    echo ""
    echo "======================================================================"
    echo "cneta Build Script"
    echo "----------------------------------------------------------------------"
    echo "Usage: ./build.sh [target] [cores]"
    echo ""
    echo "Arguments:"
    echo "  target            'local' to build without loading modules,"
    echo "                    'clean' to remove the build directory,"
    echo "                    or a cluster name (e.g., Eureka2)."
    echo ""
    echo "                    Note: Cluster scripts must be located at:"
    echo "                    ./envs/setup_env_<target>.sh"
    echo ""
    echo "  cores             (Optional) Number of CPU cores to use."
    echo "                    Defaults to 1 if not specified."
    echo ""
    echo "Examples:"
    echo "  ./build.sh                # Shows this help menu"
    echo "  ./build.sh local          # Build locally with 1 core"
    echo "  ./build.sh local 4        # Build locally with 4 cores"
    echo "  ./build.sh Eureka2 8      # Build on Eureka2 with 8 cores"
    echo "  ./build.sh clean          # Delete the build directory"
    echo "======================================================================"
}

# 2. Check for empty args or help flags
if [ $# -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    show_help
    exit 0
fi

# 3. Check for too many arguments
if [ $# -gt 2 ]; then
    echo "Error: Too many arguments provided."
    echo ""
    show_help
    exit 1
fi

# 4. Handle 'clean' command
if [ "$1" == "clean" ]; then
    echo "Cleaning build directory..."
    rm -rf build/
    exit 0
fi

# 5. Target Selection Logic
if [ "$1" != "local" ]; then
    ENV_SCRIPT="envs/setup_env_${1}.sh"
    if [ -f "$ENV_SCRIPT" ]; then
        echo "Loading environment from $ENV_SCRIPT..."
        source "$ENV_SCRIPT"
    else
        echo "Error: Cannot find environment script '$ENV_SCRIPT'."
        echo "Did you spell the cluster name correctly? Or did you mean 'local'?"
        exit 1
    fi
else
    echo "Building locally (no HPC modules loaded)..."
fi

# 6. Validate Core Count (Graceful Failure)
CORES=${2:-1}
if ! [[ "$CORES" =~ ^[0-9]+$ ]]; then
    echo "Error: Core count must be a positive integer. You provided: '$CORES'"
    exit 1
fi

# 7. The Build Process
echo "Configuring cneta..."
mkdir -p build
cd build
cmake -Wno-dev ..

echo "Compiling with $CORES core(s)..."
make -j"$CORES"

echo "Build complete! Executables are in build/"
