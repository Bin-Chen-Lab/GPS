# RCL Training Docker Image

This Docker image runs the RCL traning model script for GPS model creation.  
It comes with all dependencies pre-installed and organized project directories.

---

## 🚀 Quick Start

mkdir output

Run the container with this command:

docker run --rm --gpus all \
    -v $(pwd)/output:/app/code/results \
    leshchi4/rcl:latest \
    python code/main_multi.py --cl VCAP_t1 --num_networks 4 --forget_rate 0.1 --stop_fuzzy 2 --fuzzy_exponent 3 --seed 6





