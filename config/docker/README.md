This directory contains the code needed to build and push the Docker image for this workflow. This code was written to ensure that the workflow will run on ARM-64 systems (e.g., Apple Silicon Macs) or on AMD-64 systems.

Most users will not need to re-build and push the image, but for the sake of documentation, we did so by changing into this directory and running:

```
bash build-and-push-docker-image.sh
```
