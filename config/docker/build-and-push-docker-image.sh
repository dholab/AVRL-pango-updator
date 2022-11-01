#!/bin/bash

DOCKER_DESTINATION='nrminor/pango-extender:v1.0.1'

# build and push Docker image
docker buildx create --use
docker buildx build -t $DOCKER_DESTINATION --platform=linux/amd64,linux/arm64 . --push
