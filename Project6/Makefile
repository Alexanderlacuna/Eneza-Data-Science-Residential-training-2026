# Eneza — convenience Makefile for Docker-based reproducible runs
#
# Available targets:
#   make build      - build the Docker image
#   make run        - run the full pipeline in a container
#   make shell      - open an interactive bash shell in the container
#   make freeze     - generate a fully frozen requirements-freeze.txt
#   make clean      - remove local output directories and Docker artefacts

IMAGE_NAME := eneza

.PHONY: build run shell freeze clean

build:
	docker build -t $(IMAGE_NAME) .

run:
	docker run --rm -v "$$(pwd)":/app $(IMAGE_NAME)

shell:
	docker run --rm -it -v "$$(pwd)":/app --entrypoint /bin/bash $(IMAGE_NAME)

freeze:
	docker run --rm --entrypoint /bin/bash -v "$$(pwd)":/app $(IMAGE_NAME) -c \
		"pip freeze > /app/requirements-freeze.txt"

clean:
	rm -rf output output2 outputs __pycache__ scripts/__pycache__ *.pyc
	docker rmi -f $(IMAGE_NAME) 2>/dev/null || true
