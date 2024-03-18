SHELL := /bin/bash

.PHONY: help clean purge run zip

help:
	@echo
	@echo "PLEASE USE CONDA ENV FOR RUNNING ANY OF THE MAKE TARGETS."
	@echo
	@echo "Usage: make <target>"
	@echo
	@echo "Available targets:"
	@echo
	@echo "    help - show this help message."
	@echo "    clean - remove results, logs and benchmarks."
	@echo "    purge - clean and remove zipped output, .snakemake and .parallel directory."
	@echo "    run - run the pipeline."
	@echo "    zip - create a zip file with results, logs and benchmarks."
	@echo

clean:
	@rm -rf results logs benchmarks

purge:
	@rm -rf results logs benchmarks .snakemake .parallel output.zip

run:
	@snakemake --profile config

zip:
	@zip -r output.zip results logs benchmarks
	@echo "Zip file created: output.zip"