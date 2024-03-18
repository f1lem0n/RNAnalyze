SHELL := /bin/bash

.PHONY: help clean purge zip

help:
	@echo
	@echo "MAKE TARGETS ARE FOR DEVELOPMENT PURPOSES ONLY."
	@echo "PLEASE USE CONDA ENV FOR RUNNING ANY OF THE MAKE TARGETS."
	@echo
	@echo "Available targets:"
	@echo "help - show this help message."
	@echo "clean - remove results, logs and benchmarks."
	@echo "purge - clean and remove zipped output and .snakemake directory."
	@echo "zip - create a zip file with results, logs and benchmarks."

clean:
	@rm -rf results logs benchmarks

purge:
	@rm -rf results logs benchmarks .snakemake output.zip

zip:
	@zip -r output.zip results logs benchmarks
	@echo "Zip file created: output.zip"