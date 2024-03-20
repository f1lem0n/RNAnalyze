SHELL := /bin/bash

.PHONY: help blastdb clean dryrun purge run zip

help:
	@echo
	@echo "PLEASE USE CONDA ENV FOR RUNNING ANY OF THE MAKE TARGETS."
	@echo
	@echo "Usage: make <target>"
	@echo
	@echo "Available targets:"
	@echo
	@echo "    help - show this help message."
	@echo "    blastdb - download blast database."
	@echo "    clean - remove results, logs and benchmarks."
	@echo "    dryrun - run the pipeline in dryrun mode."
	@echo "    purge - clean and remove zipped output, .snakemake and .parallel directory."
	@echo "    run - run the pipeline."
	@echo "    zip - create a zip file with results, logs and benchmarks."
	@echo

blastdb:
	@for db in $$(awk '{if($$1 == "database:") print $$2}' config/params.yaml); do \
		db=$$(echo $$db | sed 's/\"//g'); \
		if [ -d database/$$db ]; then \
			if [ ! -z "$$(ls -A database/$$db)" ]; then \
				echo "Database $$db already exists."; \
				continue; \
			fi; \
		fi; \
		echo "Searching for $$db..."; \
		perl database/update_blastdb.pl $$db; \
		echo "Moving and unpacking..."; \
		mkdir database/$$db; \
		mv $$db*tar.gz* database/$$db/.; \
		cd database/$$db && \
		tar -xvzf $$db.tar.gz && \
		rm -rf $$db.tar.gz && \
		cd - > /dev/null; \
		echo "Success!"; \
	done

clean:
	@rm -rf results logs benchmarks

dryrun:
	@snakemake --profile config --dryrun

purge:
	@rm -rf results logs benchmarks .snakemake .parallel output.zip

run:
	@snakemake --profile config

zip:
	@zip -r output.zip results logs benchmarks
	@echo "Zip file created: output.zip"