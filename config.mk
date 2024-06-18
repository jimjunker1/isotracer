### Configuration file for Makefile ###

TOP_DIR=$(shell git rev-parse --show-toplevel)
MAKEFLAGS += --no-print-directory

### * Colors

RED = "\\033[31m"
BLUE = "\\033[94m"
GREEN = "\\033[92m"
NC = "\\033[0m"

### * help (default rule)

.PHONY: help
help: Makefile
	@printf "\n"
	@printf "Please use 'make <target>' where <target> is one of\n"
	@printf "\n"
	@sed -n 's/^## /    /p' $< | column -t -s ":"
