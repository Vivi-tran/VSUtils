#!/bin/bash

flake8 . --count --max-complexity=13 --max-line-length=120 \
	--per-file-ignores="__init__.py:F401, filters.py:E402" \
	--exclude venv \
	--statistics