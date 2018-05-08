.PHONY: all lint test test-cov install dev

lint:
	flake8 q2_american_gut setup.py

test: all
	py.test

test-cov: all
	py.test --cov=q2_american_gut

install: all
	python setup.py install

dev: all
	pip install -e .
