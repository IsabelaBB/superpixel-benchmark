
default: all

all: module

module: 
	pipenv install
	pipenv sync
	#pipenv run pip install numpy pillow scikit-image pysnic

clean:
	pipenv --rm
	
help:
	@echo make '{ module | clean }'

