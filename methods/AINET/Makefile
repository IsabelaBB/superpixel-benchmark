
default: all

all: module

module: 
	pipenv install
	pipenv sync
	cd third_party/cython/; pipenv run python -m pip install .; cd ../../; exit

clean:
	pipenv --rm
	
help:
	@echo make '{ module | clean }'

