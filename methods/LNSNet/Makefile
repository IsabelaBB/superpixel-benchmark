
default: all

all: module

module: 
	pipenv install
	pipenv sync

clean:
	pipenv --rm
	
help:
	@echo make '{ module | clean }'

