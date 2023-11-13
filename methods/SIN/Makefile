
default: all

all: module

module: 
	pipenv install
	# pipenv run pip install numpy cython scikit-image opencv-python matplotlib torch torchvision
	pipenv sync

clean:
	pipenv --rm
	
help:
	@echo make '{ module | clean }'

