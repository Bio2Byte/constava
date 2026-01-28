PACKAGE_NAME = constava

.PHONY: build install uninstall publish clean

build:
# 	python3 setup.py sdist bdist_wheel
	python3 -m build

install:
	python3 -m pip install dist/$(PACKAGE_NAME)-*.tar.gz

.PHONY: test
test:
	constava test

uninstall:
	-pip uninstall -y $(PACKAGE_NAME)

publish:
	twine upload dist/*

clean:
	rm -rf dist build $(PACKAGE_NAME).egg-info
