.PHONY:  test test-gui tags

all:

clean:
	rm -rf dist
	rm -rf doc/generated
	rm -rf __pycache__ */__pycache__

tags:
	etags *.py */*.py */*/*.py


test-load-all: clean
	find .. -name mfix.dat -print0 | MFIX_NO_VTK=1 xargs -0 -n 1 ./run.sh -linfo -t

thumbnails: clean
	find ../tutorials ../benchmarks -name mfix.dat -or -name *.mfx -print0 | xargs -0 -n 1 ./run.sh -linfo -t -ct

tutorial_info:
	cd ..; python -m mfixgui.tools.collect_tutorial_info

../mfix:
	cd ..; ./configure_mfix --python --smp
	make -C .. mfix mfixsolver.so

pyinst: all ../mfix
	pyinstaller gui.spec

sphinx-docs:
	cd doc; sphinx-apidoc -o ../generated/source ../..
	cp doc/index.rst doc/sphinx/conf.py doc/generated/source/
	cp doc/Makefile doc/generated/
	make -C doc/generated html

docs: sphinx-docs
