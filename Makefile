# Top-Level Makefile

.PHONY: python cpp clean

# Targets for the sub-makefiles
python:
	$(MAKE) -f Makefile-python

cpp:
	$(MAKE) -f Makefile-cpp

clean:
	$(MAKE) -f Makefile-python clean
	$(MAKE) -f Makefile-cpp clean
