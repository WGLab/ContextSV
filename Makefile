# Top-Level Makefile

.PHONY: python cpp clean

# Targets for the sub-makefiles
python:
	$(MAKE) -f Makefile-python

cpp:
	$(MAKE) -f Makefile-cpp

debug:
	$(MAKE) -f Makefile-cpp DEBUG=1

clean:
	$(MAKE) -f Makefile-python clean
	$(MAKE) -f Makefile-cpp clean
