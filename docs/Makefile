SPHINXBUILD ?= sphinx-build
SOURCEDIR    = .
BUILDDIR     = _build

.PHONY: html clean livehtml
html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html
livehtml:
	sphinx-autobuild $(SOURCEDIR) $(BUILDDIR)/html
clean:
	rm -rf $(BUILDDIR)
