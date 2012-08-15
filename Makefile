BRANCH_NAME = $(shell git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1 /')

docs:
	cd docs;make dirhtml
	git checkout gh-pages
	cp -r docs/_build/dirhtml/* .
	git add .
	git commit -am "Updating docs"
	git push -u origin gh-pages
	git checkout ${BRANCH_NAME}
	rm -rf _themes

.PHONY: docs
