docs:
	branch_name="$(git symbolic-ref HEAD 2>/dev/null)"
	cd docs;make dirhtml
	git checkout gh-pages
	cp -r docs/_build/dirhtml/* .
	git add .
	git commit -am "Updating docs"
	git push -u origin gh-pages
	git checkout $(branch_name)
	rm -rf _themes

.PHONY: docs
