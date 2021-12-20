.PHONY: doc install all
MANDOC := man/*.Rd
RCODES := R/*.R

all: README.md doc

README.md: README.Rmd
	Rscript -e "devtools::build_readme()"

# cannot stop running
# $(MANDOC): $(RCODES)
	# Rscript -e "devtools::document()"

doc: $(RCODES)
	Rscript -e "devtools::document(pkg = '.')"

install:
	Rscript -e 'remotes::install_github("beyondpie/mssc")'
