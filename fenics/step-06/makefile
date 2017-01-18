G=mesh

TARGETS = $(G).msh $(G).xml

all: $(TARGETS)

$(G).msh: $(G).geo
	gmsh -2 -optimize $(G).geo

$(G).xml: $(G).msh
	dolfin-convert -i gmsh -o xml $(G).msh $(G).xml

clean:
	rm -f $(TARGETS) *.pvd *.vtu
