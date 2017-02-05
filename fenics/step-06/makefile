G=mesh

TARGETS = msh xml

all: $(TARGETS)

msh: $(G).geo
	gmsh -2 -optimize $(G).geo

xml: $(G).msh
	dolfin-convert -i gmsh -o xml $(G).msh $(G).xml

clean:
	rm -f $(TARGETS) *.pvd *.vtu
