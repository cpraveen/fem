CXX = g++ -Wall

TARGETS = sparse_test linsol_test

all: $(TARGETS)

sparse_matrix.o: sparse_matrix.cc sparse_matrix.h
Vector.o: Vector.cc Vector.h
cg_solver.o: cg_solver.cc cg_solver.h
jacobi_solver.o: jacobi_solver.cc jacobi_solver.h
sor_solver.o: sor_solver.cc sor_solver.h
ssor_solver.o: ssor_solver.cc ssor_solver.h

sparse_test: sparse_matrix.o Vector.o sparse_test.o 
	$(CXX) -o $@ $^

linsol_test: sparse_matrix.o Vector.o cg_solver.o jacobi_solver.o \
             sor_solver.o ssor_solver.o linsol_test.o
	$(CXX) -o $@ $^

clean:
	rm -f $(TARGETS) *.o
